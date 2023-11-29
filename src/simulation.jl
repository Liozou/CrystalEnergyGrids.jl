export SimulationSetup

"""
    SimulationSetup{Trecord}

Definition of the simulation setup, which must provide the following information:
- `T` the temperature. Either a number (in K), a list of temperatures (in K, one per cycle)
  or a function. If a function, it must take two parameters: `j` the number of the current
  cycle and `n` the total number of cycles, and return a temperature (in K).
  Internally, the list of temperatures is stored as a field `temperatures` indexed by `j`,
  which varies between `1` and `n`.
- `ncycles` the number of cycles of the simulation. Each cycle consists in `N` steps, where
  a step is a Monte-Carlo movement attempt on a random species, and `N` is the number of
  species.
- `ninit` the number of initialization cycles, where no output is written and no energy
  stored. The total number of cycles `n` used in the definition of `T` corresponds to
  `ncycles + ninit`. The function `record` defined below is still called during the
  initialization cycles.
- `outdir` the path to a directory in which the result is stored. If not provided (or set
  to an empty string), no record will be stored on file.
- `printevery`: one record of the positions is stored at the end of every `printevery`
  cycles. Default to 1000. The first (at the end of initialization) and last positions
  are always recorded: use a value of 0 in printevery to only record those.
- `record::Trecord` a function which can be used to record extra information at the end of
  each cycle. `record` can return `:stop` to end the computation at that point, otherwise
  its return value is ignored. Stopping may not be immediate and can take up to one
  additional cycle. `record` is called even during the initialization cycles.

If provided, `record` must have the signature
    `record(o::SimulationStep, e::BaselineEnergyReport, k::Int, mc::MonteCarloSetup, simu::SimulationSetup)`
where
- `o` contains the information on the positions of the species at this point. Read below
  for more information about what can be done with `o` depending on `needcomplete`
- `e` is the energy at this point.
- `k` is the cycle number (between `-ninit` and `ncycles`). Negative values of `k`
  correspond to the initialization cycles. `k == -ninit` corresponds to the starting
  configuration, before the initialization cycles. `k == 0` corresponds to the state before
  the start of the production cycles. `k == ncycles` is the last state of the simulation.
  !!! tip
      The cycle number `k` corresponds to the temperature at index `j = k + ninit`.
- `mc` is the input `MonteCarloSetup`. Note that the `positions` field may be incorrect,
  use that of `o`. It must not be modified by the `record` function.
- `simu` is this `SimulationSetup`. It must not be modified by the `record` function.

By default, the `record` function is given as an argument a value of `o` which is a copy
of `mc.step` which is thread-safe to read and modify (except for the `ff`, `charges`,
`isrigid` and `ffidx` fields which are never meant to change).
If the `record` function is guaranteed to only ever read or modify the `positions` and
`atoms` fields of `o`, then it can specify it by defining a method for
`needcomplete(::Trecord)` that returns `false` (the default is `true`).

!!! warning
    If the value returned by `needcomplete` is `false` when it should be `true`, this may
    lead to race conditions and related subtle bugs, which may not be easy to identify.

Additionally, if a call to `record` can spawn asynchronous calls, a method for `Base.fetch`
should be implemented such that `Base.fetch(record)` blocks until all asynchronous calls
are done.

Upon creation of a `SimulationSetup`, a call to `initialize_record!(record, setup)` will be
realized. By default, `initialize_record!` does nothing, but you can specialize it for
for your appropriate record type if need be, i.e. by defining a method of signature:
`initialize_record!(record::T, setup::SimulationSetup{T}) where {T<:MyRecordType}`.

## Example

This is an example of a record function that stores the averages of the distance between
the particles and a given point, weighted by the current temperature:

```julia
struct AverageDistanceTemperatureRecord <: Function
    refpoint::SVector{3,typeof(1.0u"Å")}
    weightedaverages::Vector{typeof(1.0u"Å * K")}
end
AverageDistanceTemperatureRecord(refpoint) = AverageDistanceTemperatureRecord(refpoint, [])

function (x::AverageDistanceTemperatureRecord)(o, _, k, _, simu)
    k ≤ 0 && return # skip initialization cycles
    buffer, ortho, safemin = prepare_periodic_distance_computations(o.cell)
    weightedaverage = 0.0u"Å * K"
    for pos in o.positions
        buffer .= pos .- x.refpoint
        weightedaverage += periodic_distance!(buffer, o.cell, ortho, safemin)
    end
    x.weightedaverages[k] = weightedaverage * simu.temperatures[k] / length(o.atoms)
    nothing
end

# Initialize the internal `weightedaverages` vector to the right size
function CrystalEnergyGrids.initialize_record!(record::T, setup::SimulationSetup{T}) where {T<:AverageDistanceTemperatureRecord}
    resize!(record.weightedaverages, setup.ncycles)
end

# No need to give a complete copy of `o` since we only read the positions
CrystalEnergyGrids.needcomplete(::AverageDistanceTemperatureRecord) = false

# (optional since nothing is blocking here)
Base.fetch(::AverageDistanceTemperatureRecord) = nothing
```

To perform a simulation that runs for 1000 cycles with logarithmically increasing
temperatures, no trajectory stored but a record the average weighted distance to point
`(0.5, 0.9, 1.3)` (in Å) and energies stored every 2 cycles, use
`SimulationSetup((j,n)->log1p(j/n)*300u"K"/log(2), 1000, "", 2, AverageDistanceTemperatureRecord([0.5, 0.9, 1.3]u"Å"))`
"""
struct SimulationSetup{Trecord}
    ncycles::Int
    record::Trecord

    function SimulationSetup(; ncycles::Int, record=Returns(nothing))
        ret = new{typeof(record)}(ncycles, record)
        initialize_record!(record, ret)
        ret
    end
end

initialize_record!(::T, ::SimulationSetup{T}) where {T} = nothing



"""
    run_montecarlo!(mc::MonteCarloSetup, simu::SimulationSetup)

Run a Monte-Carlo simulation. Return the list of energies of the system during the
simulation.

See [`MonteCarloSetup`](@ref) for the definition of the system and
[`SimulationSetup`](@ref) for the parameters of the simulation.
"""
function run_montecarlo!(mc::MonteCarloSetup, simu::SimulationSetup)
    for idx_cycle in 1:10

        for idnummol in 1:mc.indices
            # choose the species on which to attempt a move
            idx = rand(mc.rng, 1:mc.indices)

            currentposition = mc.step.positions[idx:idx]

            mc.step.positions[idx:idx] .= random_translation(mc.rng, currentposition, 1.2)
        end
        # end of cycle
        put!(simu.record.lb, deepcopy(mc.step))
        yield()
    end
    fetch(simu.record)
    nothing
end
