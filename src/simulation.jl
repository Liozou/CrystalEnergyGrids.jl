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
Base.@kwdef struct SimulationSetup{Trecord}
    temperatures::Vector{TK}
    ncycles::Int
    ninit::Int=0
    outdir::String=""
    printevery::Int=1000
    record::Trecord=Returns(nothing)

    function SimulationSetup(T, ncycles::Int, ninit::Int=0, outdir::String="", printevery::Int=1000, record=Returns(nothing))
        @assert ncycles ≥ 0 && ninit ≥ 0
        n = ncycles + ninit
        temperatures = if T isa Number
            fill(T, n + (n == 0))
        elseif T isa Vector
            convert(Vector{TK}, T)
        else
            n ≤ 1 ? [T(1,2)] : T.(1:n, n)
        end
        ret = new{typeof(record)}(temperatures, ncycles, ninit, outdir, printevery, record)
        initialize_record!(record, ret)
        ret
    end
end
function SimulationSetup(T, ncycles; kwargs...)
    SimulationSetup(; temperatures=T, ncycles, kwargs...)
end

initialize_record!(::T, ::SimulationSetup{T}) where {T} = nothing
needcomplete(::Any) = true
needcomplete(::Returns) = false

mutable struct MoveKind
    trials::Int
    accepted::Int
end
MoveKind() = MoveKind(0, 0)
attempt!(x::MoveKind) = (x.trials += 1)
accept!(x::MoveKind) = (x.accepted += 1)
(x::MoveKind)() = x.accepted / x.trials

mutable struct MoveStatistics
    dmax::TÅ # maximum distance for a translation move
    θmax::typeof(1.0u"°") # maximum angle for a rotation move
    const translation::MoveKind
    const rotation::MoveKind
    const randomtranslation::MoveKind
end
function MoveStatistics(dmax, θmax)
    MoveStatistics(dmax, θmax, MoveKind(), MoveKind(), MoveKind())
end
function accept!(statistics, move)
    # avoid using a type-unstable getfield here by explicitly splitting
    if move === :translation
        accept!(statistics.translation)
    elseif move === :rotation
        accept!(statistics.rotation)
    end
    nothing
end


"""
    choose_newpos!(statistics::MoveStatistics, pos::Vector{SVector{3,TÅ}}, randomdmax::TÅ, bead)

Choose an MC move for the molecule at positions given by `pos` and return a pair
`(newpos, movekind)` where `newpos` is the list of atom positions after the move, and
`movekind` is the kind of move that was chosen, as a `Symbol`.

The list of possible MC moves is:
- :translation : displacement of the molecule while conserving the angle.
- :rotation : rotation of the molecule around the atom given by `bead`.
- :randomtranslation : displacement of the molecule while conserving the angle up to a
  distance given by `randomdmax` (which is the maximum distance between two points in the
  unit cell).
"""
function choose_newpos!(statistics::MoveStatistics, pos::AbstractVector{SVector{3,TÅ}}, randomdmax::TÅ, bead)
    r = rand()
    if length(pos) == 1
        if r < 0.98
            attempt!(statistics.translation)
            random_translation(pos, statistics.dmax)
        else
            random_translation(pos, randomdmax)
        end, :translation
    else
        if r < 0.49
            attempt!(statistics.translation)
            random_translation(pos, statistics.dmax), :translation
        elseif r < 0.98
            attempt!(statistics.rotation)
            random_rotation(pos, statistics.θmax, bead), :rotation
        else
            random_rotation(random_translation(pos, randomdmax), 90u"°", bead), :randomtranslation
        end
    end
end


"""
    run_montecarlo!(mc::MonteCarloSetup, simu::SimulationSetup)

Run a Monte-Carlo simulation. Return the list of energies of the system during the
simulation.

See [`MonteCarloSetup`](@ref) for the definition of the system and
[`SimulationSetup`](@ref) for the parameters of the simulation.
"""
function run_montecarlo!(mc::MonteCarloSetup, simu::SimulationSetup)
    # report
    energy = baseline_energy(mc)
    reports = [energy]

    # record and outputs
    mkpath(simu.outdir)
    output, output_task = pdb_output_handler(isempty(simu.outdir) ? "" : joinpath(simu.outdir, "trajectory.pdb"), mc.step.mat)
    startstep = SimulationStep(mc.step, needcomplete(simu.record) ? :all : :output)
    if mc.step.parallel
        record_task = @spawn $simu.record(startstep, $energy, -simu.ninit, $mc, $simu)
    else
        simu.record(startstep, energy, -simu.ninit, mc, simu) === :stop && return reports
        record_task = @spawn nothing # for type stability
    end

    if simu.ncycles == 0 && simu.ninit == 0 # single-point computation
        put!(output, SimulationStep(mc.step, :zero))
        wait(record_task)
        wait(output_task[])
        return reports
    end

    # idle initialisations
    running_update = @spawn nothing
    local oldpos::Vector{SVector{3,TÅ}}
    local before::MCEnergyReport
    local after::MCEnergyReport

    # value initialisations
    nummol = max(20, length(mc.indices))
    old_idx = (0,0)
    accepted = false
    statistics = MoveStatistics(1.3u"Å", 30.0u"°")
    randomdmax = maximum(norm, eachcol(mc.step.mat))

    # main loop
    for (counter_cycle, idx_cycle) in enumerate((-simu.ninit+1):simu.ncycles)
        temperature = simu.temperatures[counter_cycle]

        for idnummol in 1:nummol
            # choose the species on which to attempt a move
            idx = choose_random_species(mc)
            ffidxi = mc.step.ffidx[idx[1]]

            # currentposition is the position of that species
            # currentposition is either a Vector or a @view, that's OK
            currentposition = (accepted&(old_idx==idx)) ? oldpos : @view mc.step.positions[mc.step.posidx[idx[1]][idx[2]]]

            # newpos is the position after the trial move
            newpos, move = choose_newpos!(statistics, currentposition, randomdmax, mc.bead[idx[1]])

            # check blocking pockets to refuse the trial early when possible
            if inblockpocket(mc.speciesblocks[idx[1]], mc.atomblocks, ffidxi, newpos)
                # note that if the move is blocked here, neither oldpos, old_idx nor
                # accepted are updated.
                # This is allows waiting for the next cycle before waiting on the
                # running_update task.
                # However it will still be taken it into account for the translation and
                # rotation trial values updates.
                continue
            end

            # if the previous move was accepted, wait for mc to be completely up to date
            # before computing
            accepted && wait(running_update)

            i, j = idx
            ij = mc.offsets[i]+j

            if old_idx == idx
                if accepted
                    before = after
                # else
                #     before = fetch(before_task) # simply keep its previous value
                end
                @spawnif mc.step.parallel begin
                    singlereciprocal = @spawn single_contribution_ewald($mc.ewald, $ij, $newpos)
                    fer = @spawn framework_interactions($mc, $i, $newpos)
                    singlevdw = single_contribution_vdw(mc.step, idx, newpos)
                    after = MCEnergyReport(fetch(fer)::FrameworkEnergyReport, singlevdw, fetch(singlereciprocal)::TK)
                end
            else
                molpos = mc.step.posidx[i][j]
                positions = @view mc.step.positions[molpos]
                @spawnif mc.step.parallel begin
                    singlevdw_before_task = @spawn single_contribution_vdw($mc.step, $(i,j), $positions)
                    singlereciprocal_after = @spawn single_contribution_ewald($mc.ewald, $ij, $newpos)
                    singlereciprocal_before = @spawn single_contribution_ewald($mc.ewald, $ij, nothing)
                    fer_before = @spawn framework_interactions($mc, $i, $positions)
                    fer_after = @spawn framework_interactions($mc, $i, $newpos)
                    singlevdw_before = fetch(singlevdw_before_task)::TK
                    singlevdw_after = single_contribution_vdw(mc.step, (i,j), newpos)
                    before = MCEnergyReport(fetch(fer_before), singlevdw_before, fetch(singlereciprocal_before))
                    after = MCEnergyReport(fetch(fer_after), singlevdw_after, fetch(singlereciprocal_after))
                end
            end

            old_idx = idx
            accepted = compute_accept_move(before, after, temperature)
            oldpos = newpos
            if accepted
                if mc.step.parallel
                    # do not use newpos since it can be changed in the next iteration before the Task is run
                    running_update = @spawn update_mc!($mc, $idx, $oldpos)
                else
                    update_mc!(mc, idx, oldpos)
                end
                diff = after - before
                accept!(statistics, move)
                if abs(Float64(diff.framework)) > 1e50 # an atom left a blocked pocket
                    wait(running_update)
                    energy = baseline_energy(mc) # to avoid underflows
                else
                    energy += diff
                end
            end
        end

        # end of cycle
        report_now = idx_cycle ≥ 0 && (idx_cycle == 0 || (simu.printevery > 0 && idx_cycle%simu.printevery == 0))
        if !(simu.record isa Returns) || report_now
            accepted && wait(running_update)
            o = SimulationStep(mc.step, :output)
            if report_now
                put!(output, o)
                push!(reports, energy)
            end
            if !(simu.record isa Returns)
                ocomplete = needcomplete(simu.record) ? SimulationStep(o, :complete_output) : o
                if mc.step.parallel
                    # ret_record cannot be inferred, but that's fine
                    ret_record = fetch(record_task)
                    ret_record === :stop && break
                    record_task = @spawn $simu.record($ocomplete, $energy, $idx_cycle, $mc, $simu)
                else
                    simu.record(ocomplete, energy, idx_cycle, mc, simu)
                end
            end
            yield()
        end

        translation_ratio = statistics.translation()
        if !isnan(translation_ratio)
            statistics.dmax *= 1 + (translation_ratio - 0.5)*(1.0 + 3*translation_ratio)*100/(99+counter_cycle)
        end
        rotation_ratio = statistics.rotation()
        if !isnan(rotation_ratio)
            statistics.θmax = (counter_cycle-1)/counter_cycle*statistics.θmax + rotation_ratio/counter_cycle*120u"°"
        end
    end
    wait(record_task)
    accepted && wait(running_update)
    push!(reports, energy)
    put!(output, SimulationStep(mc.step, :zero))
    isempty(simu.outdir) || serialize(joinpath(simu.outdir, "energies.serial"), reports)
    wait(output_task[])
    fetch(simu.record)
    reports
end
