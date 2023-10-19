export SimulationSetup

"""
    SimulationSetup{Trecord}

Definition of the simulation setup, which must provide the following information:
- `T` the temperature. Either a number (in K), a list of temperatures (in K, one per cycle)
  or a function. If a function, it must take two parameters: `k` the number of the current
  cycle and `n` the total number of cycles, and return a temperature (in K).
  Internally, the list of temperatures is stored as a field `temperatures` indexed by `k`.
- `ncycles` the number of cycles of the simulation. Each cycle consists in `N` steps, where
  a step is a Monte-Carlo movement attempt on a random species, and `N` is the number of
  species.
- `outdir` the path to a directory in which the result is stored. If not provided (or set
  to an empty string), no record will be stored on file.
- `printevery`: one record of the positions is stored at the end of every `printevery`
  cycles. Default to 1000. The initial and last positions are always recorded: use a value
  of 0 in printevery to only record those.
- `record::Trecord` a function which can be used to record extra information at the end of
  each cycle. `record` can return `:stop` to end the computation at that point, otherwise
  its return value is ignored. Stopping may not be immediate and can take up to one
  additional cycle.

If provided, `record` must have the signature
    `record(o::SimulationStep, e::BaselineEnergyReport, k::Int, mc::MonteCarloSetup, simu::SimulationSetup)`
where
- `o` contains the information on the positions of the species at this point. Read below
  for more information about what can be done with `o` depending on `needcomplete`
- `e` is the energy at this point.
- `k` is the cycle number (between 0 and `cycles`). `k == 0` corresponds to the
  initialization of the simulation.
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
`SimulationSetup((k,n)->log1p(k/n)*300u"K"/log(2), 1000, "", 2, AverageDistanceTemperatureRecord([0.5, 0.9, 1.3]u"Å"))`
"""
Base.@kwdef struct SimulationSetup{Trecord}
    temperatures::Vector{TK}
    ncycles::Int
    outdir::String=""
    printevery::Int=1000
    record::Trecord=Returns(nothing)

    function SimulationSetup(T, ncycles::Int, outdir::String="", printevery::Int=1000, record=Returns(nothing))
        temperatures = if T isa Number
            fill(T, ncycles + (ncycles == 0))
        elseif T isa Vector
            convert(Vector{TK}, T)
        else
            ncycles ≤ 1 ? [T(1,2)] : T.(1:ncycles, ncycles)
        end
        ret = new{typeof(record)}(temperatures, ncycles, outdir, printevery, record)
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

    waittime = parse(Float64, ARGS[1])

    startstep = SimulationStep(mc.step, :output)
    simu.record(deepcopy(startstep), energy, 0, deepcopy(mc), simu)


    # value initialisations
    nummol = max(20, length(mc.indices))
    accepted = false
    translation_dmax = 1.3u"Å"

    # main loop
    for idx_cycle in 1:simu.ncycles
        temperature = simu.temperatures[idx_cycle]
        accepted_translations = 0
        attempted_translations = 0

        println("<<< MAIN TASK STARTED")
        flush(stdout)

        for idnummol in 1:nummol
            # println("<<< MAIN TASK A ", idnummol, '/', nummol)
            # flush(stdout)
            # choose the species on which to attempt a move
            idx = choose_random_species(mc)

            # println("<<< MAIN TASK B ", idnummol)
            # flush(stdout)
            # currentposition is the position of that species
            currentposition = @view mc.step.positions[mc.step.posidx[idx[1]][idx[2]]]

            attempted_translations += 1
            # newpos is the position after the trial move
            newpos = random_translation(currentposition, translation_dmax)
            # if the previous move was accepted, wait for mc to be completely up to date
            # before computing

            # println("<<< MAIN TASK C ", idnummol)
            # flush(stdout)

            i, j = idx
            ij = mc.offsets[i]+j

            molpos = mc.step.posidx[i][j]
            positions = @view mc.step.positions[molpos]
            # println("<<< MAIN TASK D ", idnummol)
            # flush(stdout)
            singlevdw_before = single_contribution_vdw(mc.step, (i,j), positions)
            # println("<<< MAIN TASK E ", idnummol)
            # flush(stdout)
            singlereciprocal_after = single_contribution_ewald(mc.ewald, ij, newpos)
            # println("<<< MAIN TASK F ", idnummol)
            # flush(stdout)
            singlereciprocal_before = single_contribution_ewald(mc.ewald, ij, nothing)
            # println("<<< MAIN TASK G ", idnummol)
            # flush(stdout)
            fer_before = framework_interactions(mc, i, positions)
            # println("<<< MAIN TASK H ", idnummol)
            # flush(stdout)
            fer_after = framework_interactions(mc, i, newpos)
            # println("<<< MAIN TASK I ", idnummol)
            # flush(stdout)
            singlevdw_after = single_contribution_vdw(mc.step, (i,j), newpos)
            # println("<<< MAIN TASK J ", idnummol)
            # flush(stdout)
            before = MCEnergyReport(fer_before, singlevdw_before, singlereciprocal_before)
            after = MCEnergyReport(fer_after, singlevdw_after, singlereciprocal_after)

            accepted_translations += compute_accept_move(before, after, temperature)
        end

        # end of cycle

        @show accepted_translations
        report_now = simu.printevery > 0 && idx_cycle%simu.printevery == 0
        if !(simu.record isa Returns) || report_now
            println("<<< MAIN TASK K")
            flush(stdout)
            o = SimulationStep(mc.step, :output)
            if report_now
                push!(reports, energy)
            end
            println("<<< MAIN TASK L")
            flush(stdout)
            if !(simu.record isa Returns)
                ocomplete = needcomplete(simu.record) ? SimulationStep(o, :complete_output) : o
                t1 = time()
                println("<<< MAIN TASK ABOUT TO RECORD")
                flush(stdout)
                simu.record(deepcopy(ocomplete), energy, idx_cycle, deepcopy(mc), simu)
                println("<<< MAIN TASK RECORD SENT")
                flush(stdout)
                t2 = time()
                @show t2 - t1
                println("<<< MAIN TASK sleeping now for ", waittime, "s...")
                flush(stdout)
                sleep(waittime)
                println("<<< MAIN TASK resuming!")
                flush(stdout)
            end
        end

        translation_ratio = accepted_translations / attempted_translations
        if !isnan(translation_ratio)
            translation_dmax *= 1 + (translation_ratio - 0.5)*(1.0 + 3*translation_ratio)*100/(99+idx_cycle)
        end
    end

    println("MAIN TASK FINISHED")
    flush(stdout)
    fetch(simu.record)
    push!(reports, energy)
    reports
end
















function run_montecarlo_sub!(mc::MonteCarloSetup, simu::SimulationSetup)
    # report
    energy = baseline_energy(mc)
    reports = [energy]

    mkpath(simu.outdir)

    startstep = SimulationStep(mc.step, :output)
    simu.record(startstep, energy, 0, mc, simu)

    # value initialisations
    nummol = max(20, length(mc.indices))
    accepted = false
    translation_dmax = 1.3u"Å"

    # main loop
    for idx_cycle in 1:simu.ncycles
        temperature = simu.temperatures[idx_cycle]
        accepted_translations = 0
        attempted_translations = 0

        for idnummol in 1:nummol
            # choose the species on which to attempt a move
            idx = choose_random_species(mc)

            # currentposition is the position of that species
            currentposition = @view mc.step.positions[mc.step.posidx[idx[1]][idx[2]]]

            attempted_translations += 1
            # newpos is the position after the trial move
            newpos = random_translation(currentposition, translation_dmax)
            # if the previous move was accepted, wait for mc to be completely up to date
            # before computing

            i, j = idx
            ij = mc.offsets[i]+j

            molpos = mc.step.posidx[i][j]
            positions = @view mc.step.positions[molpos]
            singlevdw_before = single_contribution_vdw(mc.step, (i,j), positions)
            singlereciprocal_after = single_contribution_ewald(mc.ewald, ij, newpos)
            singlereciprocal_before = single_contribution_ewald(mc.ewald, ij, nothing)
            fer_before = framework_interactions(mc, i, positions)
            fer_after = framework_interactions(mc, i, newpos)
            singlevdw_after = single_contribution_vdw(mc.step, (i,j), newpos)
            before = MCEnergyReport(fer_before, singlevdw_before, singlereciprocal_before)
            after = MCEnergyReport(fer_after, singlevdw_after, singlereciprocal_after)

            accepted = compute_accept_move(before, after, temperature)
            if accepted
                # update_mc!(mc, idx, newpos)
                diff = after - before
                accepted_translations += 1
                energy += diff
            end
        end

        # end of cycle

        report_now = simu.printevery > 0 && idx_cycle%simu.printevery == 0
        if !(simu.record isa Returns) || report_now
            o = SimulationStep(mc.step, :output)
            if report_now
                push!(reports, energy)
            end
            if !(simu.record isa Returns)
                ocomplete = needcomplete(simu.record) ? SimulationStep(o, :complete_output) : o
                simu.record(ocomplete, energy, idx_cycle, mc, simu)
                idx_cycle%100 == 0 && println("sub-task reached checkpoint ", idx_cycle÷100, '/', simu.ncycles÷100)
                flush(stdout)
                yield()
            end
        end

        translation_ratio = accepted_translations / attempted_translations
        if !isnan(translation_ratio)
            translation_dmax *= 1 + (translation_ratio - 0.5)*(1.0 + 3*translation_ratio)*100/(99+idx_cycle)
        end
    end

    println("sub-task FINISHED")
    flush(stdout)
    # fetch(simu.record)
    push!(reports, energy)
    reports
end
