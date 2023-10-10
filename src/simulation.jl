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
- `o` contains the information on the positions of the species at this point. It must not
  be modified by the `record` function.
- `e` is the energy at this point.
- `k` is the cycle number (between 0 and `cycles`). `k == 0` corresponds to the
  initialization of the simulation.
- `mc` is the input `MonteCarloSetup`. Note that the `positions` field may be incorrect,
  use that of `o`. It must not be modified by the `record` function.
- `simu` is this `SimulationSetup`. It must not be modified by the `record` function.

Additionally, if a call to `record` can spawn asynchronous calls, a method for `Base.fetch`
should be implemented such that `Base.fetch(record)` blocks until all asynchronous calls
are done.

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
```

To perform a simulation that runs for 1000 cycles with logarithmically increasing
temperatures, no trajectory stored but a record the average weighted distance to point
`(0.5, 0.9, 1.3)` (in Å) and energies stored every 2 cycles, use
`SimulationSetup((k,n)->log1p(k/n)*300u"K"/log(2), 1000, "", 2, AverageDistanceTemperatureRecord([0.5, 0.9, 1.3]u"Å"))`
"""
@kwdef struct SimulationSetup{Trecord}
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
        new{typeof(record)}(temperatures, ncycles, outdir, printevery, record)
    end
end
function SimulationSetup(T, ncycles; kwargs...)
    SimulationSetup(; temperatures=T, ncycles, kwargs...)
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
    if mc.step.parallel
        record_task = @spawn $simu.record(SimulationStep($mc.step, :output), $energy, 0, $mc, $simu)
    else
        simu.record(SimulationStep(mc.step, :output), energy, 0, mc, simu) === :stop && return reports
        record_task = @spawn nothing # for type stability
    end

    if simu.ncycles == 0 # single-point computation
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
    translation_dmax = 1.3u"Å"
    rotation_θmax = 30.0u"°"

    # main loop
    for idx_cycle in 1:simu.ncycles
        temperature = simu.temperatures[idx_cycle]
        accepted_translations = 0
        accepted_rotations = 0
        attempted_translations = 0
        attempted_rotations = 0

        for idnummol in 1:nummol
            # choose the species on which to attempt a move
            idx = choose_random_species(mc)
            ffidxi = mc.step.ffidx[idx[1]]

            # currentposition is the position of that species
            # currentposition is either a Vector or a @view, that's OK
            currentposition = (accepted&(old_idx==idx)) ? oldpos : @view mc.step.positions[mc.step.posidx[idx[1]][idx[2]]]

            istranslation = false
            isrotation = false
            # newpos is the position after the trial move
            newpos = if length(ffidxi) == 1
                istranslation = true
                attempted_translations += 1
                random_translation(currentposition, translation_dmax)
            else
                if rand() < 0.5
                    istranslation = true
                    attempted_translations += 1
                    random_translation(currentposition, translation_dmax)
                else
                    isrotation = true
                    attempted_rotations += 1
                    random_rotation(currentposition, rotation_θmax, mc.bead[idx[1]])
                end
            end

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
            k = mc.offsets[i]+j

            if old_idx == idx
                if accepted
                    before = after
                # else
                #     before = fetch(before_task) # simply keep its previous value
                end
                @spawnif mc.step.parallel begin
                    singlereciprocal = @spawn single_contribution_ewald($mc.ewald, $k, $newpos)
                    fer = @spawn framework_interactions($mc, $i, $newpos)
                    singlevdw = single_contribution_vdw(mc.step, idx, newpos)
                    after = MCEnergyReport(fetch(fer), singlevdw, fetch(singlereciprocal))
            end
            else
                molpos = mc.step.posidx[i][j]
                positions = @view mc.step.positions[molpos]
                @spawnif mc.step.parallel begin
                    singlevdw_before_task = @spawn single_contribution_vdw($mc.step, $(i,j), $positions)
                    singlereciprocal_after = @spawn single_contribution_ewald($mc.ewald, $k, $newpos)
                    singlereciprocal_before = @spawn single_contribution_ewald($mc.ewald, $k, nothing)
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
                if istranslation
                    accepted_translations += 1
                elseif isrotation
                    accepted_rotations += 1
                end
                if abs(Float64(diff.framework)) > 1e50 # an atom left a blocked pocket
                    wait(running_update)
                    energy = baseline_energy(mc) # to avoid underflows
                else
                    energy += diff
                end
            end
        end
        # end of cycle
        translation_ratio = accepted_translations / attempted_translations
        report_now = simu.printevery > 0 && idx_cycle%simu.printevery == 0
        if !(simu.record isa Returns) || report_now
            accepted && wait(running_update)
            o = SimulationStep(mc.step, :output)
            if report_now
                put!(output, o)
                push!(reports, energy)
            end
            if !(simu.record isa Returns)
                if mc.step.parallel
                    # ret_record cannot be inferred, but that's fine
                    ret_record = fetch(record_task)
                    ret_record === :stop && break
                    record_task = @spawn $simu.record($o, $energy, $idx_cycle, $mc, $simu)
                else
                    simu.record(o, energy, idx_cycle, mc, simu) === :stop && break
                end
            end
        end

        rotation_ratio = accepted_rotations / attempted_rotations
        if !isnan(translation_ratio)
            translation_dmax *= 1 + (translation_ratio - 0.5)*(1.0 + 3*translation_ratio)*100/(99+idx_cycle)
        end
        if !isnan(rotation_ratio)
            rotation_θmax = (idx_cycle-1)/idx_cycle*rotation_θmax + rotation_ratio/idx_cycle*120u"°"
        end
    end
    wait(record_task)
    fetch(simu.record)
    accepted && wait(running_update)
    push!(reports, energy)
    put!(output, SimulationStep(mc.step, :zero))
    wait(output_task[])
    reports
end
