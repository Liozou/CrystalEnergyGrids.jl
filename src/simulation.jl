struct SimulationSetup
    T::typeof(1.0u"K")
    ncycles::Int
    outdir::String
    printevery::Int
end

"""
    run_montecarlo!(mc::MonteCarloSimulation, T, nsteps::Int)

Run a Monte-Carlo simulation at temperature `T` (given in K) during `nsteps`.
"""
function run_montecarlo!(mc::MonteCarloSimulation, simu::SimulationSetup)
    nummol = max(20, length(mc.indexof))
    energy = baseline_energy(mc)
    reports = [energy]
    running_update = @spawn nothing
    oldpos = SVector{3,typeof(1.0)}[]
    old_idx = (0,0)
    before = 0.0u"K"
    after = 0.0u"K"
    accepted = false
    mkpath(simu.outdir)
    output = pdb_output_handler(joinpath(simu.outdir, "trajectory.pdb"), mc.cell)
    launch_output_task = @spawn put!(output, OutputSimulationStep(mc))
    for k in 1:simu.ncycles, idnummol in 1:nummol
        if idnummol==1 && k%simu.printevery == 0
            launch_output_task = @spawn put!(output, OutputSimulationStep(mc))
            push!(reports, energy)
        end
        idx = choose_random_species(mc)
        currentposition = (accepted&(old_idx==idx)) ? oldpos : mc.positions[idx[1]][idx[2]]
        idxi = mc.idx[idx[1]]
        newpos = if length(idxi) == 1
            random_translation(currentposition, 1.3u"Å")
        else
            if rand() < 0.5
                random_translation(currentposition, 1.3u"Å")
            else
                random_rotation(currentposition, 30.0u"°", mc.bead[idx[1]])
            end
        end
        inblockpocket(mc.speciesblocks[idx[1]], mc.atomblocks, idxi, newpos) && continue
        accepted && wait(running_update)
        if old_idx == idx
            if accepted
                before = after
            # else
            #     before = fetch(before_task) # simply keep its previous value
            end
            after = movement_energy(mc, idx, newpos)
        else
            before_task = @spawn movement_energy(mc, idx)
            after = movement_energy(mc, idx, newpos)
            before = fetch(before_task)
        end
        old_idx = idx
        accepted = compute_accept_move(before, after, simu.T)
        oldpos = newpos
        if accepted
            wait(launch_output_task)
            running_update = @spawn update_mc!(mc, idx, oldpos) # do not use newpos since it can be changed in the next iteration before the Task is run
            diff = after - before
            if abs(Float64(diff.framework)) > 1e50 # an atom left a blocked pocket
                wait(running_update)
                energy = baseline_energy(mc) # to avoid underflows
            else
                energy += diff
            end

            # wait(running_update)
            # shadow = MonteCarloSimulation(mc)
            # if !isapprox(Float64(baseline_energy(shadow)), Float64(energy), rtol=1e-4)
            #     println("discrepancy ! ", Float64(baseline_energy(shadow)), " vs ", Float64(energy))
            #     @show idx, k, idnummol
            #     push!(reports, energy)
            #     return reports
            # end
            # println(Float64(energy), " vs ", Float64(baseline_energy(shadow)))
            # @show shadow.ewald.ctx.Eiks[1][30]
            # display(energy)
            # display(baseline_energy(shadow))
            # println(mc.positions)
        end
    end
    accepted && wait(running_update)
    push!(reports, energy)
    wait(launch_output_task)
    put!(output, OutputSimulationStep(mc, nothing))
    reports
end
