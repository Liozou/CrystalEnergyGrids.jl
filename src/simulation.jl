export SimulationSetup

"""
    SimulationSetup{T}

Definition of the simulation setup, which must provide the following information:
- `T` the temperature. Either a number (in K), a list of temperatures (in K, one per cycle)
  or a function. If a function, it must take two parameters: `k` the number of the current
  cycle and `n` the total number of cycles, and return a temperature (in K).
- `ncycles` the number of cycles of the simulation. Each cycle consists in `N` steps, where
  a step is a Monte-Carlo movement attempt on a random species, and `N` is the number of
  species.
- `outdir` the path to a directory in which the result is stored. If not provided (or set
  to an empty string), no record will be stored on file.
- `printevery`: one record of the positions is stored every `printevery` cycles. Default to
  1000.
"""
struct SimulationSetup
    temperatures::Vector{typeof(1.0u"K")}
    ncycles::Int
    outdir::String
    printevery::Int
end
function SimulationSetup(T, ncycles::Int, outdir::String, printevery::Int)
    temperatures = if T isa Number
        fill(T, ncycles)
    elseif T isa Vector
        convert(Vector{typeof(1.0u"K")}, T)
    else
        T.(1:ncycles, ncycles)
    end
    SimulationSetup(temperatures, ncycles, outdir, printevery)
end
SimulationSetup(T, ncycles, outdir="") = SimulationSetup(T, ncycles, outdir, 1000)

abstract type TemperatureFunction <: Function end

struct TRamp <: TemperatureFunction
    start::typeof(1.0u"K")
    Δ::typeof(1.0u"K")

    function TRamp(start::typeof(1.0u"K"), Δ::typeof(1.0u"K"), ::Nothing)
        new(start, Δ)
    end
end
function TRamp(start, stop)
    TRamp(convert(typeof(1.0u"K"), start), convert(typeof(1.0u"K"), stop-start), nothing)
end
function (tf::TRamp)(k, n)
    tf.start + tf.Δ*(k-1)/(n-1)
end

struct TPoly <: TemperatureFunction
    ramp::TRamp
    γ::Float64
end
function (tf::TPoly)(k, n)
    tf.ramp.start + tf.ramp.Δ*((k-1)/(n-1))^tf.γ
end

struct TSteps <: TemperatureFunction
    ramp::TRamp
    nstepsp1::Int
end
function (tf::TSteps)(k, n)
    tf.ramp.start + tf.Δ*floor(tf.nstepsp1*(k-1)/(n-1))/tf.nstepsp1
end

struct TParts{N,T} <: TemperatureFunction
    fractions::NTuple{N,Float64} # sorted between 0 and 1, the last element is always 1
    parts::T
end
function (tf::TParts{N})(k, n) where N
    ρ = (k-1)/(n-1)
    for i in 1:N
        frac = tf.fractions[i]
        if frac ≥ ρ
            num_start = i==1 ? -1 : floor(Int, n*tf.fractions[i-1])-1
            num_stop = floor(Int, n*frac)
            num_current = floor(Int, n*ρ)
            return tf.parts[i](num_current - num_start, num_stop - num_start)
        end
    end
    @assert false
    0.0u"K"
end

function TAnneal(bottom::typeof(1.0u"K"), top::typeof(1.0u"K"), wait::Float64)
    TParts(((1.0 - wait)/2, (1.0 + wait)/2, 1.0), (
        TRamp(bottom, top),
        Returns(top),
        TRamp(top, bottom)
    ))
end
function TAnneal(sequence, wait::Number)
    N = length(sequence)
    TParts(ntuple(i -> begin
                        j, o = fldmod(i, 2)
                        (j + o*(1-wait))/(N-1)
                       end, 2*(N-1)),
           ntuple(i -> begin
                        j, o = fldmod1(i, 2)
                        seqjp1 = sequence[j+1]
                        o == 1 ? TRamp(sequence[j], seqjp1) : Returns(seqjp1)
                       end, 2*(N-1))
    )
end



"""
    run_montecarlo!(mc::MonteCarloSetup, T, nsteps::Int)

Run a Monte-Carlo simulation at temperature `T` (given in K) during `nsteps`.
"""
function run_montecarlo!(mc::MonteCarloSetup, simu::SimulationSetup)
    nummol = max(20, length(mc.indexof))
    energy = baseline_energy(mc)
    reports = [energy]
    running_update = @spawn nothing
    local oldpos::Vector{SVector{3,typeof(1.0u"Å")}}
    old_idx = (0,0)
    local before::MCEnergyReport
    local after::MCEnergyReport
    accepted = false
    mkpath(simu.outdir)
    output = pdb_output_handler(isempty(simu.outdir) ? "" : joinpath(simu.outdir, "trajectory.pdb"), mc.cell)
    launch_output_task = @spawn put!(output, OutputSimulationStep(mc))
    for k in 1:simu.ncycles
        temperature = simu.temperatures[k]
        for idnummol in 1:nummol
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
                before = fetch(before_task)::MCEnergyReport
            end
            old_idx = idx
            accepted = compute_accept_move(before, after, temperature)
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
                # shadow = MonteCarloSetup(mc)
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
    end
    accepted && wait(running_update)
    push!(reports, energy)
    wait(launch_output_task)
    put!(output, OutputSimulationStep(mc, nothing))
    reports
end
