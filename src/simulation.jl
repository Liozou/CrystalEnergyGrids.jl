struct SimulationSetup{Trecord}
    outdir::String
    record::Trecord
end


"""
    run_montecarlo!(mc::MonteCarloSetup, simu::SimulationSetup)

Run a Monte-Carlo simulation. Return the list of energies of the system during the
simulation.

See [`MonteCarloSetup`](@ref) for the definition of the system and
[`SimulationSetup`](@ref) for the parameters of the simulation.
"""
function run_montecarlo!(mc::MonteCarloSetup, simu::SimulationSetup)
    mkpath(simu.outdir)

    for idx_cycle in 1:10
        pos = copy(mc.positions)
        simu.record(pos, rand(mc.rng), idx_cycle, mc, simu)
        yield()
    end
    fetch(simu.record)
end


function ignore_args(::Nothing)
    Float64[rand()]
end
