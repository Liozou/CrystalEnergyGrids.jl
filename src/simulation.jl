export SimulationSetup

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





"""
    run_montecarlo!(mc::MonteCarloSetup, simu::SimulationSetup)

Run a Monte-Carlo simulation. Return the list of energies of the system during the
simulation.

See [`MonteCarloSetup`](@ref) for the definition of the system and
[`SimulationSetup`](@ref) for the parameters of the simulation.
"""
function run_montecarlo!(mc::MonteCarloSetup, simu::SimulationSetup)
    # record and outputs
    mkpath(simu.outdir)

    # main loop
    for idx_cycle in 1:10
        # end of cycle
        if idx_cycle%simu.printevery == 0
            pos = copy(mc.positions)
            simu.record(pos, rand(mc.rng), idx_cycle, mc, simu)
            yield()
        end
    end
    fetch(simu.record)
end


function ignore_args(mc::MonteCarloSetup, simu::SimulationSetup)
    # energy initialization
    Float64[rand()]
end
