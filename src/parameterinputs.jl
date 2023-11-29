using Base.Threads


struct ShootingStarMinimizer{N,T}
    every::Int
    length::Int
    lb::LoadBalancer{MonteCarloSetup{N,T}}
end
function ShootingStarMinimizer{N}(; length::Int=100, every::Int=1) where {N}
    T = typeof_psystem(Val(N))
    # lb = LoadBalancer{Tuple{Int,MonteCarloSetup{N,T},SimulationSetup{RMinimumEnergy{N,T}}}}(nthreads()) do (ik, newmc, newsimu)
    lb = LoadBalancer{MonteCarloSetup{N,T}}(output_restart, 7)
    ShootingStarMinimizer(every, length, lb)
end

function (star::ShootingStarMinimizer)(o::SimulationStep, e::Float64, k::Int, mc::MonteCarloSetup, _)
    k ≤ 0 && return
    _, r = divrem(k, star.every)
    r == 0 || return
    newmc = MonteCarloSetup(mc, o)
    put!(star.lb, newmc)
    nothing
end

Base.fetch(x::ShootingStarMinimizer) = wait(x.lb)
