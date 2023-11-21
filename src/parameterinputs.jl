## Record inputs
abstract type RecordFunction <: Function end

struct ShootingStarMinimizer <: RecordFunction
    every::Int
    length::Int
    outdir::String
    lb::LoadBalancer{Nothing}
end
function ShootingStarMinimizer(; length::Int=100, every::Int=1, outdir="")
    lb = LoadBalancer{Nothing}(7) do newmc
        let newmc=newmc
            ignore_args(newmc)
        end
    end
    ShootingStarMinimizer(every, length, outdir, lb)
end

function (star::ShootingStarMinimizer)(o::Vector{SVector{3,Float64}}, k::Int, mc::MonteCarloSetup, _)
    k ≤ 0 && return
    k % star.every == 0 || return
    put!(star.lb, nothing)
    nothing
end

Base.fetch(x::ShootingStarMinimizer) = wait(x.lb)







struct RainfallMinimizer <: RecordFunction
    every::Int
    length::Int
    positions::Vector{Vector{SVector{3,Float64}}}
    energies::Vector{Float64}
    tasks::Vector{Task}
end

function RainfallMinimizer(; length::Int=100, every::Int=1)
    positions = Vector{Vector{SVector{3,Float64}}}(undef, 0)
    energies = Vector{Float64}(undef, 0)
    tasks = Vector{Task}(undef, 0)
    RainfallMinimizer(every, length, positions, energies, tasks)
end

function (rain::RainfallMinimizer)(o::Vector{SVector{3,Float64}}, e::Float64, k::Int, mc::MonteCarloSetup, _)
    k ≤ 0 && return
    ik, r = divrem(k, rain.every)
    r == 0 || return
    newmc = MonteCarloSetup(mc, o)
    newsimu = SimulationSetup("", Returns(nothing))
    task = let newmc=newmc, newsimu=newsimu
        Task(() -> begin
            run_montecarlo!(newmc, newsimu)
        end)
    end
    task.sticky = false
    rain.tasks[ik] = task
    errormonitor(task)

    schedule(task)
    nothing
end

function Base.fetch(x::RainfallMinimizer)
    retry = false
    for (i, task) in enumerate(x.tasks)
        if !isassigned(x.tasks, i)
            retry = true
            continue
        end
        wait(task)
    end
    if retry
        foreach(wait, x.tasks)
    end
end
