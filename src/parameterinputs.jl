using Base.Threads


## Record inputs
abstract type RecordFunction <: Function end

mutable struct RMinimumEnergy{N,T} <: RecordFunction
    mine::Float64
    minpos::SimulationStep{N,T}

    RMinimumEnergy{N,T}() where {N,T} = new{N,T}(Inf)
    function RMinimumEnergy(mine::Float64, minpos::SimulationStep{N,T}) where {N,T}
        new{N,T}(mine, minpos)
    end
end
function (record::RMinimumEnergy)(o::SimulationStep, e::Float64, k::Int,
                                  mc::MonteCarloSetup, simu::SimulationSetup)
    if k == simu.ncycles
        output_restart(record.minpos)
    end
    nothing
end


struct ShootingStarMinimizer{N,T} <: RecordFunction
    every::Int
    length::Int
    positions::Vector{SimulationStep{N,T}}
    energies::Vector{Float64}
    lb::LoadBalancer{Tuple{Int,MonteCarloSetup{N,T},SimulationSetup{RMinimumEnergy{N,T}}}}
end
function ShootingStarMinimizer{N}(; length::Int=100, every::Int=1) where {N}
    T = typeof_psystem(Val(N))
    positions = Vector{SimulationStep{N,T}}(undef, 0)
    energies = Vector{Float64}(undef, 0)
    # lb = LoadBalancer{Tuple{Int,MonteCarloSetup{N,T},SimulationSetup{RMinimumEnergy{N,T}}}}(nthreads()) do (ik, newmc, newsimu)
    lb = LoadBalancer{Tuple{Int,MonteCarloSetup{N,T},SimulationSetup{RMinimumEnergy{N,T}}}}(7) do (ik, newmc, newsimu)
        let ik=ik, newmc=newmc, newsimu=newsimu, positions=positions, energies=energies
            current_task().storage = ik
            run_montecarlo_sub!(newmc, newsimu)
        end
    end
    ShootingStarMinimizer(every, length, positions, energies, lb)
end
function initialize_record!(star::T, simu::SimulationSetup{T}) where {T <: ShootingStarMinimizer}
    n = simu.ncycles ÷ star.every
    resize!(star.positions, n)
    resize!(star.energies, n)
end

function (star::ShootingStarMinimizer)(o::SimulationStep, e::Float64, k::Int, mc::MonteCarloSetup, _)
    k ≤ 0 && return
    ik, r = divrem(k, star.every)
    r == 0 || return
    newmc = MonteCarloSetup(mc, o)
    recordminimum = RMinimumEnergy(e, o)
    newsimu = SimulationSetup(; T=300, ncycles=star.length, printevery=0, record=recordminimum)
    put!(star.lb, (ik, newmc, newsimu))
    nothing
end

Base.fetch(x::ShootingStarMinimizer) = wait(x.lb)


struct RainfallMinimizer{N,T} <: RecordFunction
    every::Int
    length::Int
    positions::Vector{SimulationStep{N,T}}
    energies::Vector{Float64}
    tasks::Vector{Task}
end

function RainfallMinimizer{N}(; length::Int=100, every::Int=1) where {N}
    T = typeof_psystem(Val(N))
    positions = Vector{SimulationStep{N,T}}(undef, 0)
    energies = Vector{Float64}(undef, 0)
    tasks = Vector{Task}(undef, 0)
    RainfallMinimizer{N,T}(every, length, positions, energies, tasks)
end

function initialize_record!(rain::T, simu::SimulationSetup{T}) where {T<:RainfallMinimizer}
    n = simu.ncycles ÷ rain.every
    resize!(rain.positions, n)
    resize!(rain.energies, n)
    resize!(rain.tasks, n)
end

function (rain::RainfallMinimizer)(o::SimulationStep, e::Float64, k::Int, mc::MonteCarloSetup, _)
    k ≤ 0 && return
    ik, r = divrem(k, rain.every)
    r == 0 || return
    newmc = MonteCarloSetup(mc, o)
    recordminimum = RMinimumEnergy(e, o)
    newsimu = SimulationSetup(; T=300, ncycles=rain.length, printevery=0, record=recordminimum)
    task = let newmc=newmc, newsimu=newsimu, rain=rain, ik=ik
        Task(() -> begin
            run_montecarlo!(newmc, newsimu)
            rain.positions[ik] = newsimu.record.minpos
            rain.energies[ik] = newsimu.record.mine
            nothing
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


function reconstitute_trace(path::AbstractString, skip, keep)
    numdirs = length(readdir(path)) - 1
    energies = [deserialize(joinpath(path, string(i), "energies.serial")) for i in 1:numdirs if begin
        x = ispath(joinpath(path, string(i), "energies.serial"))
        x || @warn "Subfolder $i does not contain energies.serial"
        x
    end]
    numdirs = length(energies)
    trace = Vector{Vector{Float64}}(undef, numdirs)
    for (i, energy) in enumerate(energies)
        n = length(energy)
        start = skip isa Integer ? skip : floor(Int, skip*n)
        start += (start == 0)
        stop = keep isa Integer ? start + keep - 1 : start + round(Int, keep*(n-start+1))
        stop -= (stop == length(energy) + 1)
        trace[i] = Float64.(@view energy[start:stop])
    end
    trace
end
