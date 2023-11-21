module CrystalEnergyGrids

using Base.Threads

struct LoadBalancer{T}
    channel::Channel{T}
    tasks::Vector{Task}
    event::Event
end

function LoadBalancer{T}(f, n::Integer=nthreads()-1) where T
    event::Event = Event(true)
    channel::Channel{T} = Channel{T}(Inf)
    tasks = [errormonitor(@spawn while true
        x = take!($channel)
        $f(x)
        notify($event)
    end) for _ in 1:n]
    LoadBalancer{T}(channel, tasks, event)
end

Base.put!(lb::LoadBalancer, x) = put!(lb.channel, x)

function Base.wait(lb::LoadBalancer)
    while !isempty(lb.channel)
        wait(lb.event)
    end
end

struct Setup
    positions::Vector{Float64}
end

function run(setup::Setup, lb::LoadBalancer)
    for _ in 1:10
        pos = copy(setup.positions)
        put!(lb, pos)
        yield()
    end
    wait(lb)
end


function ignore_args(_)
    Float64[rand()]
end

function main()
    for _ in 1:100000
        setup = Setup(rand(300))
        lb = LoadBalancer{Vector{Float64}}(ignore_args, 9)
        run(setup, lb)
    end
end

end
