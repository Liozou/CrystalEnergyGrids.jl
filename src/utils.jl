"""
    LoadBalancer{T}

Channel-like abstraction to balance function calls across a fixed number of tasks.
See [`LoadBalancer{T}(f, n::Integer)`](@ref) to create an instance.
"""
struct LoadBalancer{T}
    channel::Channel{T}
    tasks::Vector{Task}
    busy::Atomic{Int}
    event::Event
end

using Serialization

"""
    LoadBalancer{T}(f, n::Integer=nthreads())

Given a function `f`, create `n` tasks that wait on an input `x::T` to execute `f(x)`.
With `lb = LoadBalancer(f, n)`, use `put!(lb, x)` to send `x` to one of the tasks: this
call is not blocking and will always return immediately, but the `f(x)` call will not occur
until one of the `n` tasks is free. Use `wait(lb)` to wait until all inputs have been
processed.

## Example

```julia
lb = LoadBalancer{Tuple{Int,String}}() do (i,s)
    some_complicated_function(i, s)
end

for i in 1:1000
    put!(lb, (i, "")) # setup works
end

do_something_else()

wait(lb)
```
"""
function LoadBalancer{T}(f, n::Integer=nthreads()-1) where T
    busy::Atomic{Int} = Atomic{Int}(0)
    event::Event = Event(true)
    channel::Channel{T} = Channel{T}(Inf)
    tasks = [(@spawn while true
        x = take!($channel)
        atomic_add!($busy, 1)
        $f(x)
        atomic_sub!($busy, 1)
        notify($event)
    end) for i in 1:n]
    foreach(errormonitor, tasks)
    LoadBalancer{T}(channel, tasks, busy, event)
end
Base.put!(lb::LoadBalancer, x) = put!(lb.channel, x)
function Base.wait(lb::LoadBalancer)
    while !isempty(lb.channel) || lb.busy[] > 0
        wait(lb.event)
    end
end
function Base.close(lb::LoadBalancer)
    if lb.busy[] != 0
        @error "Forcibly closing the LoadBalancer may yield errors from the still-working tasks"
    end
    close(lb.channel)
    lb.busy[] = -max(0, lb.busy[]) # make lb.busy zero or negative to signal the close
    notify(lb.event)
end
