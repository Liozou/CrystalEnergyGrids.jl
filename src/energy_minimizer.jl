using Base.Threads

struct PeriodicNeighbors{N}
    a::NTuple{N,Int}
    i::NTuple{N,Int}

    function PeriodicNeighbors{N}(a::NTuple{N,Int}, i::NTuple{N,Int}) where {N}
        if any(≤(2), a)
            error(lazy"Cannot construct a PeriodicNeighbors with dims $a")
        end
        new{N}(a,i)
    end
end
PeriodicNeighbors(a::NTuple{N,Int}, i::NTuple{N,Int}) where {N} = PeriodicNeighbors{N}(a, i)
Base.length(::PeriodicNeighbors{N}) where {N} = 2*N
Base.eltype(::PeriodicNeighbors{N}) where {N} = NTuple{N,Int}
function Base.iterate(x::PeriodicNeighbors{N}, state=(0,true)) where N
    i, ispos = state
    if ispos
        i += 1
        i == N+1 && return nothing
    end
    return (ntuple(N) do j
        j == i ? mod1(x.i[j] + ifelse(ispos, 1, -1), x.a[j]) : x.i[j]
    end, (i, !ispos))
end

function local_minima(grid::Array{Float64,3})
    a1, a2, a3 = size(grid)
    ret = CartesianIndex{3}[]
    @threads for i3 in 1:a3
        for i2 in 1:a2, i1 in 1:a1
            val = grid[a1,a2,a3]
            if val < grid[mod1(i1-1,a1),a2,a3] &&
               val < grid[mod1(i1+1,a1),a2,a3] &&
               val < grid[a1,mod1(i2-1,a2),a3] &&
               val < grid[a1,mod1(i2+1,a2),a3] &&
               val < grid[a1,a2,mod1(i3-1,a3)] &&
               val < grid[a1,a2,mod1(i3+1,a3)]
                push!(ret, CartesianIndex(i1, i2, i3))
            end
        end
    end
    ret
end

function attraction_basins(grid::Array{Float64,3}, locals=local_minima(grid))
    basins = Vector{CartesianIndex{3}}[]
    a1, a2, a3 = size(grid)
    tovisit = falses(a1, a2, a3)
end