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

"""
    local_minima(grid::Array{Float64,3}, tolerance=1e-2)

Find the positions of the local minima of an energy grid.

If `tolerance ≥ 0`, only return the smallest minima within the relative tolerance compared
to the energy minimum. `tolerance == 0` corresponds to only returning the minimum while
`tolerance == 1` corresponds to returning all the local minima with the same sign as the
minimum.

Use a negative `tolerance` to return all local minima.
"""
function local_minima(grid::Array{Float64,3}, tolerance=1e-2)
    a1, a2, a3 = size(grid)
    localmins_t = [CartesianIndex{3}[] for _ in 1:nthreads()]
    @threads :static for i3 in 1:a3
        for i2 in 1:a2, i1 in 1:a1
            val = grid[i1,i2,i3]
            if val < grid[mod1(i1-1,a1),i2,i3] &&
               val < grid[mod1(i1+1,a1),i2,i3] &&
               val < grid[i1,mod1(i2-1,a2),i3] &&
               val < grid[i1,mod1(i2+1,a2),i3] &&
               val < grid[i1,i2,mod1(i3-1,a3)] &&
               val < grid[i1,i2,mod1(i3+1,a3)]
                push!(localmins_t[threadid()], CartesianIndex(i1, i2, i3))
            end
        end
    end
    localmins = reduce(vcat, localmins_t)
    min_energy = minimum(Base.Fix1(getindex, grid), localmins)
    tolerance < 0 && return localmins
    @assert min_energy < 0
    max_threshold = min_energy + abs(min_energy*tolerance)
    kept = CartesianIndex{3}[]
    for localmin in localmins
        grid[localmin] > max_threshold && continue
        push!(kept, localmin)
    end
    return kept
end

function attraction_basins(grid::Array{Float64,3}, locals=local_minima(grid))
    basins = Vector{CartesianIndex{3}}[]
    a1, a2, a3 = size(grid)
    tovisit = falses(a1, a2, a3)
end