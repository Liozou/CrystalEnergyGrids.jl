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

"""
    energy_levels(egrid::Array{Float64,4}, n)

Given an angular energy grids, return the energies and relative volumes of the `n` lowest
levels (between energy -∞ and 0 K)
"""
function energy_levels(egrid::Array{Float64,4}, n)
    kept = filter(<(0.0), vec(egrid))
    m = minimum(kept)
    hist = zeros(Float64, n)
    γ = n/m
    for x in kept
        idx = floor(Int, x*γ)
        hist[end-idx+(idx==n)] += 1
    end
    λ = -m/(2n)
    hist ./= (length(egrid))
    [m+(2k-1)*λ for k in 1:n], map(x -> iszero(x) ? missing : x, hist)
end

function compute_basins(grid::Array{Float64,4}, mine, maxe, T=300)
    basins = Vector{NTuple{3,Int}}[]
    _, a1, a2, a3 = size(grid)
    visited = falses(a1, a2, a3)
    A = (a1, a2, a3)
    for i3 in 1:a3, i2 in 1:a2, i1 in 1:a1
        visited[i1,i2,i3] && continue
        visited[i1,i2,i3] = true
        ei = meanBoltzmann(view(grid, :, i1, i2, i3), T)
        mine ≤ ei ≤ maxe || continue
        basin = [(i1,i2,i3)]
        Q = [(i1,i2,i3)]
        for I in Q
            for J in PeriodicNeighbors(A, I)
                j1, j2, j3 = J
                visited[j1,j2,j3] && continue
                visited[j1,j2,j3] = true
                ej = meanBoltzmann(view(grid, :, j1, j2, j3), T)
                mine < ej < maxe || continue
                push!(basin, J)
                push!(Q, J)
            end
        end
        push!(basins, basin)
    end
    basins
end

# function export_basins(path, grid::Array{Float64,4}, mine, maxe, framework, T=300)
#     basins = compute_basins(grid, mine, maxe, T)
#     newgrid = zeros(size(grid)[2:end]...)
#     newgrid .= 0.0
#     for l in basins, x in l
#         newgrid[x...] = 1.0
#     end
#     output_cube(path, newgrid, framework)
# end
