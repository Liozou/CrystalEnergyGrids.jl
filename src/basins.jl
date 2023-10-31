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
    local_minima(grid::Array{Float64,3}, tolerance=1e-2; lt=(<))

Find the positions of the local minima of an energy grid.

If `tolerance ≥ 0`, only return the smallest minima within the relative tolerance compared
to the energy minimum. `tolerance == 0` corresponds to only returning the minimum while
`tolerance == 1` corresponds to returning all the local minima with the same sign as the
minimum.

Use a negative `tolerance` to return all local minima.

`lt` is the "lower than" operator: use `lt=(>)` to compute local maxima.
"""
function local_minima(grid::Array{Float64,3}, tolerance=1e-2; lt=(<))
    a1, a2, a3 = size(grid)
    localmins_t = [CartesianIndex{3}[] for _ in 1:nthreads()]
    @threads :static for i3 in 1:a3
        for i2 in 1:a2, i1 in 1:a1
            val = grid[i1,i2,i3]
            if lt(val, grid[mod1(i1-1,a1),i2,i3]) &&
               lt(val, grid[mod1(i1+1,a1),i2,i3]) &&
               lt(val, grid[i1,mod1(i2-1,a2),i3]) &&
               lt(val, grid[i1,mod1(i2+1,a2),i3]) &&
               lt(val, grid[i1,i2,mod1(i3-1,a3)]) &&
               lt(val, grid[i1,i2,mod1(i3+1,a3)])
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
    circular_distance(a::Int, (m, n))

Smallest distance between the images of two integers `n` and `m` placed on a periodic axe
of length `a`.
"""
function circular_distance(a::Int, (m, n))
    b = a ÷ 2
    i, j = minmax(m, n)
    i < b < j && return min(j - i, i - j + a)
    return j - i
end


"""
    local_basins(grid::Array{Float64,3}, minima::Vector{CartesianIndex{3}}; lt=(<), merge_manhattan=0, skip=Returns(false))
    local_basins(grid::Array{Float64,3}, tolerance::Float64=-Inf; lt=(<), merge_manhattan=0, skip=Returns(false))

Return a list `basinsets` whose elements are the local attraction basins around each energy
minimum given in `minima`. Each element is a `Set` of grid point `(i,j,k)` such that there
is a path that only decreases in energy when going from the corresponding energy minimum to
`(i,j,k)`.

For each pair of minima whose Manhattan distance on the grid is lower or equal to
`merge_manhattan`, their corresponding basins are merged into a single one.

Grid points whose energy `e` is such that `skip(e)` are not included in the basins. If such
a point is an element of `minima`, the corresponding basin will be skipped as well.

If `minima` is not provided, it will be computed from [`local_minima`](@ref) with the given
`tolerance` and keyword arguments.
"""
function local_basins(grid::Array{Float64,3}, minima::Vector{CartesianIndex{3}};
                      lt=(<), merge_manhattan=0, skip=Returns(false))
    a, b, c, = dims = size(grid)
    nminima = length(minima)
    basinsets = Vector{Set{NTuple{3,Int}}}(undef, nminima)
    @threads for itask in 1:nminima
        start = minima[itask]
        visited = falses(dims)
        istart, jstart, kstart = Tuple(start)
        if skip(grid[istart, jstart, kstart])
            basinsets[itask] = Set{NTuple{3,Int}}() # filtered below
            continue
        end
        visited[istart, jstart, kstart] = true
        Q = [Tuple(start)]
        for u in Q
            iu, ju, ku = u
            gu = grid[iu,ju,ku]
            for v in PeriodicNeighbors(dims, u)
                iv, jv, kv = v
                visited[iv, jv, kv] && continue
                gv = grid[iv, jv, kv]
                skip(gv) && continue
                if lt(gu, grid[iv, jv, kv])
                    push!(Q, v)
                    visited[iv, jv, kv] = true
                end
            end
        end
        basinsets[itask] = Set(Q)
    end
    filter!(!isempty, basinsets)
    if merge_manhattan > 0
        unions = collect(1:nminima)
        for i1 in 1:nminima
            a1, b1, c1 = Tuple(minima[i1])
            for i2 in (i1+1):nminima
                a2, b2, c2 = Tuple(minima[i2])
                if circular_distance(a, (a1,a2)) + circular_distance(b, (b1,b2)) + circular_distance(c, (c1,c2)) ≤ merge_manhattan
                    unions[i2] = unions[i1]
                end
            end
        end
        todelete = Int[]
        for i in 1:nminima
            j = unions[i]
            i == j && continue
            union!(basinsets[j], basinsets[i])
            push!(todelete, i)
        end
        deleteat!(basinsets, todelete)
    end
    basinsets
end

function local_basins(grid::Array{Float64,3}, tolerance::Float64=-Inf;
                      lt=(<), merge_manhattan=0, skip=Returns(false))
    local_basins(grid, local_minima(grid, tolerance; lt); lt, merge_manhattan, skip)
end


# struct FrontierTree{T}
#     children::Union{Vector{T}, Vector{Tuple{Int,FrontierTree}}}
# end
# Base.getindex(x::FrontierTree) = x.children
# function Base.getindex(x::FrontierTree{T}, i, I...) where T
#     children = x.children
#     if children isa Vector{Tuple{Int,FrontierTree}}
#         for (k, f) in children
#             k == i && return f[I...]
#             k > i && break
#         end
#     end
#     return T[]
# end

# struct BasinDecompositon{U}
#     tree::FrontierTree{NTuple{3,Int}}
#     nodes::Array{Vector{Int},3}
# end
# Base.getindex(x::BasinDecompositon, I...) = x.tree[I...]

"""
    decompose_basins(grid::Array{Float64,3}, basinsets::Vector{Set{NTuple{3,Int}}})
    decompose_basins(grid::Array{Float64,3}, tolerance::Float64=-Inf; lt=(<), merge_manhattan=0, skip=Returns(false))

Decompose the provided `basinsets` into a grid `nodes` where `nodes[i,j,k]` is the list of
basins to which the grid point `(i,j,k)` belongs to. Return `(points, nodes)` where
`points` is the list of grid point for which `nodes` is not empty.

If `basinsets` is not provided, it will be computed from [`local_basins`](@ref) with the
given `tolerance` and keyword arguments.
"""
function decompose_basins(grid::Array{Float64,3}, basinsets::Vector{Set{NTuple{3,Int}}})
    nodes = fill(Int[], size(grid))
    n = length(basinsets)
    refvectors = [[i] for i in 1:n]
    refidx = Dict{Vector{Int},Int}([[i] => i for i in 1:n])
    pointsset = Set{NTuple{3,Int}}()
    for (iset, set) in enumerate(basinsets)
        union!(pointsset, set)
        for (i, j, k) in set
            node = nodes[i,j,k]
            push!(node, iset)
            ref = get(refidx, node, nothing)
            if ref isa Nothing
                newref = copy(node)
                push!(refvectors, newref)
                refidx[newref] = length(refvectors)
                ref = length(refvectors)
            end
            pop!(node)
            nodes[i,j,k] = refvectors[ref]
        end
    end
    points = collect(pointsset)
    sort!(points; lt=((i1,j1,k1),(i2,j2,k2)) -> (k1<k2 || (k1==k2&&(j1<j2 || (j1==j2 && i1<i2)))))
    points, nodes
end

function decompose_basins(grid::Array{Float64,3}, tolerance::Float64=-Inf;
                          lt=(<), merge_manhattan=0, skip=Returns(false))
    decompose_basins(grid, local_basins(grid, tolerance; lt, merge_manhattan, skip))
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

function compute_levels(grid::Array{Float64,4}, mine, maxe, T=300)
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

# function export_levels(path, grid::Array{Float64,4}, mine, maxe, framework, T=300)
#     basins = compute_levels(grid, mine, maxe, T)
#     newgrid = zeros(size(grid)[2:end]...)
#     newgrid .= 0.0
#     for l in basins, x in l
#         newgrid[x...] = 1.0
#     end
#     output_cube(path, newgrid, framework)
# end
