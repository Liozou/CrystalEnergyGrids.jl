using Base.Threads: nthreads, @threads
using ImageFiltering: Kernel, imfilter

struct GridNeighbors{N}
    i::NTuple{N,Int}
end
Base.length(::GridNeighbors{N}) where {N} = 2*N
Base.eltype(::GridNeighbors{N}) where {N} = NTuple{N,Int}

struct _GridNeighborIterator{N}
    is::NTuple{N,Int}
    i::Int
    sign::Int
end
(g::_GridNeighborIterator)(j::Int) = j == g.i ? g.is[j] + g.sign : g.is[j]
function Base.iterate(x::GridNeighbors{N}, state=(0,1)) where N
    i, sign = state
    if sign == 1
        i += 1
        i == N+1 && return nothing
    end
    return (ntuple(_GridNeighborIterator{N}(x.i, i, sign), Val{N}()), (i, -sign))
end

"""
    local_minima(grid::Array{<:Any,3}, tolerance=1e-2; lt=(<))

Find the positions of the local minima of an energy grid, i.e. the indices I such that
`grid[I]` is lower than the 26 grid points forming the smallest cube around it.

If `tolerance ≥ 0`, only return the smallest minima within the relative tolerance compared
to the energy minimum. `tolerance == 0` corresponds to only returning the minimum while
`tolerance == 1` corresponds to returning all the local minima with the same sign as the
minimum.

Use a negative `tolerance` to return all local minima.

`lt` is the "lower than" operator: use `lt=(>)` to compute local maxima.
"""
function local_minima(grid::Array{<:Any,3}, tolerance=1e-2; lt=(<))
    a1, a2, a3 = size(grid)
    localmins = CartesianIndex{3}[]
    # for i3 in 1:a3
    N = nthreads()
    localmins_t = [CartesianIndex{3}[] for _ in 1:N]
    # @loadbalance 0.1 for I in CartesianIndices((a1, a2, a3))
    #     i1, i2, i3 = Tuple(I)
    @loadbalance 0.5 for i3 in 1:a3
        for i2 in 1:a2, i1 in 1:a1
            val = grid[i1,i2,i3]
            m1 = mod1(i1-1,a1)
            lt(val, grid[m1,i2,i3]) || continue
            p1 = mod1(i1+1, a1)
            lt(val, grid[p1,i2,i3]) || continue
            m2 = mod1(i2-1,a2)
            lt(val, grid[i1,m2,i3]) || continue
            lt(val, grid[m1,m2,i3]) || continue
            lt(val, grid[p1,m2,i3]) || continue
            p2 = mod1(i2+1,a2)
            lt(val, grid[i1,p2,i3]) || continue
            lt(val, grid[m1,p2,i3]) || continue
            lt(val, grid[p1,p2,i3]) || continue
            m3 = mod1(i3-1,a3)
            lt(val, grid[i1,i2,m3]) || continue
            lt(val, grid[m1,i2,m3]) || continue
            lt(val, grid[p1,i2,m3]) || continue
            lt(val, grid[i1,m2,m3]) || continue
            lt(val, grid[m1,m2,m3]) || continue
            lt(val, grid[p1,m2,m3]) || continue
            lt(val, grid[i1,p2,m3]) || continue
            lt(val, grid[m1,p2,m3]) || continue
            lt(val, grid[p1,p2,m3]) || continue
            p3 = mod1(i3+1,a3)
            lt(val, grid[i1,i2,p3]) || continue
            lt(val, grid[m1,i2,p3]) || continue
            lt(val, grid[p1,i2,p3]) || continue
            lt(val, grid[i1,m2,p3]) || continue
            lt(val, grid[m1,m2,p3]) || continue
            lt(val, grid[p1,m2,p3]) || continue
            lt(val, grid[i1,p2,p3]) || continue
            lt(val, grid[m1,p2,p3]) || continue
            lt(val, grid[p1,p2,p3]) || continue
            push!(localmins_t[taskid], CartesianIndex(i1, i2, i3))
            # push!(localmins, CartesianIndex(i1, i2, i3))
        end
    end
    # close(lb) # this causes spurious errors to be printed from the `take!`ing tasks
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

# Perform union-find to group together basins whose origins are closer than `cluster`
function clean_cluster!(basinsets::Vector{Set{NTuple{3,Int}}}, minima::Vector{CartesianIndex{3}}, cluster::Real, (a,b,c)::NTuple{3,Int})
    cluster > 0 || return
    nminima = length(minima)
    unions = collect(1:nminima)
    for i1 in 1:nminima
        a1, b1, c1 = Tuple(minima[i1])
        for i2 in (i1+1):nminima
            a2, b2, c2 = Tuple(minima[i2])
            if circular_distance(a, (a1,a2)) + circular_distance(b, (b1,b2)) + circular_distance(c, (c1,c2)) ≤ cluster
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
    nothing
end

function keep_shortest_distance!(protobasinsets::Vector{Tuple{Vector{NTuple{3,Int}},Vector{Int}}}, dims::NTuple{3,Int})
    n = length(protobasinsets)
    basinsets = [Set(Q) for (Q, _) in protobasinsets]
    old = Set{CartesianIndex{3}}()
    new = Set{CartesianIndex{3}}()
    indices = fill(1, n)
    while true
        noadvance = true
        for itask in 1:n
            basin = basinsets[itask]
            Q, dists = protobasinsets[itask]
            start = indices[itask]
            if isempty(dists)
                start == length(Q) && continue
                stop = length(Q)
            else
                noadvance = false
                stop = popfirst!(dists) - 1
            end
            indices[itask] = stop + 1
            for i in start:stop
                u = Q[i]
                iu = CartesianIndex(mod1.(u, dims))
                iu in old && delete!(basin, u)
                push!(new, iu)
            end
        end
        noadvance && break
        union!(old, new)
        empty!(new)
    end
    basinsets
end

"""
    local_basins(grid::Array{<:Any,3}, minima::Vector{CartesianIndex{3}}; lt=(<), cluster=0, skip=Returns(false), maxdist=0, inplace=false)

Return a list `basinsets` whose elements are the local attraction basins around each energy
minimum `p` given in `minima`. Each element is a `Set` of grid points `(i,j,k)` such that
`p` is the closest local minimum for which there is a path that only decreases in energy
when going from `p` to `(i,j,k)`.
If `maxdist > 0`, the path length must be below `maxdist`.

The indices `i`, `j` and `k` may not be in the bounds of `grid` if the point is in a
periodic image of the unit cell.

For each pair of minima whose Manhattan distance on the grid is lower or equal to
`cluster`, their corresponding basins are merged into a single one.

Grid points whose energy `e` is such that `skip(e)` are not included in the basins. If such
a point is an element of `minima`, the corresponding basin will be skipped as well.

If `inplace` is set, the list `minima` will be filtered in-place by removing all indices
that should be `skip`ped.
"""
function local_basins(grid::Array{<:Any,3}, minima::Vector{CartesianIndex{3}};
                      lt=(<), cluster=0, skip=Returns(false), maxdist=0, inplace=false)
    nminima = length(minima)
    dims = size(grid)
    protobasinsets = Vector{Tuple{Vector{NTuple{3,Int}},Vector{Int}}}(undef, nminima)
    @threads for itask in 1:nminima
        start = minima[itask]
        visited = falses(dims)
        istart, jstart, kstart = Tuple(start)
        if skip(grid[istart, jstart, kstart])
            protobasinsets[itask] = (NTuple{3,Int}[], Int[]) # filtered below
            continue
        end
        visited[istart, jstart, kstart] = true
        Q = [Tuple(start)]
        current_dist = 0
        start_dist = 1
        dists = Int[1]
        mark_dist = false
        for (idx, u) in enumerate(Q)
            if idx == start_dist
                current_dist += 1
                mark_dist = true
            end
            iu = CartesianIndex(mod1.(u, dims))
            gu = grid[iu]
            for v in GridNeighbors(u)
                iv = CartesianIndex(mod1.(v, dims))
                visited[iv] && continue
                gv = grid[iv]
                skip(gv) && continue
                if (maxdist == 0 || current_dist < maxdist) && lt(gu, grid[iv])
                    push!(Q, v)
                    if mark_dist
                        start_dist = length(Q)
                        push!(dists, start_dist)
                        mark_dist = false
                    end
                    visited[iv] = true
                end
            end
        end
        protobasinsets[itask] = (Q, dists)
    end
    emptysets1 = findall(isempty∘first, protobasinsets)
    deleteat!(protobasinsets, emptysets1)
    newminima = deleteat!(inplace || isempty(emptysets1) ? minima : copy(minima), emptysets1)

    basinsets = keep_shortest_distance!(protobasinsets, dims)
    emptysets2 = findall(isempty, basinsets)
    deleteat!(basinsets, emptysets2)
    deleteat!(newminima, emptysets2)

    # clean_cluster!(basinsets, newminima, cluster, dims)
    # basinsets = [Set([Tuple(m)]) for m in minima]
    clean_cluster!(basinsets, minima, cluster, dims)
    basinsets
end

"""
    local_basins(grid::Array{<:Any,3}, tolerance::Float64=-Inf; lt=(<), cluster=0, skip=Returns(false), maxdist=0)

Equivalent to [`local_basins(grid::Array{<:Any,3}, minima::Vector{CartesianIndex{3}}; lt=(<), smooth=0, skip=Returns(false))`](@ref)
where `minima` is computed from [`local_minima`](@ref) with the given `tolerance` and
keyword arguments.
"""
function local_basins(grid::Array{<:Any,3}, tolerance::Float64=-Inf;
                      lt=(<), cluster=0, skip=Returns(false), maxdist=0)
    local_basins(grid, local_minima(grid, tolerance; lt); lt, cluster, skip, maxdist, inplace=true)
end


"""
    local_components(grid::Array{T,3}; atol=0.0, rtol=1.0, ntol=Inf) where T

Return a list of `Vector{Tuple{NTuple{3,Int},T}`, each pair being made of a sublist
representing a connected component of `grid` point coordinates that have a non-zero value
and the corresponding sum of its values.

The returned list is sorted by decreasing sum of the value of its sublists. Any sublist of
total value below `atol` is removed. The list is has at most `ntol` elements and it is also
truncated to the last sublist such that the sum of the kept values are below the fraction
`rtol` of the total value of `grid`.
"""
function local_components(grid::Array{T,3}; atol=0.0, rtol=1.0, ntol=Inf) where T
    a, b, c = dims = size(grid)
    visited = falses(a, b, c)
    ret = Tuple{Vector{NTuple{3,Int}},T}[]
    for k in 1:c, j in 1:b, i in 1:a
        visited[i,j,k] && continue
        tot = grid[i,j,k]
        iszero(tot) && continue
        visited[i,j,k] = true
        Q = [(i,j,k)]
        for I in Q
            for J in GridNeighbors(I)
                iv, jv, kv = mod1.(J, dims)
                visited[iv,jv,kv] && continue
                val = grid[iv,jv,kv]
                iszero(val) && continue
                tot += val
                visited[iv,jv,kv] = true
                push!(Q, J)
            end
        end
        tot < atol && continue
        push!(ret, (Q, tot))
    end
    sort!(ret; by=last, rev=true)
    max_tot = rtol < 1 ? sum(last, grid)*rtol : Inf
    total = 0.0
    for (i, (_, subtot)) in enumerate(ret)
        if i > ntol || (total += subtot) > max_tot
            resize!(ret, i-1)
            break
        end
    end
    ret
end

"""
    closest_to_sites(grid::Array{T,3}, sites::Union{Vector{SVector{3,Float64}},Vector{NTuple{3,Float64}}}, mat) where T

Return `(Qs, weights, deleted_sites)` where:
- `Qs[i]` is the list of points closest to site `i`.
- `weights[i]` is the list of corresponding `grid` values.
- `deleted_sites` is the list of sites that were deleted compared to the input `sites`. The
  index `i` before refers to the kept sites.
"""
function closest_to_sites(grid::Array{T,3}, sites::Union{Vector{SVector{3,Float64}},Vector{NTuple{3,Float64}}}, mat) where T
    a, b, c = dims = size(grid)
    _, ortho, safemin = prepare_periodic_distance_computations(mat)
    safemin2 = safemin^2
    protoskin = norm(mat * SVector{3,Float64}(inv(a), inv(b), inv(c)))
    Qs = [NTuple{3,Int}[(1+round(Int, a*i), 1+round(Int, b*j), 1+round(Int, c*k))] for (i,j,k) in sites]
    deletesite = Int[]
    knownsites = Set{NTuple{3,Int}}()
    for (idxQ, Q) in enumerate(Qs)
        ijk = only(Q)
        if ijk in knownsites
            push!(deletesite, idxQ)
        else
            push!(knownsites, ijk)
        end
    end
    deleteat!(Qs, deletesite)
    deleteat!(sites, deletesite)
    n = length(sites)

    # First, record all grid points which are unambiguously closer to one site than to any
    # other because they are closer to the site than the half distance between that site
    # and all other sites.
    @threads for isite in 1:n
        site = sites[isite]
        Q = Qs[isite]
        todelete = Int[]
        bufferA = MVector{3,typeof(safemin)}(undef)
        bufferB = MVector{3,Float64}(undef)
        maxdist2 = Inf*oneunit(safemin^2) # squared half smallest distance between a site and its neighbor sites)
        for i2 in 1:n
            isite == i2 && continue
            bufferB .= site .- sites[i2]
            dsites2 = periodic_distance2!(bufferA, mat, ortho, safemin2, bufferB)
            maxdist2 = min(maxdist2, dsites2)
        end
        maxdist2 /= 4
        skin2 = (protoskin + sqrt(maxdist2))^2
        visited = falses(a, b, c)
        visited[mod1.(Q[1], dims)...] = true
        ρ = grid[mod1.(Q[1], dims)...]
        for I in Q
            for J in GridNeighbors(I)
                j1, j2, j3 = mod1.(J, dims)
                visited[j1,j2,j3] && continue
                visited[j1,j2,j3] = true
                bufferB .= site .- (J.-1)./dims
                mul!(bufferA, mat, bufferB)
                d2 = norm2(bufferA) # no periodic_distance! since the offset is known
                d2 ≥ skin2 && continue
                d2 ≥ maxdist2 && push!(todelete, length(Q)+1)
                push!(Q, J)
                ρ += grid[j1,j2,j3]
                ρ > 1 && break
            end
            ρ > 1 && break
        end
        deleteat!(Q, todelete)
    end

    allvisited = falses(a, b, c)
    counter = 0
    for Q in Qs
        for pos in Q
            i,j,k = mod1.(pos, dims)
            allvisited[i,j,k] && @show i,j,k
            @assert !allvisited[i,j,k]
            allvisited[i,j,k] = true
            counter += 1
        end
    end

    # Then, for all non-zero grid points that have not been attributed their closest site,
    # compute the distance between each of these points and all the sites and identify the
    # closest.

    tovisit = CartesianIndex{3}[] # grid points to visit
    sizehint!(tovisit, length(allvisited) - counter)
    for k in 1:c, j in 1:b, i in 1:a
        if !allvisited[i,j,k] && !iszero(grid[i,j,k])
            push!(tovisit, CartesianIndex{3}(i,j,k))
        end
    end
    m = length(tovisit)
    alldists = fill((Inf*oneunit(safemin), zero(SVector{3,Int})), m, n)

    @threads for isite in 1:n
        site = sites[isite]
        bufferX = MVector{3,Float64}(undef)
        ofs = MVector{3,Int}(undef)
        for ivisit in 1:m
            ijk = Tuple(tovisit[ivisit]) ./ dims
            bufferX .= site .- ijk
            d = periodic_distance_with_ofs!(bufferX, ofs, mat, ortho, safemin)
            ofs .*= dims
            alldists[ivisit, isite] = (d, SVector{3,Float64}(ofs))
        end
    end
    tosites = findmin(alldists; dims=2)[2]
    weights = [[grid[mod1(i, a), mod1(j, b), mod1(k, c)] for (i,j,k) in Q] for Q in Qs]
    sumweights = sum.(weights)
    for jvisit in 1:m
        jsite = Tuple(tosites[jvisit])[2]
        s = sumweights[jsite]
        visit = tovisit[jvisit]
        ρ = grid[visit]
        s + ρ > 1 && continue
        sumweights[jsite] = s + ρ
        push!(weights[jsite], ρ)
        push!(Qs[jsite], Tuple(visit) .+ Tuple(alldists[jvisit,jsite][2]))
    end

    Qs, weights, deletesite
end

"""
    updated_sites_with_closest(Qs, weights, dims)

Takes arguments `Qs` and `weights` from [`closest_to_sites`](@ref) and return the new list
of sites along with their updated weight.

`dims` is `size(grid)` with `grid` the first argument to `closest_to_sites`.
"""
function updated_sites_with_closest(Qs, weights, dims)
    # Finally, return the list of pairs (p, center) where p is the sum of all values
    # closest to the site whose position is given by center (in fractional coordinates),
    # the barycenter of the positions of the belonging points weighted by their value.
    [begin
        s = sum(ws)
        iszero(s) && error("Encountered empty weight: please use a non-zero atol value.")
        (s, sum((SVector{3,Int}(ijk).-1)./dims.*w for (ijk, w) in zip(Q, ws))/s)
    end for (Q, ws) in zip(Qs, weights)]
end

"""
    site_proximity_map!(sites::AbstractVector, dims, mat::AbstractMatrix, maxdist=2.0*oneunit(eltype(mat)))

Return a `map` of the dimension `dims` such that `map[i,j,k] == (s, ofs)` where `s` is the
index of the site closest to point `(i,j,k)` if the distance is lower than `maxdist` or
`s == 0` otherwise, and `ofs` is the cell offset of the site with respect to point
`(i,j,k)`.
"""
function site_proximity_map!(sites::AbstractVector, dims, mat::AbstractMatrix, maxdist=2.0*oneunit(eltype(mat)))
    a, b, c = dims
    _, ortho, safemin = prepare_periodic_distance_computations(mat)
    safemin2 = safemin^2
    protoskin = norm(mat * SVector{3,Float64}(inv(a), inv(b), inv(c)))
    Qs = [NTuple{3,Int}[(1+round(Int, a*i), 1+round(Int, b*j), 1+round(Int, c*k))] for (i,j,k) in sites]
    nextQs = [Tuple{NTuple{3,Int},typeof(safemin2)}[] for _ in sites]
    deletesite = Int[]
    knownsites = Set{NTuple{3,Int}}()
    for (idxQ, Q) in enumerate(Qs)
        ijk = only(Q)
        if ijk in knownsites
            push!(deletesite, idxQ)
        else
            push!(knownsites, ijk)
        end
    end
    deleteat!(Qs, deletesite)
    deleteat!(sites, deletesite)
    n = length(sites)

    @threads for isite in 1:n
        site = sites[isite]
        Q = Qs[isite]
        nextQ = nextQs[isite]
        todelete = Int[]
        bufferA = MVector{3,typeof(safemin)}(undef)
        bufferB = MVector{3,Float64}(undef)
        maxdist2 = absmaxdist2 = maxdist^2
        for i2 in 1:n
            isite == i2 && continue
            bufferB .= site .- sites[i2]
            dsites2 = periodic_distance2!(bufferA, mat, ortho, safemin2, bufferB)
            maxdist2 = min(maxdist2, dsites2)
        end
        maxdist2 /= 4
        visited = falses(a, b, c)
        visited[mod1.(Q[1], dims)...] = true

        #= The two iterations of the following loop correspond to the two parts of the
        algorithm. In both cases, we loop on the points around each site by increasing
        distance. Since it is easy to explore grid neighbors, we loop until reaching a
        target distance plus an extra "skin" distance. The points between the target and
        the skin will removed after, but they are a necessary steps to reach other points
        that fall below the target.
        In the first iteration, the target distance is the half distance between the target
        site and its closest neighbor. As a consequence, all encountered points are
        guaranteed to have the target site as their closest site.
        In the second iteration, the target distance is the `maxdist` given as input. The
        points encountered in the second iteration are kept along their distance to the
        site to be later sorted by closest site.
        =#
        for (queue, nextqueue) in ((Q, nextQ), (nextQ, nothing))
            skin2 = (protoskin + sqrt(maxdist2))^2
            for _I in queue
                I = _I isa NTuple{3,Int} ? _I : (_I::Tuple{NTuple{3,Int},typeof(safemin2)})[1]
                for J in GridNeighbors(I)
                    j1, j2, j3 = mod1.(J, dims)
                    visited[j1,j2,j3] && continue
                    visited[j1,j2,j3] = true
                    bufferB .= site .- (J.-1)./dims
                    mul!(bufferA, mat, bufferB)
                    d2 = norm2(bufferA) # no periodic_distance! since the offset is known
                    if d2 ≥ maxdist2
                        if nextqueue isa Vector{Tuple{NTuple{3,Int},typeof(safemin2)}}
                            d2 < absmaxdist2 && push!(nextqueue, (J, d2))
                        end
                        d2 ≥ skin2 && continue
                        push!(todelete, length(queue)+1)
                    end
                    push!(queue, if queue isa Vector{NTuple{3,Int}}
                        J
                    else
                        (J, d2)
                    end)
                end
            end
            deleteat!(queue, todelete)
            empty!(todelete)
            maxdist2 = maxdist^2
        end
    end

    ret = fill((zero(Int32), zero(SVector{3,Int8})), a, b, c)
    dist2s = fill(Inf*oneunit(safemin2), a, b, c)

    # Handle the points encountered in the first iteration, which are guaranteed to have
    # their designated site as closest site
    for (s, Q) in enumerate(Qs)
        for pos in Q
            (oi, i), (oj, j), (ok, k) = fldmod1.(pos, dims)
            if ret[i,j,k][1] != 0
                @error "Assertion failure: point $((i,j,k)) seen in $(ret[i,j,k][1]) and $s"
            end
            ret[i,j,k] = (s, (oi-1, oj-1, ok-1))
            dist2s[i,j,k] = zero(safemin2) # do not overwrite this in the next loop
        end
    end

    # Handle the points encountered in the second iteration and only keep for each point
    # the closest site
    for (s, Q) in enumerate(nextQs)
        for (pos, d2) in Q
            (oi, i), (oj, j), (ok, k) = fldmod1.(pos, dims)
            if d2 < dist2s[i,j,k]
                dist2s[i,j,k] = d2
                ret[i,j,k] = (s, (oi-1, oj-1, ok-1))
            end
        end
    end

    ret
end

function closest_from_proximity(grid, proximity, numsites)
    ret = Vector{Tuple{Float64, SVector{3,Float64}}}(undef, numsites+1)
    a, b, c = size(grid)
    @assert (a, b, c) == size(proximity)
    @threads for isite in 0:numsites
        ρ = 0.0
        pos = zero(SVector{3,Int})
        for k in 1:c, j in 1:b, i in 1:a
            s, o = proximity[i,j,k]
            s == isite || continue
            val = grid[i,j,k]
            ρ += val
            pos += val.*(o .+ SVector{3}((i-1)/a, (j-1)/b, (k-1)/c))
        end
        ret[ifelse(iszero(isite), numsites+1, isite)] = (ρ, pos./ρ)
    end
    ret
end

"""
    proximal_map(grid, sites, mat, maxdist=2.0u"Å")

Return a list `l` such that `l[i]` is the total density closest to `sites[i]` and no
further than `maxdist`.

`l[end]` is the remaining non-attributed density.
"""
function proximal_map(grid, sites, mat, maxdist=2*oneunit(eltype(mat)))
    proximity = site_proximity_map!(sites, size(grid), mat, maxdist);
    first.(closest_from_proximity(grid, proximity, length(sites)))
end


function cluster_closer_than!(ret::Vector{Tuple{Float64,SVector{3,Float64}}}, mat, maxdist)
    n = length(ret)
    todelete = falses(length(ret))
    buffer, ortho, safemin = prepare_periodic_distance_computations(mat)
    ofs = MVector{3,Float64}(undef)
    for i1 in 1:n
        todelete[i1] && continue
        ρ1, pos1 = ret[i1]
        for i2 in (i1+1):n
            todelete[i2] && continue
            ρ2, pos2 = ret[i2]
            ρ = ρ1 + ρ2
            ρ > 1 && continue
            buffer .= pos1 .- pos2
            d = try
                periodic_distance_with_ofs!(buffer, ofs, mat, ortho, safemin)
            catch
                @show pos1, pos2
                rethrow()
            end
            d ≥ maxdist && continue
            todelete[i2] = true
            pos1 = (ρ1.*pos1 .+ ρ2.*(pos2.+ofs))./ρ
            ρ1 = ρ
            ret[i1] = (ρ1, pos1)
        end
    end
    deleteat!(ret, todelete)
    nothing
end

"""
    coalescing_basins(grid::Array{<:Any,3}, maxdist, mat; smooth=norm(mat*(inv.(SVector{3,Int}(size(grid))))), atol=0.1, crop=atol/1000, maxbasins=Inf)

Return basins found by a proximity algorithm looking at the values of the `grid` in
decreasing order.

The value of `smooth` is the width of the smoothing Gaussian applied to `grid` before
sorting the values to choose which ones to visit in decreasing order.
The default value corresponds to the diagonal of a voxel.

After smoothing, the grid is cropped so that any value lower than `crop` is considered
equal to zero when first identifying basin centers.

`atol` is the minimum density of any returned site.

`maxbasins` is a limit on the maximal number of basins explored.
"""
function coalescing_basins(grid::Array{<:Any,3}, maxdist, mat; smooth=nothing, atol=0.1, crop=nothing, maxbasins=Inf)
    a, b, c = dims = size(grid); ab = a*b
    ret = Tuple{Float64,SVector{3,Float64}}[] # total value and basin center
    smoothed = if smooth == 0
        grid
    else
        _smooth = smooth isa Nothing ? norm(mat*(inv.(SVector{3,Int}(size(grid))))) : smooth
        imfilter(grid, Kernel.gaussian(NoUnits.((_smooth.*dims)./norm.(eachcol(mat)))), "circular")
    end
    sortedidx = sortperm(vec(smoothed), rev=true)
    buffer, ortho, safemin = prepare_periodic_distance_computations(mat)
    ofs = MVector{3,Float64}(undef)
    in_distance = Tuple{Int,SVector{3,Float64},typeof(safemin)}[]
    js = Int[]
    newcenter = MVector{3,Float64}(undef)
    _crop = crop isa Nothing ? atol*maximum(grid)/100 : crop
    for newidx in sortedidx
        smval = smoothed[newidx]
        smval ≤ _crop && break
        val = grid[newidx]
        k, retk = fldmod1(newidx, ab)
        j, i = fldmod1(retk, a)
        pos = SVector{3,Float64}(i/a, j/b, k/c)
        empty!(in_distance)
        for (idx, (_, center)) in enumerate(ret)
            buffer .= pos .- center
            d = periodic_distance_with_ofs!(buffer, ofs, mat, ortho, safemin)
            d > maxdist && continue
            push!(in_distance, (idx, SVector{3,Float64}(ofs), d))
        end
        if isempty(in_distance) && smval ≥ atol
            push!(ret, (val, pos))
            length(ret) ≥ maxbasins && break
        else
            if length(in_distance) > 1
                sort!(in_distance; by=last)
                submid = val < atol ? 1 : something(findfirst(x -> last(x) < maxdist/2, in_distance), 1)
                resize!(in_distance, submid)
            end
            newcenter .= pos .* val
            newval = val
            for (j2, o2, _) in in_distance
                val2, center2 = ret[j2]
                newcenter .+= (center2.+o2).*val2
                newval += val2
            end
            iszero(newval) && continue
            newcenter ./= iszero(newval) ? 1.0 : newval
            if length(in_distance) == 1
                ret[in_distance[1][1]] = (newval, newcenter)
            else
                resize!(js, length(in_distance))
                js .= first.(in_distance)
                sort!(js)
                deleteat!(ret, js)
                push!(ret, (newval, newcenter))
            end
        end
    end
    n1 = -1
    while length(ret) != n1
        n1 = length(ret)
        Qs, weights, _ = closest_to_sites(grid, last.(ret), mat)
        ret = updated_sites_with_closest(Qs, weights, size(grid))
        cluster_closer_than!(ret, mat, maxdist)
        newret = filter(≥(atol)∘first, ret)
        isempty(newret) && return ret
        ret = newret
    end
    return ret
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
    decompose_basins(grid::Array{<:Any,3}, basinsets::Vector{Set{NTuple{3,Int}}})
    decompose_basins(grid::Array{<:Any,3}, tolerance::Float64=-Inf; lt=(<), smooth=0, skip=Returns(false))

Decompose the provided `basinsets` into a grid `nodes` where `nodes[i,j,k]` is the list of
basins to which the grid point `(i,j,k)` belongs to. Return `(points, nodes)` where
`points` is the list of grid point for which `nodes` is not empty.

If `basinsets` is not provided, it will be computed from [`local_basins`](@ref) with the
given `tolerance` and keyword arguments.
"""
function decompose_basins(grid::Array{<:Any,3}, basinsets::Vector{Set{NTuple{3,Int}}})
    nodes = fill(Int[], size(grid))
    a, b, c = size(nodes)
    n = length(basinsets)
    refvectors = [[i] for i in 1:n]
    refidx = Dict{Vector{Int},Int}([[i] => i for i in 1:n])
    pointsset = Set{NTuple{3,Int}}()
    for (iset, set) in enumerate(basinsets)
        union!(pointsset, set)
        for (i, j, k) in set
            u, v, w = mod1(i, a), mod1(j, b), mod1(k, c)
            node = nodes[u, v, w]
            push!(node, iset)
            ref = get(refidx, node, nothing)
            if ref isa Nothing
                newref = copy(node)
                push!(refvectors, newref)
                refidx[newref] = length(refvectors)
                ref = length(refvectors)
            end
            pop!(node)
            nodes[u, v, w] = refvectors[ref]
        end
    end
    points = collect(pointsset)
    sort!(points; lt=((i1,j1,k1),(i2,j2,k2)) -> (k1<k2 || (k1==k2&&(j1<j2 || (j1==j2 && (@assert(i1!=i2); i1<i2))))))
    points, nodes
end

function decompose_basins(grid::Array{<:Any,3}, tolerance::Float64=-Inf;
                          lt=(<), smooth=0, skip=Returns(false))
    decompose_basins(grid, local_basins(grid, tolerance; lt, smooth, skip))
end


"""
    energy_levels(egrid::Array{Float64,4}, n::Integer)
    energy_levels(egrid::Array{Float64,4}, step::AbstractFloat)
    energy_levels(egrid::Array{Float64,4}, steps::AbstractRange{<:AbstractFloat})

Given an angular energy grids, return the energies and relative volumes of the lowest
levels (between energy -∞ and 0 K).

If `n` is given, split the energy equally between `n` levels.
If `step` is given, split the energy between levels of size `step` until reaching 0 K.
If `steps` is given, split the energy at the given steps.
"""
function energy_levels(egrid::Array{Float64,4}, x::Union{Integer,AbstractFloat,AbstractRange{<:AbstractFloat}})
    kept = filter(<(0.0), vec(egrid))
    if x isa Integer
        n = x
        m = minimum(kept)
        λ = -m/(2n)
        rnge = (m+λ):(2λ):0.0
    elseif x isa AbstractFloat
        n = round(Int, cld(abs(m), x))
        m = minimum(kept)
        λ = x/2
        rnge = (m+λ):x:(m+(2n-1)*λ)
    else
        n = length(x)
        λ = step(x)/2
        m = first(x)
        rnge = x
    end
    ofsm = m - 2λ
    γ = n/abs(m)
    hist = zeros(Float64, n)
    for x in kept
        idx = floor(Int, (x-ofsm)*γ)
        hist[clamp(idx, 1, n)] += 1
    end
    hist ./= (length(egrid))
    rnge, hist
end


function compute_levels(grid::Array{Float64,4}, mine, maxe, T=300)
    basins = Vector{NTuple{3,Int}}[]
    _, a1, a2, a3 = size(grid)
    dims = (a1, a2, a3)
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
            for J in GridNeighbors(I)
                j1, j2, j3 = mod1.(J, dims)
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
