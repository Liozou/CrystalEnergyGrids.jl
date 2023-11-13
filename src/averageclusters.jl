using ImageFiltering: Kernel, mapwindow
using StatsBase: median, AnalyticWeights

struct PDBModelIterator
    io::IOStream
end
PDBModelIterator(path::AbstractString) = PDBModelIterator(open(path))
Base.IteratorSize(::Type{PDBModelIterator}) = Base.SizeUnknown()
Base.eltype(::Type{PDBModelIterator}) = Tuple{Int,SMatrix{3,3,TÅ,9},Vector{Tuple{Int,SubString{String},SVector{3,TÅ}}}}
Base.isdone(stream::PDBModelIterator, _=nothing) = eof(stream.io)
function Base.iterate(stream::PDBModelIterator, _=nothing)
    record = Tuple{Int,SubString{String},SVector{3,TÅ}}[]
    modelid = 0
    mat = SMatrix{3,3,Float64,9}(NaN,0,0,0,NaN,0,0,0,NaN)
    while true
        if eof(stream.io)
            close(stream.io)
            isempty(record) && return nothing
            return (modelid, stream.mat, record), nothing
        end
        l = readline(stream.io)
        if startswith(l, "MODEL")
            modelid = parse(Int, @view l[7:end])
        elseif startswith(l, "CRYST1")
            a, b, c, α, β, γ = parse.(Float64, l[rge] for rge in (7:15, 16:24, 25:33, 34:40, 41:47, 48:54))::Vector{Float64}
            mat = mat_from_parameters((a, b, c), (α, β, γ))
        elseif startswith(l, "ATOM")
            num = parse(Int, @view l[7:11])
            el = @view l[77:78]
            pos = SVector{3,Float64}(parse(Float64, @view l[31:38]), parse(Float64, @view l[39:46]), parse(Float64, @view l[47:54]))
            push!(record, (num, el, pos*u"Å"))
        elseif startswith(l, "ENDMDL")
            return (modelid, mat, record), nothing
        end
    end
end


"""
    mat_from_output(path::AbstractString, ::Val{EXT}) where EXT

Given the `path` to an output file, return the unit cell extracted from it.
`EXT` is either `:pdb` for .pdb files, or `:stream` for `.stream` or `.zst` output.
"""
function mat_from_output(path::AbstractString, ::Val{EXT}) where EXT
    if EXT === :pdb
        open(path) do io
            l = readline(io)
            while !startswith(l, "CRYST1")
                l = readline(io)
            end
            a, b, c, α, β, γ = parse.(Float64, l[rge] for rge in (7:15, 16:24, 25:33, 34:40, 41:47, 48:54))::Vector{Float64}
            mat_from_parameters((a, b, c), (α, β, γ)).*u"Å"
        end
    else
        @assert EXT === :stream
        stream = StreamSimulationStep(path)
        mat = first(stream).mat
        close(stream)
        mat
    end
end

function find_output_file(path)
    if isfile(path)
        basename(path) == "steps.stream" && return path, :stream, Val(:stream)
        basename(path) == "steps.zst" && return path, :zst, Val(:stream)
        basename(path) == "trajectory.pdb" && return path, :pdb, Val(:pdb)
    else
        isdir(path) || error(lazy"Input $path is neither a directory nor a file.")
        file = joinpath(path, "steps.stream")
        isfile(file) && return file, :stream, Val(:stream)
        file = joinpath(path, "steps.zst")
        isfile(file) && return file, :zst, Val(:stream)
        file = joinpath(path, "trajectory.pdb")
        isfile(file) && return file, :pdb, Val(:pdb)
    end
    return "", :none, Val(:stream)
end

"""
    bin_trajectory(path; step=0.15u"Å", bins=nothing)

Given the `path` to a directory containing the outputs of a simulation, or directly to the
output file, return a pair `(bins, total)` such that `bins` is a grid of specified `step`
size, where each point is mapped to the number of atoms found in that position, and `total`
is the number of frames of the output.

As a consequence, `bins ./ total` is a map of the average atomic density.

A preallocated `bins` array can be given if it has the correct size. It will be filled by
the appropriate values and returned instead of a new array.
"""
function bin_trajectory(path; step=0.15u"Å", bins=nothing)
    output, kind, valkind = find_output_file(path)
    kind === :none && error(lazy"No available output among \"steps.stream\", \"steps.zst\" or \"trajectory.pdb\" at $dir.")
    _bin_trajectory(output, valkind; step, _bins=bins)
end

function _bin_trajectory(path, ::Val{EXT}; step, _bins) where EXT
    @assert EXT === :pdb || EXT === :stream
    steps = (EXT === :pdb ? PDBModelIterator : StreamSimulationStep)(path)
    mat = mat_from_output(path, Val(EXT))
    invmat = inv(NoUnits.(mat./u"Å"))
    if _bins isa Nothing
        na, nb, nc = floor.(Int, NoUnits.(cell_lengths(mat) ./ step))
        bins = zeros(Int, na, nb, nc)
    else
        na, nb, nc = size(_bins)
        bins = _bins
        bins .= 0
    end

    stepcounter = 0
    for step in steps
        stepcounter += 1
        for x in (EXT === :stream ? step.positions : step[3])
            pos = EXT === :stream ? x : x[3]
            y = invmat*NoUnits.(pos./u"Å")
            i, j, k = y .- floor.(y)
            bins[ceil(Int, i*na)+(i==0), ceil(Int, j*nb)+(j==0), ceil(Int, k*nc)+(k==0)] += 1
        end
    end
    bins, stepcounter
end

function bin_trajectories(path; step=0.15u"Å", except=())
    dirs = filter(x -> isdir(joinpath(path, x)) && x ∉ except, readdir(path; join=false))
    isempty(dirs) && error(lazy"Directory $path does not contain any valid folder")
    output, kind, valkind = find_output_file(joinpath(path, first(dirs)))
    kind === :none && error(lazy"Could not find any suitable output at $(joinpath(path, first(dirs)))")
    mat = mat_from_output(output, valkind)
    na, nb, nc = floor.(Int, NoUnits.(cell_lengths(mat) ./ step))
    bins = zeros(Int, na, nb, nc)
    totalbins = zeros(Int, na, nb, nc)
    totalcounter = 0
    for dir in dirs
        _, counter = bin_trajectory(joinpath(path, dir); step, bins)
        totalbins .+= bins
        totalcounter += counter
        yield()
    end
    totalbins, totalcounter
end


"""
    average_clusters(bins::Array{Float64,3}; clustering=3, smoothing=3, threshold=0.01)

Takes a grid `bins` of atomic densities and group them into sites, each site corresponding
to a local maximum in density. If two sites are closer than `clustering` (in terms of
Manhattan distance on the grid), they are merged. Sites whose total density is lower than
`threshold` are discared.

Return a list of pairs `(p, center)` where `p` is the probability of having an atom on site
at position `center`. The coordinates of `center` are in [0, 1).
"""
function average_clusters(bins::Array{Float64,3}; clustering=3, smoothing=5, threshold=0.01)
    if smoothing > 0
        # grid = imfilter(bins, Kernel.gaussian((smoothing, smoothing, smoothing)), "circular")
        weights = [inv(1+sqrt((i-(smoothing+1)/2)^2+(j-(smoothing+1)/2)^2+(k-(smoothing+1)/2)^2)) for i in 1:smoothing, j in 1:smoothing, k in 1:smoothing]
        grid = mapwindow(bins, (smoothing, smoothing, smoothing), "circular") do buf
            median(vec(buf), weights)
        end
        grid .*= sum(bins)/sum(grid)
        # skip = ≤(maximum(grid[I] for I in eachindex(bins) if iszero(bins[I])))
    else
        grid = bins
    end
    basinsets = if smoothing > 0
        local_basins(grid, -Inf; lt=(>), clustering, maxdist=0, skip=Returns(false))
    else
        local_basins(grid, -Inf; lt=(>), clustering, maxdist=0, skip=iszero)
    end
    _average_clusters(grid, basinsets, threshold)
end

function _average_clusters(grid, basinsets, threshold)
    n = length(basinsets)
    @assert !any(isempty, basinsets)
    points, nodes = decompose_basins(grid, basinsets)
    _average_clusters(grid, points, nodes, n, threshold)
end

function _average_clusters(grid, points, nodes, n, threshold)
    a, b, c = size(nodes)
    probabilities = zeros(n)
    centers = [zero(SVector{3,Float64}) for _ in 1:n]
    for (i,j,k) in points
        u, v, w = mod1(i, a), mod1(j, b), mod1(k, c)
        belongs = nodes[u,v,w]
        @assert !isempty(belongs)
        density = grid[u,v,w]/length(belongs)
        for belong in belongs
            probabilities[belong] += density
            centers[belong] += density*SVector{3,Float64}(i,j,k) # not u,v,w, keep the offset!
        end
    end
    ret = Tuple{Float64,SVector{3,Float64}}[]
    for (p, center) in zip(probabilities, centers)
        p > threshold || continue
        normalized_c = center ./ (p*a, p*b, p*c)
        push!(ret, (p, normalized_c .- floor.(Int, normalized_c)))
    end
    ret
end

"""
    average_clusters(dir; smoothing=3, step=0.15u"Å", threshold=0.01)

Given the path to a folder `dir` containing the outputs of a simulation, group atomic
densities into sites.

See [`bin_trajectory`](@ref) for the meaning of `step` and
[`average_clusters`](@ref) for the return type and the meaning of `smoothing` and
`threshold`.
"""
function average_clusters(dir; smoothing=3, step=0.15u"Å", threshold=0.01)
    average_clusters(bin_trajectory(dir; step)[1]; smoothing, threshold)
end


"""
    output_sites(path::AbstractString, framework_path::AbstractString, sites, atomsymbol::Union{AbstractString,Symbol})

Output to `path` a .cif file representing the framework (whose .cif file is stored at
`framework_path`) as well as a number of `sites` with their respective occupancies, each
corresponding to a species of symbol `atomsymbol`.

The `sites` can be obtained from [`average_clusters`](@ref).
"""
function output_sites(path::AbstractString, framework_path::AbstractString, sites, atomsymbol::Union{AbstractString,Symbol})
    if isdir(path) || isdirpath(path)
        name = string(splitext(basename(framework_path))[1], '_', atomsymbol)
        path = joinpath(path, name*".cif")
    else
        name = basename(path)
        if splitext(path)[2] != ".cif"
            path = path*".cif"
        end
    end
    mkpath(dirname(path))

    open(framework_path) do input; open(path, "w") do output
        inloop = inatomloop = loopstarted = false
        counter = 0
        order = Symbol[]
        hasoccupancy = false
        for l in eachline(input)
            if isempty(l)
                if inloop && inatomloop && loopstarted
                    for (p, (x, y, z)) in sites
                        counter += 1
                        for o in order
                            if o === :?
                                print(output, "? ")
                            elseif o === :label
                                print(output, atomsymbol, '_', counter, '\t')
                            elseif o === :symbol
                                print(output, atomsymbol, '\t')
                            elseif o === :x
                                @printf output "%12g" x
                            elseif o === :y
                                @printf output "%12g" y
                            elseif o === :z
                                @printf output "%12g" z
                            elseif o === :occ
                                @printf output "%12g" min(p, 1.0)
                            end
                        end
                        println(output)
                    end
                end
                inloop = false
                inatomloop = loopstarted = false
            elseif l == "loop_"
                inloop = true
                empty!(order)
                inatomloop = loopstarted = false
            elseif inloop
                if !loopstarted && l[1] != '_'
                    loopstarted = true
                    counter = 0
                    if inatomloop && !hasoccupancy
                        println(output, "_atom_site_occupancy")
                        push!(order, :occ)
                    end
                end
                if loopstarted
                    counter += 1
                    if inatomloop && !hasoccupancy
                        l = l * " 1.0"
                    end
                else
                    if l == "_atom_site_occupancy"
                        hasoccupancy = true
                        push!(order, :occ)
                    elseif l == "_atom_site_fract_x"
                        inatomloop = true
                        push!(order, :x)
                    elseif l == "_atom_site_fract_y"
                        push!(order, :y)
                    elseif l == "_atom_site_fract_z"
                        push!(order, :z)
                    elseif l == "_atom_site_label"
                        push!(order, :label)
                    elseif l == "_atom_site_type_symbol"
                        push!(order, :symbol)
                    else
                        push!(order, :?)
                    end
                end
            end
            println(output, l)
        end
    end end # open(input) and open(output)
end
