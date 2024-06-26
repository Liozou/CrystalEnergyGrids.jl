# using ImageFiltering: mapwindow
# using StatsBase: median, AnalyticWeights

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
    mat_from_output(path::AbstractString, ::Val{EXT}, supercell=(1,1,1)) where EXT

Given the `path` to an output file, return the unit cell extracted from it.
`EXT` is either `:pdb` for .pdb files, or `:stream` for `.stream` or `.zst` output.

`supercell` is triplet of factors representing the extracted matrix. The lengths will be
divided by these values to retrieve the actual unit cell.
"""
function mat_from_output(path::AbstractString, ::Val{EXT}, supercell=(1,1,1)) where EXT
    if EXT === :pdb
        open(path) do io
            l = readline(io)
            while !startswith(l, "CRYST1")
                l = readline(io)
            end
            a, b, c, α, β, γ = parse.(Float64, l[rge] for rge in (7:15, 16:24, 25:33, 34:40, 41:47, 48:54))::Vector{Float64}
            mat_from_parameters((a, b, c)./supercell, (α, β, γ)).*u"Å"
        end
    else
        @assert EXT === :stream
        stream = StreamSimulationStep(path)
        mat = try
            first(stream).mat
        catch
            @error "Could not read from stream at $path"
            rethrow()
        end
        close(stream)
        hcat((SVector{3,TÅ}.(eachcol(mat))./supercell)...)::SMatrix{3,3,TÅ,9}
    end
end

function mat_from_output(path; supercell=(1,1,1))
    output, kind, valkind = find_output_file(path)
    kind === :none && error(lazy"No available output among \"steps.stream\", \"steps.zst\" or \"trajectory.pdb\" at $path")
    mat_from_output(output, valkind, supercell)
end

function find_output_file(path)
    if isfile(path)
        ext = splitext(path)[2]
        ext == ".stream" && return path, :stream, Val(:stream)
        ext == ".zst" && return path, :zst, Val(:stream)
        ext == ".pdb" && return path, :pdb, Val(:pdb)
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
    bin_trajectory(path::AbstractString; step=0.15u"Å", skip=0, count=typemax(Int), bins=nothing, symmetries=[], supercell=(1,1,1))

Given the `path` to a directory containing the outputs of a simulation, or directly to the
output file, return a pair `(bins, total)` such that `bins` is a grid of specified `step`
size, where each point is mapped to the number of atoms found in that position, and `total`
is the number of frames of the output.

As a consequence, `bins ./ total` is a map of the average atomic density.

The first `skip` frames are discarded. Up to `count` frames are kept.

A preallocated `bins` array can be given if it has the correct size. The bins corresponding
to `path` will be added on the values existing in the preallocated array. Calling
`bin_trajectory` without a preallocated `bins` or with one filled with zeros is equivalent.

`symmetries` can be given as a list of `EquivalentPositions`: these operations will be
applied to each encountered atom position before binning.
`total` is multiplied by `1 + length(symmetries)` so that `bins ./ total` still represents
the density map.

Similarly, `supercell` can be provided to signify that the actual unit cell is smaller than
the one from the file, and binning should be done on the unit cell. Each position will then
be mapped to its image in the unit cell before binning. This operation happens before
applying the `symmetries` (which are taken with respect to the unit cell).
`total` is multiplied by `prod(supercell)`  so that `bins ./ total` still represents the
density map.
"""
function bin_trajectory(path::AbstractString; step=0.15u"Å", skip=0, count=typemax(Int), bins=nothing, symmetries=[], supercell=(1,1,1))
    output, kind, valkind = find_output_file(path)
    kind === :none && error(lazy"No available output among \"steps.stream\", \"steps.zst\" or \"trajectory.pdb\" at $path.")
    _bin_trajectory(output, bins, valkind; step, skip, count, symmetries, supercell)
end

function _pathofsavedgrid(path::AbstractString; skip, count, symmetries)
    file, kind, _ = find_output_file(path)
    suff = string(count==typemax(Int) ? "" : "-c$count", skip==0 ? "" : "-s$skip")
    preff = string("grid_", isempty(symmetries) ? "raw_" : "symm_")
    if kind === :none
        while path != "/" && last(path) == '/'
            path = path[1:end-1]
        end
        joinpath(path, string(preff, basename(path), suff))
    else
        d = dirname(file)
        bd = basename(d)
        if startswith(bd, "System_") && basename(dirname(d)) == "Movies"
            # In RASPA folder
            sitesdir = joinpath(dirname(dirname(d)), "sites")
            if !isdir(sitesdir)
                mkdir(sitesdir)
            end
            joinpath(sitesdir, string(preff, first(splitext(basename(file))), suff))
        else
            joinpath(d, string(preff, bd, suff))
        end
    end
end

function _bin_trajectory(path, _bins, ::Val{EXT}; step, skip, count, symmetries, supercell) where EXT
    @assert EXT === :pdb || EXT === :stream
    @assert skip ≥ 0 && count ≥ 0
    pathofsave = _pathofsavedgrid(path; skip, count, symmetries)
    if isfile(pathofsave)
        @info "Retrieved saved grid from $pathofsave"
        savedbins, savedtotal = deserialize(pathofsave)::Tuple{Array{Int,3},Int}
        _bins isa Nothing && return (savedbins, savedtotal)
        _bins .+= savedbins
        return _bins, savedtotal
    end
    steps = (EXT === :pdb ? PDBModelIterator : StreamSimulationStep)(path)
    mat = mat_from_output(path, Val(EXT), supercell)
    invmat = inv(ustrip.(u"Å", mat))
    if _bins isa Nothing
        na, nb, nc = floor.(Int, NoUnits.(cell_lengths(mat) ./ step))
        bins = zeros(Int, na, nb, nc)
    else
        na, nb, nc = size(_bins)
        bins = _bins
    end

    stepcounter = 0
    for step in steps
        stepcounter += 1
        stepcounter ≤ skip && continue
        for x in (EXT === :stream ? step.positions : step[3])
            pos = EXT === :stream ? x : x[3]
            y = invmat*ustrip.(u"Å", pos)
            u, v, w = y .- floor.(y)
            bins[ceil(Int, u*na)+(u==0), ceil(Int, v*nb)+(v==0), ceil(Int, w*nc)+(w==0)] += 1
            for eq in symmetries
                eqy = eq(y)
                equ, eqv, eqw = eqy .- floor.(eqy)
                bins[ceil(Int, equ*na)+(equ==0), ceil(Int, eqv*nb)+(eqv==0), ceil(Int, eqw*nc)+(eqw==0)] += 1
            end
        end
        if (stepcounter - skip) ≥ count
            EXT === :stream && close(steps)
            break
        end
    end
    total = (stepcounter > skip)*(1 + length(symmetries))*prod(supercell)*(stepcounter - skip)
    if _bins isa Nothing
        @info "Saved grid at $pathofsave."
        serialize(pathofsave, (bins, total))
    end
    return bins, total
end

"""
    bin_trajectories(dirs; except=(), step=0.15u"Å", skip=0, count=typemax(Int), symmetries=[], supercell=(1,1,1))

Recursively call `bin_trajectory` on the hierarchy of folder starting with those in `dirs`.
Return the global `(bins, total)` where each value is the sum of those obtained on each
subdirectoy.

Directories whose names are in `except` are skipped, as well as all their subfolders.
For the other keyword arguments, see [`bin_trajectory`](@ref).
"""
function bin_trajectories(dirs; except=(), step=0.15u"Å", skip=0, count=typemax(Int), symmetries=[], supercell=(1,1,1))
    bins = nothing
    total = 0
    for dir in dirs
        basename(dir) in except && continue
        bins, total = _bin_trajectories(bins, total, dir; except, step, skip, count, symmetries, supercell)
    end
    bins, total
end

function _bin_trajectories(bins, total, dir; except, kwargs...)
    output, kind, valkind = find_output_file(dir)
    if kind === :none
        for x in readdir(dir; join=true)
            basename(x) in except && continue
            isdir(x) || continue
            bins, total = _bin_trajectories(bins, total, x; except, kwargs...)
        end
    else
        bins, addtototal = _bin_trajectory(output, bins, valkind; kwargs...)
        total += addtototal
        yield()
    end
    bins, total
end

"""
    bin_trajectories(path::AbstractString; except=(), step=0.15u"Å", skip=0, count=typemax(Int), symmetries=[], supercell=(1,1,1))

Call [`bin_trajectories`](@ref) `path` and all its subfolders.
"""
function bin_trajectories(path::AbstractString; step=0.15u"Å", skip=0, count=typemax(Int), except=(), symmetries=[], supercell=(1,1,1))
    pathofsave = _pathofsavedgrid(path; skip, count, symmetries)
    if isfile(pathofsave)
        @info "Retrieved saved grid from $pathofsave"
        return deserialize(pathofsave)::Tuple{Array{Int,3},Int}
    end
    ret = bin_trajectories((path,); except, step, skip, count, symmetries, supercell)
    @info "Saved grid at $pathofsave."
    serialize(pathofsave, ret)
    return ret
end

"""
    average_clusters(bins::Array{Float64,3}, dist, mat; smooth=0, crop=sum(bins)/(1000length(bins)), atol=0.0, rtol=1.0, ntol=Inf)

Takes a grid `bins` of atomic densities and group them into sites, each site corresponding
to a local density cluster. `dist` is the minimal separation distance between two sites,
`mat` is the unit cell matrix.

The grid is smoothed by a Gaussian of variance `smooth`. Use `smooth=0` to skip it. The
default is to use the diagonal of a voxel as variance.
Any value of `bins` lower than `crop` (after smoothing if so) is replaced by `0.0`. The
default is computed depending on `atol`.

There are different ways to remove extraneous sites:
- Any site whose density is lower than `atol` is removed.
- Keeping the sites by decreasing order of density, once the ratio of the kept density on
  the total density is below `rtol`, all remaining sites are removed. In other words,
  `rtol=1` does not remove any site, `rtol=0` only keeps the single most populated site,
  and any value in between may remove some sites, low-density first.
- Keeping the sites by decreasing order of density, no more than `ntol` sites are kept.

Return a list of pairs `(p, center)` where `p` is the probability of having an atom on site
at position `center`. The coordinates of `center` are in [0, 1).
"""
function average_clusters(bins::Array{Float64,3}, dist, mat; smooth=nothing, crop=nothing, atol=0.0, rtol=1, ntol=Inf)
    ret = coalescing_basins(bins, dist, mat; smooth, atol, crop, maxbasins=2*ntol)
    sort!(ret; rev=true)
    rtol == 1 && return ret
    max_total = rtol * sum(bins)
    total = 0.0
    for (i, (p, _)) in enumerate(ret)
        total += p
        if total ≥ max_total || i+1 > ntol
            resize!(ret, i)
            break
        end
    end
    return ret
end


"""
    average_clusters(dir, dist; step=0.15u"Å", symmetries=[], smooth=nothing, crop=nothing, atol=0.0, rtol=1, ntol=Inf, skip=0, count=typemax(Int), supercell=(1,1,1))

Given the path to a folder `dir` containing the outputs of a simulation, group atomic
densities into sites.

See [`bin_trajectory`](@ref) for the meaning of `step`, `symmetries`, `skip`, `count` and
`supercell`.
See [`average_clusters`](@ref) for the return type and the meaning of the other arguments.
"""
function average_clusters(dir, dist; step=0.15u"Å", symmetries=[], smooth=nothing, crop=nothing, atol=0.0, rtol=1, ntol=Inf, skip=0, count=typemax(Int), supercell=(1,1,1))
    bins, counter = bin_trajectory(dir; step, symmetries, skip, count, supercell)
    mat = mat_from_output(dir; supercell)
    average_clusters(bins./counter, dist, mat; smooth, crop, atol, rtol, ntol)
end

function _print_sites(io, counter, sites, order, atomsymbol)
    for (p, (x, y, z)) in sites
        counter += 1
        for o in order
            if o === :?
                print(io, "? ")
            elseif o === :label
                print(io, atomsymbol, '_', counter, '\t')
            elseif o === :symbol
                print(io, atomsymbol, '\t')
            elseif o === :x
                @printf io "%12g" x
            elseif o === :y
                @printf io "%12g" y
            elseif o === :z
                @printf io "%12g" z
            elseif o === :occ
                @printf io "%12g" min(p, 1.0)
            end
        end
        println(io)
    end
    counter
end


"""
    output_sites(path::AbstractString, framework_path::AbstractString, sites, atomsymbol::Union{AbstractString,Symbol}, keepsource=true)

Output to `path` a .cif file representing the framework (whose .cif file is stored at
`framework_path`) as well as a number of `sites` with their respective occupancies, each
corresponding to a species of symbol `atomsymbol`.

Unset `keepsource` to only output the sites (not the atoms of the framework).

The `sites` can be obtained from [`average_clusters`](@ref).
"""
function output_sites(path::AbstractString, framework_path::AbstractString, sites, atomsymbol::Union{AbstractString,Symbol}, keepsource=true)
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
        write_symmetries = false
        symmetries_written = false
        skiploop = false
        writtensites = false
        for l in eachline(input)
            if isempty(l)
                skiploop = false
                if inloop && inatomloop && loopstarted && !writtensites
                    writtensites = true
                    counter = _print_sites(output, counter, sites, order, atomsymbol)
                elseif write_symmetries && !symmetries_written
                    symmetries_written = true
                    write_symmetries = false
                    println(output)
                    println(output, "_symmetry_space_group_name_H-M	'P 1'")
                    println(output, "_symmetry_Int_Tables_number	1")
                    println(output, "_symmetry_cell_setting	triclinic")
                end
                inloop = false
                inatomloop = loopstarted = false
            elseif l == "loop_"
                if inloop && inatomloop && loopstarted && !writtensites
                    writtensites = true
                    counter = _print_sites(output, counter, sites, order, atomsymbol)
                end
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
                    skiploop && continue
                    inatomloop && !keepsource && continue # do not print framework atoms
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
                    elseif l == "_symmetry_equiv_pos_as_xyz" || l == "_space_group_symop_operation_xyz"
                        println(output, l)
                        println(output, "x,y,z")
                        skiploop = true
                        continue
                    else
                        push!(order, :?)
                    end
                end
            elseif startswith(l, "_symmetry")
                write_symmetries = true
                continue
            end
            println(output, l)
        end
        inloop && inatomloop && loopstarted && !writtensites && _print_sites(output, counter, sites, order, atomsymbol)
    end end # open(input) and open(output)
    nothing
end

function output_sites(path::AbstractString, framework_path::AbstractString, atomsymbol::Union{AbstractString,Symbol}; dist=1.0u"Å", step=0.15u"Å", symmetries=[], smooth=nothing, crop=nothing, atol=0.0, rtol=1, ntol=Inf, except=(), skip=0, count=typemax(Int), supercell=(1,1,1))
    isdir(path) || error("Expected a directory, given $path")
    sites = if isdir(joinpath(path, "1"))
        bins, counter = bin_trajectories(path; symmetries, step, except, skip, count, supercell)
        mat = mat_from_output(joinpath(path, "1"); supercell)
        average_clusters(bins./counter, dist, mat; smooth, crop, atol, rtol, ntol)
    else
        average_clusters(path, dist; step, symmetries, smooth, crop, atol, rtol, ntol, skip, count, supercell)
    end
    output_sites(joinpath(path, "sites.cif"), framework_path, sites, atomsymbol)
    sites
end
