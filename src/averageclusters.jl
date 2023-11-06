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

function find_output_file(dir)
    path = joinpath(dir, "steps.stream")
    isfile(path) && return path, :stream, Val(:stream)
    path = joinpath(dir, "steps.zst")
    isfile(path) && return path, :zst, Val(:stream)
    path = joinpath(dir, "trajectory.pdb")
    isfile(path) && return path, :pdb, Val(:pdb)
    return "", :none, Val(:stream)
end

"""
    bin_trajectory(dir; refT::TK=300.0u"K", step=0.15u"Å")

Given the path to a directory `dir` containing the outputs of a simulation,
return a grid of specified `step` where each point is mapped to the Boltzmann average (at
temperature `refT`) of the atomic density.

Also return `(mine, λ)` where `mine` is the minimum of `energies` and `λ` is the sum of
`exp((mine-e)/refT)` for `e` in `energies`.
"""
function bin_trajectory(dir; refT::TK=300.0u"K", step=0.15u"Å", mat=nothing, bins=nothing)
    epath = joinpath(dir, "energies.serial")
    isfile(epath) || error(lazy"Cannot bin trajectory: missing file $epath.")
    energies = Float64.(deserialize(epath)::Vector{BaselineEnergyReport})*u"K"
    output, kind, valkind = find_output_file(dir)
    kind === :none && error(lazy"No available output among \"steps.stream\", \"steps.zst\" or \"trajectory.pdb\" at $dir.")
    _bin_trajectory(output, energies, valkind; refT, step, _mat=mat, _bins=bins)
end

function _bin_trajectory(path, energies::Vector{TK}, ::Val{EXT}; refT, step, _mat, _bins) where EXT
    @assert EXT === :pdb || EXT === :stream
    steps = (EXT === :pdb ? PDBModelIterator : StreamSimulationStep)(path)
    mat = _mat isa Nothing ? mat_from_output(path, Val(EXT)) : _mat
    invmat = inv(NoUnits.(mat./u"Å"))
    if _bins isa Nothing
        na, nb, nc = floor.(Int, NoUnits.(cell_lengths(mat) ./ step))
        bins = zeros(Float64, na, nb, nc)
    else
        na, nb, nc = size(_bins)
        bins = _bins
    end
    mine::TK = minimum(energies)

    stepcounter = 0
    for step in steps
        stepcounter += 1
        energy = energies[stepcounter]
        for x in (EXT === :stream ? step.positions : step[3])
            pos = EXT === :stream ? x : x[3]
            y = invmat*NoUnits.(pos./u"Å")
            i, j, k = y .- floor.(y)
            val = exp((mine-energy)/refT)
            bins[ceil(Int, i*na)+(i==0), ceil(Int, j*nb)+(j==0), ceil(Int, k*nc)+(k==0)] += val
        end
    end
    @assert stepcounter == length(energies)
    λ = sum(exp((mine-e)/refT) for e in energies; init=0.0)
    bins ./= λ
    bins, (mine, λ)
end

function bin_trajectories(path; refT::TK=300.0u"K", step=0.15u"Å", except=())
    dirs = filter(∉(except), readdir(path; join=false))
    isempty(dirs) && error(lazy"Directory $path does not contain any valid folder")
    output, kind, valkind = find_output_file(joinpath(path, first(dirs)))
    kind === :none && error(lazy"Could not find any suitable output at $(joinpath(path, first(dirs)))")
    mat = mat_from_output(output, valkind)
    na, nb, nc = floor.(Int, NoUnits.(cell_lengths(mat) ./ step))
    bins = zeros(Float64, na, nb, nc)
    totalbins = zeros(Float64, na, nb, nc)
    totalmine = Inf*u"K"
    μ = 0.0
    for dir in dirs
        _, (mine, λ) =  bin_trajectory(joinpath(path, dir); refT, step, mat, bins)
        if mine < totalmine
            α = exp((mine - totalmine)/refT)
            β = λ
            totalmine = mine
        else
            α = 1.0
            β = λ*exp((totalmine - mine)/refT)
        end
        μ = α*μ + β
        for k in 1:nc, j in 1:nb, i in 1:na
            totalbins[i,j,k] = α*totalbins[i,j,k] + β*bins[i,j,k]
        end
        bins .= 0.0
        yield()
    end
    totalbins ./= μ
    totalbins, (totalmine, μ)
end

function _add_cyclic((x, y, z), (i, j, k), (a, b, c), λ=true, μ=true)
    ma, mb, mc = (a÷2)*μ, (b÷2)*μ, (c÷2)*μ
    u = x ≤ ma && i > ma ? (i-a) : x > ma && i ≤ ma ? (i+a) : i
    v = y ≤ mb && j > mb ? (j-b) : y > mb && j ≤ mb ? (j+b) : j
    w = z ≤ mc && k > mc ? (k-c) : z > mc && k ≤ mc ? (k+c) : k
    SVector{3,Float64}(λ*u, λ*v, λ*w)
end

"""
    average_clusters(bins::Array{Float64,3}; smoothing=3, threshold=0.01)

Takes a grid `bins` of atomic densities and group them into sites, each site corresponding
to a local maximum in density. If two sites are closer than `smoothing` (in terms of
Manhattan distance on the grid), they are merged. Sites whose total density is lower than
`threshold` are discared.

Return a list of pairs `(p, center)` where `p` is the probability of having an atom on site
at position `center`. The coordinates of `center` are in [0, 1).
"""
function average_clusters(bins::Array{Float64,3}, smoothing=3, threshold=0.01)
    if smoothing > 0
        grid = imfilter(bins, Kernel.gaussian((smoothing, smoothing, smoothing)), "circular")
        # skip = ≤(maximum(grid[I] for I in eachindex(bins) if iszero(bins[I])))
    else
        grid = bins
    end
    basinsets = local_basins(grid, -Inf; lt=(>), smoothing, skip=iszero)
    n = length(basinsets)
    @assert !any(isempty, basinsets)
    points, nodes = decompose_basins(grid, basinsets)
    probabilities = zeros(n)
    center = [zero(SVector{3,Float64}) for _ in 1:n]
    numnodes = zeros(Float64, n)
    for (i,j,k) in points
        belongs = nodes[i,j,k]
        @assert !isempty(belongs)
        λ = inv(length(belongs))
        for b in belongs
            probabilities[b] += grid[i,j,k]*λ
            center[b] += _add_cyclic(center[b], SVector{3,Float64}(i,j,k), size(grid), λ, numnodes[b])
            numnodes[b] += λ
        end
    end
    ret = Tuple{Float64,SVector{3,Float64}}[]
    for (p, c, λ) in zip(probabilities, center, numnodes)
        p > threshold || continue
        normalized_c = c ./ λ
        push!(ret, (p, normalized_c .- floor.(Int, normalized_c)))
    end
    ret
end

"""
    average_clusters(dir; smoothing=3, refT::TK=300.0u"K", step=0.15u"Å", threshold=0.01)

Given the path to a folder `dir` containing the outputs of a simulation, group atomic
densities into sites.

See [`bin_trajectory`](@ref) for the meaning of `step` and `refT` and
[`average_clusters`](@ref) for the return type and the meaning of `smoothing` and
`threshold`.
"""
function average_clusters(dir; smoothing=3, refT::TK=300.0u"K", step=0.15u"Å", threshold=0.01)
    average_clusters(bin_trajectory(dir; refT, step)[1]; smoothing, threshold)
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
