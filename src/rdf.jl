using Chemfiles
using Base.Threads

const CHEMFILES_LOCK = ReentrantLock()


function compute_angle(posA, posB)
    pos = posB - posA
    p = pos/norm(pos)
    θ = acos(p[3])
    1 - abs(p[3]) < 1e-13 && return (θ, 0.0)
    q = SVector{2,Float64}(pos[1], pos[2])
    q /= norm(q)
    ϕ = acos(q[1])
    if q[2] < 0
        ϕ = 2π - ϕ
    end
    return θ, ϕ
end

function get_molecules(frame, mat; skipmolid)
    invmat = inv(mat)
    orig_poss = [SVector{3,Cdouble}(p) for p in eachcol(positions(frame))]
    poss = [invmat*x for x in orig_poss]
    buffer, ortho, safemin = prepare_periodic_distance_computations(mat)
    n = length(poss)
    groupmap = collect(1:n) # groupmap[i] is the id of the group to which atom i belongs
    groups = [[i] for i in 1:n] # groups[k] is the list of atoms l in group k, i.e. such that groupmap[l] == k
    offsets = zeros(SVector{3,Int}, n) # offsets[i] is the offset of atom i in molecule groupmap[i]
    ofs = zero(MVector{3,Int})
    for i in 1:n, j in (i+1):n
        buffer .= poss[i] .- poss[j]
        d = periodic_distance_with_ofs!(buffer, ofs, mat, ortho, safemin)
        @assert d > 0.7
        if 0.7 < d < 1.5
            ki = groupmap[i]
            kj = groupmap[j]
            ofs .+= offsets[i]
            for l in groups[kj]
                groupmap[l] = ki
                offsets[l] += ofs
            end
            append!(groups[ki], groups[kj])
            empty!(groups[kj])
        end
    end
    for (i, o) in enumerate(offsets)
        poss[i] += o
    end
    filter!(!isempty, groups)
    # from this point, groupmap is invalid. groups[l] is the list of atom numbers of molecule l
    m = length(groups)
    moleculesID = Vector{Int}(undef, m) # For each molecule, an integer that identifies its formula (e.g. 1 for CO2, 2 for O2, 3 for N2)
    IDmap = Vector{Int}[] # for each molecule ID, the list of individual molecules with this composition
    IDdict = Dict{Vector{Int},Int}() # map a list of atom identifier to a molecule ID
    key = Int[]
    moleculepos = Vector{SVector{3,Float64}}(undef, m) # in fractional coordinates
    moleculeangle = Vector{NTuple{2,Float64}}(undef, m)
    pos = MVector{3,Float64}(undef)
    for (i, mol) in enumerate(groups)
        pos[1] = pos[2] = pos[3] = 0
        numatoms = length(mol)
        resize!(key, numatoms)
        for (j, atom) in enumerate(mol)
            pos .+= poss[atom]
            key[j] = skipmolid ? 0 : Chemfiles.atomic_number(view(frame, j))
        end
        moleculepos[i] = SVector{3,Float64}(pos ./ numatoms)
        θ, ϕ = compute_angle(orig_poss[mol[1]], orig_poss[mol[2]]) # TODO: only works for linear symmetric molecules
        if ϕ > π
            ϕ -= π
            θ = π - θ
        end
        moleculeangle[i] = (θ, ϕ)
        molID = get!(IDdict, sort!(key), length(IDmap)+1)
        if molID > length(IDmap)
            push!(IDmap, copy(key))
        end
    end
    moleculepos, moleculeangle, moleculesID, IDmap, groups
end

function parse_rdf(file, maxdist=12.0, nbins=500; skipmolid=true)
    trajectory = Trajectory(file)
    threadbins = zeros(nthreads(), nbins)
    mat = SMatrix{3,3,Cdouble,9}(matrix(UnitCell(read_step(trajectory, 0))))
    ϵ = maxdist / nbins
    @threads for i_frame in 1:Int(length(trajectory))
        buffer, ortho, safemin = prepare_periodic_distance_computations(mat)
        frame = @lock CHEMFILES_LOCK read_step(trajectory, i_frame-1)
        moleculepos, moleculeangle, moleculesID, IDmap, groups = get_molecules(frame, mat; skipmolid)
        n = length(moleculepos)
        for i in 1:n, j in (i+1):n
            buffer .= moleculepos[i] .- moleculepos[j]
            idx = 1 + floor(Int, periodic_distance!(buffer, mat, ortho, safemin)/ϵ)
            idx > nbins && continue
            threadbins[threadid(), idx] += 1
        end
    end
    bins = dropdims(sum(threadbins; dims=1); dims=1)
    for i in 1:nbins
        bins[i] /= 4π*((i*ϵ)*(1 + i*ϵ)) # δV between sphere of radius i*ϵ and (i+1)*ϵ
    end
    bins
end
