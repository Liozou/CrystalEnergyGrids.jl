# Computation of energy resulting from a system with a given force field

export SimulationStep

struct SimulationStep{N,T<:AbstractSystem{N}}
    ff::ForceField
    systemkinds::Vector{T}
    positions::Vector{Vector{Vector{SVector{N,typeof(1.0u"Å")}}}}
    # positions[i][j][k] is the position of the k-th atom in the j-th system of kind i
    isrigid::BitVector # isrigid[i] applies to all systems of kind i
    idx::Vector{Vector{Int}}
    # idx[i][k] is ff.sdict[atomic_symbol(systemkinds[i], k)], i.e. the numeric index
    # in the force field for the k-th atom in a system of kind i
end
function SimulationStep(ff::ForceField, systemkinds::Vector{T} where T<:AbstractSystem,
                        positions::Vector{Vector{Vector{SVector{N,typeof(1.0u"Å")}}}} where N,
                        isrigid::BitVector=trues(length(systemkinds)))
    idx = [[ff.sdict[atomic_symbol(s, k)] for k in 1:length(s)] for s in systemkinds]
    SimulationStep(ff, systemkinds, positions, isrigid, idx)
end
function SimulationStep(ff::ForceField, systems::Vector{T},
                        isrigid::BitVector=trues(length(systems))) where T<:AbstractSystem
    length(isrigid) > length(system) && error("Less systems provided than `isrigid` specifications")
    if length(isrigid) < length(systems)
        @info "All systems beyond the $(length(isrigid)) first ones are implicitly considered rigid."
        append!(isrigid, true for _ in length(systems)-length(isrigid))
    end
    kindsdict = Dict{Tuple{Vector{Symbol},Bool},Int}()
    systemkinds = T[]
    U = Vector{SVector{N,typeof(1.0u"Å")}} # positions of the atoms of a system
    positions = Vector{U}[]
    newisrigid = BitVector()
    for (system, thisrigid) in zip(systems, isrigid)
        n = length(kindsdict)+1
        kind = get!(kindsdict, (atomic_symbol(systems), thisrigid), n)
        if kind === n
            push!(systemkinds, system)
            push!(positions, U[])
            push!(newisrigid, thisrigid)
        end
        push!(positions[kind], positions(system))
    end
    SimulationStep(ff, systemkinds, positions, isrigid)
end

function norm2(u::S, v::T) where {S,T}
    r2 = zero(eltype(S))*zero(eltype(T)) # == zero(u)^2 == zero(v)^2
    for (x, y) in zip(u, v)
        r2 += (x - y)^2
    end
    r2
end

function energy_intra(ff::ForceField, system::AbstractSystem)
    positions = position(system)
    symbols = atomic_symbol(system)
    n = length(positions)
    energy = 0.0u"K"
    for i in 1:n, j in (i+1):n
        r2 = norm2(positions[i], positions[j])
        if isinf(ff.cutoff)
            energy += ff[symbols[i], symbols[j]](r2)
        else
            energy += ff(symbols[i], symbols[j], r2)
            error("not implemented")
        end
    end
    energy
end

function energy_nocutoff(step::SimulationStep)
    energy = 0.0u"K"
    nkinds = length(step.systemkinds)
    for (i1, kind1) in enumerate(step.systemkinds)
        poskind1 = step.positions[i1]
        idx1 = step.idx[i1]
        rigid1 = step.isrigid[i1]
        for (j1, pos1) in enumerate(poskind1)
            rigid1 || (energy += energy_intra(ff, ChangePositionSystem(kind1, pos1)))
            for i2 in i1:nkinds
                poskind2 = step.positions[i2]
                idx2 = step.idx[i2]
                for j2 in (j1*(i1==i2)+1):length(poskind2)
                    pos2 = poskind2[j2]
                    for (k1, p1) in enumerate(pos1), (k2, p2) in enumerate(pos2)
                        energy += step.ff[idx1[k1], idx2[k2]](norm2(p1, p2))
                    end
                end
            end
        end
    end
    energy
end
