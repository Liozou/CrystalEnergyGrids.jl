# Computation of energy resulting from a system with a given force field

export SimulationStep, make_step, update_position, update_position!

"""
    SimulationStep{N,T<:AbstractSystem{N}}

Represent a set of systems in space along with a force-field, i.e. a single frame in time.

To compute the energy of such a step, see [`energy_nocutoff`](@ref).

The systems are factored into kinds: one kind may designate a molecule of water, another a
zeolite framework. There can be any number of systems for each kind: simulating the wetting
of a zeolite will for instance require one system of the "zeolite framework" kind and many
systems of kind "water".
"""
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

"""
    SimulationStep(ff::ForceField, systemkinds::Vector{T} where T<:AbstractSystem,
                   positions::Vector{Vector{Vector{SVector{N,typeof(1.0u"Å")}}}} where N,
                   isrigid::BitVector=trues(length(systemkinds)))

Create a `SimulationStep` from a force field `ff`, a list of system kinds, and the
positions of all atoms for all systems.

`systemkinds[i]` is a system that represents all systems of kind `i`. This means that all
its properties are shared with other systems of kind `i`, except the positions of its atoms.

`positions[i][j]` is the list of positions of all atoms in the `j`-th system of kind `i`.
The `k`-th element of such a list should correspond to the same atom as the `k`-th one of
`systemkinds[i]`.

`isrigid[i]` specifies whether system `i` is considered rigid, i.e. whether interactions
between atoms within the structure should be considerd (flexible) or not (rigid). If
unspecified, all systems are considered rigid (so only inter-system interactions are taken
into account).

See also [`make_step`](@ref) for an alternative way of constructing a `SimulationStep`
without having to specify the system kinds.
"""
function SimulationStep(ff::ForceField, systemkinds::Vector{T} where T<:AbstractSystem,
                        positions::Vector{Vector{Vector{SVector{N,typeof(1.0u"Å")}}}} where N,
                        isrigid::BitVector=trues(length(systemkinds)))
    @assert length(systemkinds) == length(positions) == length(isrigid)
    idx = [[ff.sdict[atomic_symbol(s, k)] for k in 1:length(s)] for s in systemkinds]
    SimulationStep(ff, systemkinds, positions, isrigid, idx)
end

"""
    make_step(ff::ForceField, systems::Vector{T}, isrigid::BitVector=trues(length(systems))) where T<:AbstractSystem

Create a [`SimulationStep`](@ref) from a given force field `ff` and a list of interacting
`systems`. `isrigid[i]` specifies whether `systems[i]` is considered rigid. If unspecified,
all systems are considered rigid (so only inter-system interactions are taken into account).

Return a tuple `(step, indices)` where `step` is the `SimulationStep` and `indices[i]` is
the index of `systems[i]` in `step`, for use with [`update_position`](@ref) and such.

Systems sharing the same ordered list of atom symbols and the same rigidity will be
considered as having the same system kind (e.g. several molecules of water).
"""
function make_step(ff::ForceField, systems::Vector{T}, isrigid::BitVector=trues(length(systems))) where T<:AbstractSystem
    length(isrigid) > length(systems) && error("Less systems provided than `isrigid` specifications")
    if length(isrigid) < length(systems)
        @info "All systems beyond the $(length(isrigid)) first ones are implicitly considered rigid."
        append!(isrigid, true for _ in length(systems)-length(isrigid))
    end
    kindsdict = Dict{Tuple{Vector{Symbol},Bool},Int}()
    systemkinds = T[]
    U = Vector{SVector{n_dimensions(systems[1]),typeof(1.0u"Å")}} # positions of the atoms of a system
    poss = Vector{U}[]
    newisrigid = BitVector()
    indices = Tuple{Int,Int}[]
    for (system, thisrigid) in zip(systems, isrigid)
        n = length(kindsdict)+1
        kind = get!(kindsdict, (atomic_symbol(system), thisrigid), n)
        if kind === n
            push!(systemkinds, system)
            push!(poss, U[])
            push!(newisrigid, thisrigid)
        end
        push!(poss[kind], position(system))
        push!(indices, (kind, length(poss[kind])))
    end
    SimulationStep(ff, systemkinds, poss, newisrigid), indices
end

"""
    update_position!(step::SimulationStep, idx, newpos)

Modifies the position of the system of index `idx` in `step` to `newpos`

See also [`update_position`](@ref)
"""
function update_position!(step::SimulationStep, (i,j), newpos)
    step.positions[i][j] = newpos
    nothing
end

"""
    update_position!(step::SimulationStep, idx, op, arg)

Modifies the position `p` of system of index `idx` in `step` into `op.(p, arg)`

See also [`update_position`](@ref)
"""
function update_position!(step::SimulationStep, (i,j), op, arg)
    step.positions[i][j] .= op.(step.positions[i][j], arg)
    nothing
end

"""
    update_position(step::SimulationStep, idx, newpos)

Return a new `SimulationStep` where the system of index `idx` has a new position `newpos`.

`idx` can be the index returned by [`make_step`](@ref) or a tuple `(i,j)` which designated
the `j`-th system of kind `i` in `step`.

See also [`update_position!`](@ref) to modify `step` in-place.
"""
function update_position(step::SimulationStep, (i,j), args...)
    newpositions = copy(step.positions)
    newpositions[i] = copy(step.positions[i])
    length(args) == 2 && (newpositions[i][j] = copy(newpositions[i][j]))
    x = SimulationStep(step.ff, step.systemkinds, newpositions, step.isrigid, step.idx)
    update_position!(x, (i,j), args...)
    x
end

"""
    update_position(step::SimulationStep, idx, op, arg)

Return a new `SimulationStep` where the system of index `idx` has a new position
`op.(p, arg)` where `p` is the current position.

See also [`update_position!`](@ref), [`unalias_position`](@ref)
"""
update_position

"""
    unalias_position(step, idx)

Return a new `SimulationStep` identical to `step` except that the system of index `idx` can
have its position mutated by [`update_position!`](@ref) without modifying `step`.
"""
function unalias_position(step, (i,j))
    newpositions = copy(step.positions)
    newpositions[i] = copy(step.positions[i])
    newpositions[i][j] = copy(newpositions[i][j])
    SimulationStep(step.ff, step.systemkinds, newpositions, step.isrigid, step.idx)
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
    isinf(ff.cutoff) || error("finite cutoff not implemented")
    for i in 1:n, j in (i+1):n
        r2 = norm2(positions[i], positions[j])
        rule = ff[symbols[i], symbols[j]]
        if r2 < ff.cutoff2
            energy += rule(r2) + tailcorrection(rule, ff.cutoff)
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

# function energy_withcutoff(step::SimulationStep)
#     energy = 0.0u"K"
#     nkinds = length(step.systemkinds)
#     for (i1, kind1) in enumerate(step.systemkinds)
#         poskind1 = step.positions[i1]
#         idx1 = step.idx[i1]
#         rigid1 = step.isrigid[i1]
#         for (j1, pos1) in enumerate(poskind1)
#             rigid1 || (energy += energy_intra(ff, ChangePositionSystem(kind1, pos1)))
#             for i2 in i1:nkinds
#                 poskind2 = step.positions[i2]
#                 idx2 = step.idx[i2]
#                 for j2 in (j1*(i1==i2)+1):length(poskind2)
#                     pos2 = poskind2[j2]
#                     for (k1, p1) in enumerate(pos1), (k2, p2) in enumerate(pos2)
#                         energy += step.ff(idx1[k1], idx2[k2], norm2(p1, p2))
#                     end
#                 end
#             end
#         end
#     end
#     energy
# end
