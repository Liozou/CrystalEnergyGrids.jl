# Computation of energy resulting from a system with a given force field

export SimulationStep, make_step, update_position, update_position!

"""
    SimulationStep{N}

Represent a set of systems in space along with a force-field, i.e. a single frame in time.

To compute the energy of such a step, see [`energy_nocutoff`](@ref) and [`compute_vdw`](@ref).

The systems are factored into kinds: one kind may designate a molecule of water, another a
zeolite framework. There can be any number of systems for each kind: simulating the wetting
of a zeolite will for instance require one system of the "zeolite framework" kind and many
systems of kind "water".
"""
struct SimulationStep{N}
    ff::ForceField
    charges::Vector{Te_au}
    # charges[ix] is the charge of the atom of index ix in ff, i.e. the charge of the k-th
    # atom in a system of kind i is charges[ffidx[i][k]].
    positions::Vector{SVector{N,TÅ}} # one position per atom
    atoms::Vector{Tuple{Int,Int,Int}} # one index per atom
    posidx::Vector{Vector{Vector{Int}}}
    # the position of the k-th atom in the j-th system of kind i is positions[l] where
    # atoms[l] == (i,j,k) and posidx[i][j][k] == l
    freespecies::Vector{Vector{Int}}
    # freespecies[i] is a list of indices j that indicate that posidx[i][j] are not to be
    # taken into account. For example, if there used to be 3 species of kind i, and the 2nd
    # was removed, then 2 is added to the list freespecies[i]. If a new species of kind i
    # is added, the last element j of freespecies[i] is removed from the list and the
    # positions of the new species are given by posidx[i][j]. If freespecies[i] is empty,
    # take j = length(posidx[i]) + 1.
    isrigid::BitVector # isrigid[i] applies to all systems of kind i
    ffidx::Vector{Vector{Int}}
    # ffidx[i][k] is ff.sdict[atomic_symbol(systemkinds[i], k)], i.e. the numeric index
    # in the force field for the k-th atom in a system of kind i
    cell::CellMatrix
end

function SimulationStep(ff::ForceField, charges::Vector{Te_au},
                        inputpos::Vector{Vector{Vector{SVector{N,TÅ}}}},
                        isrigid::BitVector, ffidx::Vector{Vector{Int}}, cell::CellMatrix) where N

    numatoms = sum(x -> sum(length, x; init=0), inputpos; init=0)
    positions = Vector{SVector{N,TÅ}}(undef, numatoms)
    atoms = Vector{Tuple{Int,Int,Int}}(undef, numatoms)
    posidx = Vector{Vector{Vector{Int}}}(undef, length(inputpos))
    l = 0
    for (i, posi) in enumerate(inputpos)
        posidxi = Vector{Vector{Int}}(undef, length(posi))
        posidx[i] = posidxi
        for (j, poss) in enumerate(posi)
            posidxi[j] = collect(l+1:l+length(poss))
            for (k, pos) in enumerate(poss)
                l += 1
                atoms[l] = (i,j,k)
                positions[l] = pos
            end
        end
    end
    freespecies = [Int[] for _ in inputpos]
    SimulationStep(ff, charges, positions, atoms, posidx, freespecies, isrigid, ffidx, cell)
end

"""
    SimulationStep(ff::ForceField, systemkinds::Vector{T} where T<:AbstractSystem,
                   inputpos::Vector{Vector{Vector{SVector{N,TÅ}}}},
                   cell::CellMatrix,
                   isrigid::BitVector=trues(length(systemkinds))) where N

Create a `SimulationStep` from a force field `ff`, a list of system kinds, the cell matrix
and the positions of all atoms for all systems.

`systemkinds[i]` is a system that represents all systems of kind `i`. This means that all
its properties are shared with other systems of kind `i`, except the positions of its atoms.

`inputpos[i][j]` is the list of positions of all atoms in the `j`-th system of kind `i`.
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
                        inputpos::Vector{Vector{Vector{SVector{N,TÅ}}}},
                        cell::CellMatrix,
                        isrigid::BitVector=trues(length(systemkinds))) where N
    @assert length(systemkinds) == length(inputpos) == length(isrigid)
    ffidx = [[ff.sdict[atomic_symbol(s, k)] for k in 1:length(s)] for s in systemkinds]
    charges = [[uconvert(u"e_au", s[k,:atomic_charge])::Te_au for k in 1:length(s)] for s in systemkinds]
    SimulationStep(ff, charges, inputpos, isrigid, ffidx, cell)
end

"""
    make_step(ff::ForceField, systems::Vector{T}, cell::CellMatrix=CellMatrix(), isrigid::BitVector=trues(length(systems))) where T<:AbstractSystem

Create a [`SimulationStep`](@ref) from a given force field `ff` and a list of interacting
`systems`. `isrigid[i]` specifies whether `systems[i]` is considered rigid. If unspecified,
all systems are considered rigid (so only inter-system interactions are taken into account).

Return a tuple `(step, indices)` where `step` is the `SimulationStep` and `indices[i]` is
the index of `systems[i]` in `step`, for use with [`update_position`](@ref) and such.

Systems sharing the same ordered list of atom symbols and the same rigidity will be
considered as having the same system kind (e.g. several molecules of water).

If unspecified, `cell` is set to an open (i.e. aperiodic) box.
"""
function make_step(ff::ForceField, systems::Vector{T}, cell::CellMatrix=CellMatrix(),
                   isrigid::BitVector=trues(length(systems))) where T<:AbstractSystem
    length(isrigid) > length(systems) && error("Less systems provided than `isrigid` specifications")
    if length(isrigid) < length(systems)
        @info "All systems beyond the $(length(isrigid)) first ones are implicitly considered rigid."
        append!(isrigid, true for _ in length(systems)-length(isrigid))
    end
    kindsdict = Dict{Tuple{Vector{Symbol},Bool},Int}()
    systemkinds = T[]
    U = Vector{SVector{n_dimensions(systems[1]),TÅ}} # positions of the atoms of a system
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
    SimulationStep(ff, systemkinds, poss, cell, newisrigid), indices
end

"""
    SimulationStep(step::SimulationStep, mode=:all)

Copy of the input on which the positions and number of species can be modified without
affecting the original `step`.

If `mode === :output`, only the `positions` and `atoms` fields are copied.
If `mode === :zero`, create a `SimulationStep` with an empty `ffidx` field.
"""
function SimulationStep(step::SimulationStep, mode=:all)
    if mode === :all
        SimulationStep(step.ff, step.charges, copy(step.positions), copy(step.atoms),
                       [[copy(js) for js in is] for is in step.posidx],
                       [copy(x) for x in step.freespecies], step.isrigid, step.ffidx, step.cell)
    elseif mode === :output
        SimulationStep(step.ff, step.charges, copy(step.positions), copy(step.atoms),
                       step.posidx, step.freespecies, step.isrigid, step.ffidx, step.cell)
    elseif mode === :zero
        SimulationStep(step.ff, step.charges, step.positions, step.atoms, step.posidx,
                       step.freespecies, step.isrigid, Vector{Int}[], step.cell)
    else
        error("Please use either :all, :output or :zero as value for argument mode")
    end
end

"""
    update_position!(step::SimulationStep, idx, newpos)

Modifies the position of the system of index `idx` in `step` to `newpos`.
Return `newpos`.

See also [`update_position!(step::SimulationStep, idx, op, arg)`](@ref) and [`update_position`](@ref).
"""
function update_position!(step::SimulationStep, (i,j), newpos)
    molpos = step.posidx[i][j]
    for (k, pos) in enumerate(newpos)
        step.positions[molpos[k]] = pos
    end
end

"""
    update_position(step::SimulationStep, idx, newpos)

Return a new `SimulationStep` where the system of index `idx` has a new position `newpos`.

`idx` can be the index returned by [`make_step`](@ref) or a tuple `(i,j)` which designates
the `j`-th system of kind `i` in `step`.

See also [`update_position!`](@ref) to modify `step` in-place.
"""
function update_position(step::SimulationStep, (i,j), newpos)
    newpositions = copy(step.positions)
    x = SimulationStep(step.ff, step.charges, newpositions, step.atoms, step.posidx, step.freespecies, step.isrigid, step.ffidx, step.cell)
    update_position!(x, (i,j), newpos)
    x
end

# function energy_intra(step::SimulationStep, i::Int, positions::Vector{SVector{N,TÅ}}) where N
#     n = length(positions)
#     energy = 0.0u"K"
#     cutoff2 = step.ff.cutoff^2
#     isinf(step.ff.cutoff) || error("finite cutoff not implemented")
#     indices = step.ffidx[i]
#     for k1 in 1:n
#         ix1 = indices[k1]
#         for k2 in (i+1):n
#             ix2 = indices[k2]
#             r2 = norm2(positions[k1], positions[k2])
#             r2 > cutoff2 && error("Flexible molecule cannot be larger than the cutoff")
#             rule = step.ff[ix1, ix2]
#             energy += rule(r2)# + tailcorrection(rule, step.ff.cutoff)
#         end
#     end
#     energy
# end

# function energy_nocutoff(step::SimulationStep)
#     energy = 0.0u"K"
#     nkinds = length(step.ffidx)
#     for i1 in 1:nkinds
#         poskind1 = step.positions[i1]
#         idx1 = step.ffidx[i1]
#         rigid1 = step.isrigid[i1]
#         for (j1, pos1) in enumerate(poskind1)
#             rigid1 || (energy += energy_intra(step, i1, pos1))
#             for i2 in i1:nkinds
#                 poskind2 = step.positions[i2]
#                 idx2 = step.ffidx[i2]
#                 for j2 in (j1*(i1==i2)+1):length(poskind2)
#                     pos2 = poskind2[j2]
#                     for (k1, p1) in enumerate(pos1), (k2, p2) in enumerate(pos2)
#                         energy += step.ff[idx1[k1], idx2[k2]](norm2(p1, p2))
#                     end
#                 end
#             end
#         end
#     end
#     energy
# end

function compute_vdw(step::SimulationStep)
    energy = 0.0u"K"
    cutoff2 = step.ff.cutoff^2
    buffer = MVector{3,TÅ}(undef)
    buffer2 = MVector{3,Float64}(undef)
    for (l1, pos1) in enumerate(step.positions)
        i1, j1, k1 = step.atoms[l1]
        ix1 = step.ffidx[i1][k1]
        for l2 in (l1+1):length(step.positions)
            i2, j2, k2 = step.atoms[l2]
            i1 == i2 && j1 == j2 && continue # no intra energy considered here
            pos2 = step.positions[l2]
            buffer .= pos2 .- pos1
            d2 = unsafe_periodic_distance2!(buffer, buffer2, step.cell)
            if d2 < cutoff2
                energy += step.ff[ix1, step.ffidx[i2][k2]](d2)
            end
        end
    end
    for (i, rigid) in enumerate(step.isrigid)
        rigid && continue
        posidxi = step.posidx[i]
        for molpos in posidxi
            energy += energy_intra(step, i, @view step.positions[molpos])
        end
    end
    energy
end

function single_contribution_vdw(step::SimulationStep, (i1,j1)::Tuple{Int,Int}, poss1::AbstractVector{SVector{N,TÅ}}) where N
    energy = step.isrigid[i1] ? 0.0u"K" : energy_intra(step, i1, poss1)
    cutoff2 = step.ff.cutoff^2
    buffer = MVector{3,TÅ}(undef)
    buffer2 = MVector{3,Float64}(undef)
    idx1 = step.ffidx[i1]
    for (k1, pos1) in enumerate(poss1)
        ix1 = idx1[k1]
        for (l2, pos2) in enumerate(step.positions)
            i2, j2, k2 = step.atoms[l2]
            i1 == i2 && j1 == j2 && continue
            buffer .= pos2 .- pos1
            d2 = unsafe_periodic_distance2!(buffer, buffer2, step.cell)
            if d2 < cutoff2
                energy += step.ff[ix1, step.ffidx[i2][k2]](d2)
            end
        end
    end
    energy
end
