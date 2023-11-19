# Computation of energy resulting from a system with a given force field

export SimulationStep, make_step, update_position, update_position!
import LinearAlgebra

"""
    SimulationStep{N,T}

Represent a set of systems in space along with a force-field, i.e. a single frame in time.

To compute the energy of such a step, see [`energy_nocutoff`](@ref) and [`compute_vdw`](@ref).


## Nomenclature

The systems are factored into kinds: one kind may designate a molecule of water, another a
zeolite framework. There can be any number of systems for each kind: simulating the wetting
of a zeolite will for instance require one system of the "zeolite framework" kind and many
systems of kind "water".

Each kind of system is referred by a number `i`.
For each such kind, each instance of the kind (i.e. each molecule) is referred by another
number `j`.
Finally, for each kind, each atom is referred by a last number `k`.

As a consequence, each atom of the `SimulationStep` is uniquely determined by its
corresponding triplet `(i,j,k)`.

The force field associates properties to each atom according to its category. The category
is usually referred by a name in the force field definition (like "Ow" to designate the
oxygens of water molecule), and is also internally represented by a number `ix`.
Each pair `(i,k)` designating a specific atom in a given species kind corresponds to a
unique category `ix` in the forcefield.


## Fields (part of the API)
- `ff`: the [`ForceField`](@ref).
- `charges`: `charges[ix]` is the charge of an atom of category `ix`.
- `mat`: the 3×3 matrix representing the unit cell.
- `positions`: a `Vector` of `SVector{N,typeof(1.0u"Å")}` representing the coordinates of
  all the atoms in the system.
- `parallel`: true (by default) if computations on the step should be parallelized.
- `atoms`: the list of triplet `(i,j,k)` defining the atoms, in the same order than their
  corresponding coordinates in `positions`.
- `posidx`: the converse mapping, i.e. `atoms[l] == (i,j,k)` implies `posidx[i][j][k] == l`.
- `freespecies`: `freespecies[i]` is a list of indices `j` that indicate that `posidx[i][j]`
  is an invalid sublist, not to be taken into account.
- `ffidx`: a `Vector{Vector{Int}}` such that `ffidx[i][k] == ix` where `ix` is the category
  of the `k`-th atom of a species of kind `i` in the force field.

Note that `atoms[l] == (i,j,k)` implies that `posidx[i][j][k] == l` but the converse does
not hold: `posidx[i][j][k]` may be defined even though there is no atom `l` such that
`atoms[l] == (i,j,k)`. The value of `posidx[i][j][k]` is then meaningless.
In particular, the list of valid indices `j` for kind `i` is
`setdiff(eachindex(posidx[i]), freespecies[i])` and the number of species of kind `i` is
`length(posidx[i]) - length(freespecies[i])`.
"""
struct SimulationStep{N,T}
    positions::Vector{SVector{3,TÅ}}
end


function SimulationStep(inputpos::Vector{Vector{Vector{SVector{N,TÅ}}}}) where N

    numatoms = sum(x -> sum(length, x; init=0), inputpos; init=0)
    positions = Vector{SVector{N,TÅ}}(undef, numatoms)
    l = 0
    for (i, posi) in enumerate(inputpos)
        for (j, poss) in enumerate(posi)
            for (k, pos) in enumerate(poss)
                l += 1
                positions[l] = pos
            end
        end
    end

    SimulationStep{N,TÅ}(positions)
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
    SimulationStep(poss), indices
end

"""
    SimulationStep(step::SimulationStep, mode=:all)

Copy of the input on which the positions and number of species can be modified without
affecting the original `step`.

If `mode === :output`, only the `positions` and `atoms` fields are copied.
If `mode === :zero`, create a `SimulationStep` with an empty `ffidx` field.

The `parallel` field is passed on to the created copy (except with `mode === :zero`)
"""
function SimulationStep(step::SimulationStep{N,T}, mode=:all) where {N,T}
    mode === :zero ? step : deepcopy(step)
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
        l = molpos[k]
        step.positions[l] = pos
    end
end

"""
    update_position(step::SimulationStep, idx, newpos)

Return a new `SimulationStep` where the system of index `idx` has a new position `newpos`.

`idx` can be the index returned by [`make_step`](@ref) or a tuple `(i,j)` which designates
the `j`-th system of kind `i` in `step`.

See also [`update_position!`](@ref) to modify `step` in-place.
"""
function update_position(step::SimulationStep{N,T}, (i,j), newpos) where {N,T}
    x = SimulationStep{N,T}(copy(step.positions))
    update_position!(x, (i,j), newpos)
    x
end

function compute_vdw(step::SimulationStep)
    energy = 0.0u"K"
    cutoff2 = step.ff.cutoff^2
    buffer = MVector{3,TÅ}(undef)
    buffer2 = MVector{3,Float64}(undef)
    cell = CellMatrix(step.mat)
    for (l1, pos1) in enumerate(step.positions)
        i1, j1, k1 = step.atoms[l1]
        ix1 = step.ffidx[i1][k1]
        for l2 in (l1+1):length(step.positions)
            i2, j2, k2 = step.atoms[l2]
            i1 == i2 && j1 == j2 && continue # no intra energy considered here
            pos2 = step.positions[l2]
            buffer .= pos2 .- pos1
            d2 = unsafe_periodic_distance2!(buffer, buffer2, cell)
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


function single_contribution_vdw(step::SimulationStep, idx2::Tuple{Int,Int}, poss2::AbstractVector{SVector{N,TÅ}} where N)
    return rand()*1e5*u"K"
end
