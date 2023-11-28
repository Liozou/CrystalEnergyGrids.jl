# Computation of energy resulting from a system with a given force field

export SimulationStep, make_step, update_position, update_position!
using CellListMap.PeriodicSystems
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
    # charges[ix] is the charge of the atom of index ix in ff, i.e. the charge of the k-th
    # atom in a system of kind i is charges[ffidx[i][k]].
    psystem::T
    atoms::Vector{Tuple{Int,Int,Int}} # one index per atom
    # the position of the k-th atom in the j-th system of kind i is positions[l] where
    # atoms[l] == (i,j,k) and posidx[i][j][k] == l
end

@inline function Base.getproperty(s::SimulationStep, x::Symbol)
    if x === :positions
        return getfield(s, :psystem).xpositions
    elseif x === :mat
        return getfield(s, :psystem).unitcell
    elseif x === :parallel
        return getfield(s, :psystem).parallel
    else
        return getfield(s, x)
    end
end

function Base.show(io::IO, step::SimulationStep)
    n = length(step.atoms)
    print(io, "Simulation step with ", n , " atoms in ", m, " molecule kind")
    m > 1 && print(io, 's')
end


function SimulationStep(inputpos::Vector{Vector{Vector{SVector{N,TÅ}}}},
                        cell::CellMatrix;
                        parallel::Bool=true) where N

    numatoms = sum(x -> sum(length, x; init=0), inputpos; init=0)
    positions = Vector{SVector{N,TÅ}}(undef, numatoms)
    atoms = Vector{Tuple{Int,Int,Int}}(undef, numatoms)
    l = 0
    for (i, posi) in enumerate(inputpos)
        for (j, poss) in enumerate(posi)
            for (k, pos) in enumerate(poss)
                l += 1
                atoms[l] = (i,j,k)
                positions[l] = pos
            end
        end
    end
    psystem = PeriodicSystem(; xpositions=positions,
                               ypositions=SVector{3,TÅ}[],
                               unitcell=cell.mat,
                               cutoff=12.0u"Å",
                               parallel,
                               output=0.0u"K")

    SimulationStep{N,typeof(psystem)}(psystem, atoms)
end


"""
    update_position!(step::SimulationStep, idx, newpos)

Modifies the position of the system of index `idx` in `step` to `newpos`.
Return `newpos`.

See also [`update_position!(step::SimulationStep, idx, op, arg)`](@ref) and [`update_position`](@ref).
"""
function update_position!(step::SimulationStep, (i,j), newpos)
    for pos in newpos
        step.positions[i] = pos
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
    psystem = PeriodicSystem(; xpositions=copy(step.positions),
                               ypositions=SVector{3,TÅ}[],
                               unitcell=step.mat,
                               parallel=step.parallel,
                               cutoff=12.0u"Å", output=0.0u"K")

    x = SimulationStep{N,T}(psystem, step.atoms)
    update_position!(x, (i,j), newpos)
    x
end
