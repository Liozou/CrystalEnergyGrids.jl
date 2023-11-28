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
    ff::ForceField
    charges::Vector{Te_au}
    # charges[ix] is the charge of the atom of index ix in ff, i.e. the charge of the k-th
    # atom in a system of kind i is charges[ffidx[i][k]].
    psystem::T
    atoms::Vector{Tuple{Int,Int,Int}} # one index per atom
    posidx::Vector{Vector{Vector{Int}}}
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


function SimulationStep(ff::ForceField, charges::Vector{Te_au},
                        inputpos::Vector{Vector{Vector{SVector{N,TÅ}}}},
                        cell::CellMatrix;
                        parallel::Bool=true) where N

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
    psystem = PeriodicSystem(; xpositions=positions,
                               ypositions=SVector{3,TÅ}[],
                               unitcell=cell.mat,
                               cutoff=ff.cutoff,
                               parallel,
                               output=0.0u"K")

    SimulationStep{N,typeof(psystem)}(ff, charges, psystem, atoms, posidx)
end

"""
    SimulationStep(ff::ForceField, systemkinds::Vector{T} where T<:AbstractSystem,
                   inputpos::Vector{Vector{Vector{SVector{N,TÅ}}}},
                   cell::CellMatrix,
                   isrigid::BitVector=trues(length(systemkinds));
                   parallel::Bool=true) where N

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

`parallel` should be unset to do all computations on the step on a single thread.

See also [`make_step`](@ref) for an alternative way of constructing a `SimulationStep`
without having to specify the system kinds.
"""
function SimulationStep(ff::ForceField, systemkinds::Vector{T} where T<:AbstractSystem,
                        inputpos::Vector{Vector{Vector{SVector{N,TÅ}}}},
                        cell::CellMatrix;
                        parallel::Bool=true) where N
    @assert length(systemkinds) == length(inputpos)
    charges = [[uconvert(u"e_au", s[k,:atomic_charge])::Te_au for k in 1:length(s)] for s in systemkinds]
    SimulationStep(ff, charges, inputpos, cell; parallel)
end

"""
    make_step(ff::ForceField, systems::Vector{T}, cell::CellMatrix=CellMatrix(), isrigid::BitVector=trues(length(systems)); parallel::Bool=true) where T<:AbstractSystem

Create a [`SimulationStep`](@ref) from a given force field `ff` and a list of interacting
`systems`. `isrigid[i]` specifies whether `systems[i]` is considered rigid. If unspecified,
all systems are considered rigid (so only inter-system interactions are taken into account).

Return a tuple `(step, indices)` where `step` is the `SimulationStep` and `indices[i]` is
the index of `systems[i]` in `step`, for use with [`update_position`](@ref) and such.

Systems sharing the same ordered list of atom symbols and the same rigidity will be
considered as having the same system kind (e.g. several molecules of water).

If unspecified, `cell` is set to an open (i.e. aperiodic) box.
"""
function make_step(ff::ForceField, systems::Vector{T}, cell::CellMatrix=CellMatrix(); parallel::Bool=true) where T<:AbstractSystem
    kindsdict = Dict{Tuple{Vector{Symbol},Bool},Int}()
    systemkinds = T[]
    U = Vector{SVector{n_dimensions(systems[1]),TÅ}} # positions of the atoms of a system
    poss = Vector{U}[]
    indices = Tuple{Int,Int}[]
    for system in systems
        n = length(kindsdict)+1
        kind = get!(kindsdict, atomic_symbol(system), n)
        if kind === n
            push!(systemkinds, system)
            push!(poss, U[])
        end
        push!(poss[kind], position(system))
        push!(indices, (kind, length(poss[kind])))
    end
    SimulationStep(ff, systemkinds, poss, cell; parallel), indices
end

"""
    SimulationStep(step::SimulationStep, mode=:all; parallel=step.parallel)

Copy of the input on which the positions and number of species can be modified without
affecting the original `step`.

If `mode === :output`, only the `positions` and `atoms` fields are copied.
If `mode === :zero`, create a `SimulationStep` with an empty `ffidx` field.

The `parallel` field is passed on to the created copy (except with `mode === :zero`)
"""
function SimulationStep(step::SimulationStep{N,T}, mode=:all; parallel=step.parallel) where {N,T}
    if mode === :all
        return deepcopy(step)
        psystem = PeriodicSystem(; xpositions=copy(step.positions),
                                   ypositions=SVector{3,TÅ}[],
                                   unitcell=step.mat,
                                   parallel,
                                   cutoff=step.ff.cutoff, output=0.0u"K")
        SimulationStep{N,T}(step.ff, step.charges, psystem, copy(step.atoms),
                       [[copy(js) for js in is] for is in step.posidx])
    elseif mode === :output
        return deepcopy(step)
        psystem = PeriodicSystem(; xpositions=copy(step.positions),
                                   ypositions=SVector{3,TÅ}[],
                                   unitcell=step.mat,
                                   parallel,
                                   cutoff=step.ff.cutoff, output=0.0u"K")
        SimulationStep{N,T}(step.ff, step.charges, psystem, copy(step.atoms),
                       step.posidx)
    elseif mode === :complete_output
        return deepcopy(step)
        SimulationStep{N,T}(step.ff, step.charges, step.psystem, step.atoms,
                            [[copy(js) for js in is] for is in step.posidx])
    elseif mode === :zero
        SimulationStep{N,T}(step.ff, step.charges, step.psystem, step.atoms,
                       step.posidx)
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
    psystem = PeriodicSystem(; xpositions=copy(step.positions),
                               ypositions=SVector{3,TÅ}[],
                               unitcell=step.mat,
                               parallel=step.parallel,
                               cutoff=step.ff.cutoff, output=0.0u"K")

    x = SimulationStep{N,T}(step.ff, step.charges, psystem, step.atoms, step.posidx)
    update_position!(x, (i,j), newpos)
    x
end


function single_contribution_vdw(step::SimulationStep, idx2::Tuple{Int,Int}, poss2::AbstractVector{SVector{N,TÅ}} where N)
    return rand()*1e5*u"K"
end
