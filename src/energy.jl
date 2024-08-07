# Computation of energy resulting from a system with a given force field

export SimulationStep, make_step, update_position, update_position!
using CellListMap
import LinearAlgebra

"""
    SimulationStep{T}

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
- `positions`: a `Vector` of `SVector{3,typeof(1.0u"Å")}` representing the coordinates of
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
struct SimulationStep{T}
    ff::ForceField
    charges::Vector{Te_au}
    # charges[ix] is the charge of the atom of index ix in ff, i.e. the charge of the k-th
    # atom in a system of kind i is charges[ffidx[i][k]].
    psystem::T
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
function Base.propertynames(::SimulationStep{T}, private::Bool=false) where T
    private && return fieldnames(SimulationStep{T})
    (:ff, :charges, :mat, :positions, :parallel, :atoms, :posidx, :freespecies, :ffidx)
end

function Base.:(==)(s1::SimulationStep, s2::SimulationStep)
    s1.positions == s2.positions &&
    s1.ff == s2.ff &&
    mapreduce(===, &, s1.charges, s2.charges) # consider NaN as equals
    s1.mat == s2.mat &&
    s1.atoms == s2.atoms &&
    s1.posidx == s2.posidx &&
    s1.freespecies == s2.freespecies &&
    s1.ffidx == s2.ffidx
    # do not check for equality of "parallel" field
end


function SimulationStep(ff::ForceField, charges::Vector{Te_au},
                        inputpos::Vector{Vector{Vector{SVector{3,TÅ}}}},
                        isrigid::BitVector, ffidx::Vector{Vector{Int}}, cell::CellMatrix;
                        parallel::Bool=true)

    numatoms = sum(x -> sum(length, x; init=0), inputpos; init=0)
    positions = Vector{SVector{3,TÅ}}(undef, numatoms)
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
    psystem = ParticleSystem(; xpositions=positions,
                               ypositions=SVector{3,TÅ}[],
                               unitcell=cell.mat,
                               cutoff=ff.cutoff,
                               parallel,
                               output=0.0u"K")

    SimulationStep{typeof(psystem)}(ff, charges, psystem, atoms, posidx, freespecies, isrigid, ffidx)
end

"""
    SimulationStep(ff::ForceField, systemkinds::Vector{T} where T<:AbstractSystem,
                   inputpos::Vector{Vector{Vector{SVector{3,TÅ}}}},
                   cell::CellMatrix,
                   isrigid::BitVector=trues(length(systemkinds));
                   parallel::Bool=true)

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
                        inputpos::Vector{Vector{Vector{SVector{3,TÅ}}}},
                        cell::CellMatrix,
                        isrigid::BitVector=trues(length(systemkinds));
                        parallel::Bool=true)
    @assert length(systemkinds) == length(inputpos) == length(isrigid)
    ffidx = [[ff.sdict[atomic_symbol(s, k)] for k in 1:length(s)] for s in systemkinds]
    charges = fill(NaN*u"e_au", length(ff.symbols))
    for (syst, ffidxi) in zip(systemkinds, ffidx)
        for (k, ix) in enumerate(ffidxi)
            isnan(charges[ix]) || continue
            charges[ix] = uconvert(u"e_au", syst[k,:atomic_charge])::Te_au
        end
    end
    SimulationStep(ff, charges, inputpos, isrigid, ffidx, cell; parallel)
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
function make_step(ff::ForceField, systems::Vector{T}, cell::CellMatrix=CellMatrix(),
                   isrigid::BitVector=trues(length(systems)); parallel::Bool=true) where T<:AbstractSystem
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
    SimulationStep(ff, systemkinds, poss, cell, newisrigid; parallel), indices
end

"""
    SimulationStep(step::SimulationStep, mode=:all; parallel=step.parallel)

Copy of the input on which the positions and number of species can be modified without
affecting the original `step`.

If `mode === :output`, only the `positions` and `atoms` fields are copied.
If `mode === :zero`, create a `SimulationStep` with an empty `ffidx` field.

The `parallel` field is passed on to the created copy (except with `mode === :zero`)
"""
function SimulationStep(step::SimulationStep{T}, mode=:all; parallel=step.parallel) where T
    if mode === :all
        psystem = ParticleSystem(; xpositions=copy(step.positions),
                                   ypositions=SVector{3,TÅ}[],
                                   unitcell=step.mat,
                                   parallel,
                                   cutoff=step.ff.cutoff, output=0.0u"K")
        SimulationStep{T}(step.ff, step.charges, psystem, copy(step.atoms),
                       [[copy(js) for js in is] for is in step.posidx],
                       [copy(x) for x in step.freespecies], step.isrigid, step.ffidx)
    elseif mode === :output
        psystem = ParticleSystem(; xpositions=copy(step.positions),
                                   ypositions=SVector{3,TÅ}[],
                                   unitcell=copy(step.mat),
                                   parallel,
                                   cutoff=step.ff.cutoff, output=0.0u"K")
        SimulationStep{T}(step.ff, step.charges, psystem, copy(step.atoms),
                       step.posidx, step.freespecies, step.isrigid, step.ffidx)
    elseif mode === :complete_output
        SimulationStep{T}(step.ff, step.charges, step.psystem, step.atoms,
                            [[copy(js) for js in is] for is in step.posidx],
                            [copy(x) for x in step.freespecies], step.isrigid, step.ffidx)
    elseif mode === :zero
        SimulationStep{T}(step.ff, step.charges, step.psystem, step.atoms,
                       step.posidx, step.freespecies, step.isrigid, Vector{Int}[])
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
function update_position(step::SimulationStep{T}, (i,j), newpos) where T
    psystem = ParticleSystem(; xpositions=copy(step.positions),
                               ypositions=SVector{3,TÅ}[],
                               unitcell=step.mat,
                               parallel=step.parallel,
                               cutoff=step.ff.cutoff, output=0.0u"K")

    x = SimulationStep{T}(step.ff, step.charges, psystem, step.atoms, step.posidx, step.freespecies, step.isrigid, step.ffidx)
    update_position!(x, (i,j), newpos)
    x
end

function energy_intra(step::SimulationStep, i::Int, positions::Vector{SVector{3,TÅ}})
    n = length(positions)
    energy = 0.0u"K"
    cutoff2 = step.ff.cutoff^2
    isinf(step.ff.cutoff) || error("finite cutoff not implemented")
    indices = step.ffidx[i]
    for k1 in 1:n
        ix1 = indices[k1]
        for k2 in (i+1):n
            ix2 = indices[k2]
            r2 = norm2(positions[k1], positions[k2])
            r2 > cutoff2 && error("Flexible molecule cannot be larger than the cutoff")
            rule = step.ff[ix1, ix2]
            energy += rule(r2)# + tailcorrection(rule, step.ff.cutoff)
        end
    end
    energy
end

struct TotalVdwComputation{T} <: Function
    step::SimulationStep{T}
end
function (vdw::TotalVdwComputation)(pos1, pos2, l1, l2, d², output)
    i1, j1, k1 = vdw.step.atoms[l1]
    i2, j2, k2 = vdw.step.atoms[l2]
    i1 == i2 && j1 == j2 && return output # internal energy
    output += vdw.step.ff[vdw.step.ffidx[i1][k1], vdw.step.ffidx[i2][k2]](d²)
    output
end

function compute_vdw(step::SimulationStep)
    (step.parallel*nthreads() + !step.parallel)*length(step.atoms) < 1200 && return compute_vdw_noneighbour(step)
    system = ParticleSystem(;
        xpositions=step.positions,
        unitcell=step.mat,
        cutoff=step.psystem.cutoff,
        parallel=step.parallel,
        output=0.0u"K"
    )
    map_pairwise(TotalVdwComputation(step), system)
    system.output
end


function compute_vdw_noneighbour(step::SimulationStep)
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


struct SingleVdwComputation{T} <: Function
    step::SimulationStep{T}
    i2::Int
    j2::Int
    idx2::Vector{Int}
end
function (vdw::SingleVdwComputation)(pos1, pos2, l1, k2, d², output)
    i1, j1, k1 = vdw.step.atoms[l1]
    i1 == vdw.i2 && j1 == vdw.j2 && return output
    output += vdw.step.ff[vdw.step.ffidx[i1][k1], vdw.idx2[k2]](d²)
    output
end

function single_contribution_vdw(step::SimulationStep, idx2::Tuple{Int,Int}, poss2::AbstractVector{SVector{3,TÅ}})
    (step.parallel*nthreads() + !step.parallel)*length(step.atoms) < 1200 && return single_contribution_vdw_noneighbour(step, idx2, poss2)
    system = step.psystem
    resize!(system.ypositions, length(poss2))
    system.ypositions .= poss2
    i2, j2 = idx2
    map_pairwise(SingleVdwComputation(step, i2, j2, step.ffidx[i2]), system)
    system.output
end

function single_contribution_vdw_noneighbour(step::SimulationStep, (i2,j2)::Tuple{Int,Int}, poss2::AbstractVector{SVector{3,TÅ}})
    energy = step.isrigid[i2] ? 0.0u"K" : energy_intra(step, i2, poss2)
    cutoff2 = step.ff.cutoff^2
    buffer = MVector{3,TÅ}(undef)
    buffer2 = MVector{3,Float64}(undef)
    idx2 = step.ffidx[i2]
    cell = CellMatrix(step.mat)
    for (k2, pos2) in enumerate(poss2)
        ix2 = idx2[k2]
        for (l1, pos1) in enumerate(step.positions)
            i1, j1, k1 = step.atoms[l1]
            i1 == i2 && j1 == j2 && continue
            buffer .= pos2 .- pos1
            d2 = unsafe_periodic_distance2!(buffer, buffer2, cell)
            if d2 < cutoff2
                energy += step.ff[step.ffidx[i1][k1], ix2](d2)
            end
        end
    end
    energy
end
