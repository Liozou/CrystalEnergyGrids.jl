# Computation of energy resulting from a system with a given force field

export SimulationStep, make_step, update_position, update_position!
using CellListMap.PeriodicSystems
import LinearAlgebra

struct SimulationStep{N,T}
    ff::ForceField
    charges::Vector{Te_au}
    # charges[ix] is the charge of the atom of index ix in ff, i.e. the charge of the k-th
    # atom in a system of kind i is charges[ffidx[i][k]].
    psystem::T
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

function Base.show(io::IO, ::SimulationStep)
    print(io, "Simulation step")
end


function SimulationStep(ff::ForceField, charges::Vector{Te_au},
                        inputpos::Vector{Vector{Vector{SVector{N,TÅ}}}},
                        cell::CellMatrix;
                        parallel::Bool=true) where N

    numatoms = sum(x -> sum(length, x; init=0), inputpos; init=0)
    positions = Vector{SVector{N,TÅ}}(undef, numatoms)
    posidx = Vector{Vector{Vector{Int}}}(undef, length(inputpos))
    l = 0
    for (i, posi) in enumerate(inputpos)
        posidxi = Vector{Vector{Int}}(undef, length(posi))
        posidx[i] = posidxi
        for (j, poss) in enumerate(posi)
            posidxi[j] = collect(l+1:l+length(poss))
            for (k, pos) in enumerate(poss)
                l += 1
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

    SimulationStep{N,typeof(psystem)}(ff, charges, psystem, posidx)
end

function SimulationStep(ff::ForceField, systemkinds::Vector{T} where T<:AbstractSystem,
                        inputpos::Vector{Vector{Vector{SVector{N,TÅ}}}},
                        cell::CellMatrix;
                        parallel::Bool=true) where N
    @assert length(systemkinds) == length(inputpos)
    charges = [[uconvert(u"e_au", s[k,:atomic_charge])::Te_au for k in 1:length(s)] for s in systemkinds]
    SimulationStep(ff, charges, inputpos, cell; parallel)
end

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
        SimulationStep{N,T}(step.ff, step.charges, psystem,
                       [[copy(js) for js in is] for is in step.posidx])
    elseif mode === :output
        return deepcopy(step)
        psystem = PeriodicSystem(; xpositions=copy(step.positions),
                                   ypositions=SVector{3,TÅ}[],
                                   unitcell=step.mat,
                                   parallel,
                                   cutoff=step.ff.cutoff, output=0.0u"K")
        SimulationStep{N,T}(step.ff, step.charges, psystem,
                       step.posidx)
    elseif mode === :complete_output
        return deepcopy(step)
        SimulationStep{N,T}(step.ff, step.charges, step.psystem,
                            [[copy(js) for js in is] for is in step.posidx])
    elseif mode === :zero
        SimulationStep{N,T}(step.ff, step.charges, step.psystem,
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

    x = SimulationStep{N,T}(step.ff, step.charges, psystem, step.posidx)
    update_position!(x, (i,j), newpos)
    x
end

