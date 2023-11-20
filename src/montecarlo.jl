export MonteCarloSetup, setup_montecarlo, baseline_energy, movement_energy, run_montecarlo!

struct MonteCarloSetup{N,T,Trng}
    step::SimulationStep{N,T}
    # step contains all the information that is not related to the framework nor to Ewald.
    # It contains all the information necessary to compute the species-species VdW energy.
    offsets::Vector{Int}
    # offsets[i] is the number of molecules belonging to a kind strictly lower than i.
    indices::Set{Tuple{Int,Int}}
    # indices is the set of all (i,j) with 1 ≤ j ≤ number of species of kind i.
    mcmoves::Vector{MCMoves} # mcmoves[i] is the set of MC moves for kind i
    rng::Trng
end

function Base.show(io::IO, mc::MonteCarloSetup)
    n = length(mc.indices)
    m = length(mc.offsets)
    print(io, "Monte-Carlo setup with ", n , " atoms in ", m, " molecule kind")
    m > 1 && print(io, 's')
end


struct EwaldSystem # pseudo-AbstractSystem with only positions and charges
    position::Vector{SVector{3,TÅ}}
    atomic_charge::Vector{Te_au}
end
AtomsBase.position(x::EwaldSystem) = x.position
AtomsBase.position(x::EwaldSystem, i::Int) = x.position[i]
Base.getindex(x::EwaldSystem, ::Colon, s::Symbol) = getproperty(x, s)
Base.getindex(x::EwaldSystem, i::Int, s::Symbol) = getproperty(x, s)[i]
Base.length(x::EwaldSystem) = length(x.position)

struct IdSystem # pseudo-AbstractSystem with only atomic symbols and charges
    atomic_symbol::Vector{Symbol}
    atomic_charge::Vector{Te_au}
end
IdSystem(s::AbstractSystem{3}) = IdSystem(atomic_symbol(s), s[:,:atomic_charge])
Base.length(s::IdSystem) = length(s.atomic_charge)

function setup_montecarlo(systems::Vector; rng=default_rng())
    kindsdict = Dict{Tuple{Vector{Symbol},MCMoves},Int}()
    systemkinds = IdSystem[]
    U = Vector{SVector{3,TÅ}} # positions of the atoms of a system
    poss = Vector{U}[]
    indices = Tuple{Int,Int}[]
    rev_indices = Vector{Int}[]
    newmcmoves = MCMoves[]
    for (i, s) in enumerate(systems)
        system, n = s isa Tuple ? s : (s, 1)
        m = length(kindsdict)+1
        mcmove = MCMoves(length(systems) == 1)
        kind = get!(kindsdict, (atomic_symbol(system)::Vector{Symbol}, mcmove), m)
        if kind === m
            push!(systemkinds, IdSystem(system)::IdSystem)
            push!(poss, U[])
            push!(rev_indices, Int[])
            push!(newmcmoves, mcmove)
        end
        push!(indices, (kind, length(poss[kind])+1))
        append!(rev_indices[kind], i for _ in 1:n)
        append!(poss[kind], copy(position(system)::Vector{SVector{3,TÅ}}) for _ in 1:n)
    end

    n = length(poss)
    offsets = Vector{Int}(undef, n)
    offsets[1] = 0
    indices_list = Tuple{Int,Int}[]
    for i in 1:(n-1)
        numi = length(poss[i])
        offsets[i+1] = offsets[i] + numi
        append!(indices_list, (i,j) for j in 1:length(poss[i]))
    end
    append!(indices_list, (n,j) for j in 1:length(poss[n]))

    MonteCarloSetup(SimulationStep(poss), offsets, Set(indices_list), newmcmoves, rng), indices
end


"""
    MonteCarloSetup(mc::MonteCarloSetup)

Create a copy of `mc` that does not share its modifiable internal states (the positions and
the Ewald state). For example, the copy and the original can be used to run Monte-Carlo
simulations in parallel, and the state of one will not influence that of the other.

`parallel` specifies whether the computations on the resulting `MonteCarloSetup` should be
parallelized or not.

!!! warn
    Internal states that are semantically immutable are shared, although some of them are
    technically mutable, like the value of the atom charges for instance. Modifying such a
    state on the original or the copy can thus also modify that of the other: use
    `deppcopy` or a custom copy implementation to circumvent this issue if you plan on
    modifying states outside of the API.
"""
function MonteCarloSetup(mc::MonteCarloSetup, o::SimulationStep=mc.step)
    rng = deepcopy(mc.rng)
    MonteCarloSetup(SimulationStep(o, :all), copy(mc.offsets),
                    copy(mc.indices), copy(mc.mcmoves), rng)
end


choose_random_species(mc::MonteCarloSetup) = rand(mc.rng, mc.indices)
