export MonteCarloSetup, setup_montecarlo, baseline_energy, movement_energy, run_montecarlo!

struct MonteCarloSetup{N,T,Trng}
    step::SimulationStep{N,T}
    # step contains all the information that is not related to the framework nor to Ewald.
    # It contains all the information necessary to compute the species-species VdW energy.
    offsets::Vector{Int}
    # offsets[i] is the number of molecules belonging to a kind strictly lower than i.
    indices::Set{Tuple{Int,Int}}
    # indices is the set of all (i,j) with 1 ≤ j ≤ number of species of kind i.
    bead::Vector{Int} # k = bead[i] is the number of the reference bead of kind i.
    rng::Trng
end

function Base.show(io::IO, mc::MonteCarloSetup)
    n = length(mc.indices)
    m = length(mc.offsets)
    print(io, "Monte-Carlo setup with ", n , " atoms in ", m, " molecule kind")
    m > 1 && print(io, 's')
end


function setup_montecarlo(systems)
    cell = CellMatrix(SMatrix{3,3,Float64,9}(25, 0, 0, 0, 25, 0, 0, 0, 25).*u"Å")
    parallel = true
    rng = default_rng()

    U = Vector{SVector{3,TÅ}} # positions of the atoms of a system
    poss = Vector{U}[U[[rand(SVector{3,TÅ})]] for _ in 1:systems]
    indices = Tuple{Int,Int}[(1, 1)]

    ffidx = [[1] for _ in 1:systems]
    charges = [NaN*u"e_au"]

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

    beads = fill(1, n)


    MonteCarloSetup(SimulationStep(ForceField(), charges, poss, trues(length(ffidx)), ffidx, cell; parallel),
                    offsets, Set(indices_list), beads, rng), indices
end


"""
    MonteCarloSetup(mc::MonteCarloSetup; parallel::Bool=mc.step.parallel)

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
function MonteCarloSetup(mc::MonteCarloSetup, o::SimulationStep=mc.step; parallel::Bool=mc.step.parallel)
    rng = deepcopy(mc.rng)
    MonteCarloSetup(SimulationStep(o, :all; parallel), copy(mc.offsets),
                    copy(mc.indices), copy(mc.bead), rng)
end

function set_position!(mc::MonteCarloSetup, (i, j), newpositions, newEiks=nothing)
    molpos = mc.step.posidx[i][j]
    for (k, newpos) in enumerate(newpositions)
        mc.step.positions[molpos[k]] = if eltype(newpositions) <: AbstractVector{<:AbstractFloat}
            newpos
        else
            NoUnits(newpos/u"Å")
        end
    end
end


choose_random_species(mc::MonteCarloSetup) = rand(mc.rng, mc.indices)

function random_translation(rng, positions::AbstractVector{SVector{3,TÅ}}, dmax::TÅ)
    r = SVector{3}(((2*rand(rng)-1)*dmax) for _ in 1:3)
    [poss + r for poss in positions]
end
