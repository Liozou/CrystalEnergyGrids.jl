export MonteCarloSetup, setup_montecarlo, baseline_energy, movement_energy, run_montecarlo!

struct MonteCarloSetup{N,T,Trng}
    step::SimulationStep{N,T}
    # step contains all the information that is not related to the framework nor to Ewald.
    # It contains all the information necessary to compute the species-species VdW energy.
    indices::Int
    # indices is the set of all (i,j) with 1 ≤ j ≤ number of species of kind i.
    rng::Trng
end

function Base.show(io::IO, ::MonteCarloSetup)
    print(io, "Monte-Carlo setup")
end


function setup_montecarlo(systems)
    cell = SMatrix{3,3,Float64,9}(25, 0, 0, 0, 25, 0, 0, 0, 25)
    rng = default_rng()

    U = Vector{SVector{3,Float64}} # positions of the atoms of a system
    poss = Vector{U}[U[[rand(SVector{3,Float64})]] for _ in 1:systems]

    n = length(poss)

    MonteCarloSetup(SimulationStep(poss, cell), n, rng)
end


function random_translation(rng, positions::AbstractVector{SVector{3,Float64}}, dmax::Float64)
    r = SVector{3}(((2*rand(rng)-1)*dmax) for _ in 1:3)
    [poss + r for poss in positions]
end
