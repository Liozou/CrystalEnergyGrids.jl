export MonteCarloSetup, setup_montecarlo, baseline_energy, movement_energy, run_montecarlo!

struct MonteCarloSetup{Trng}
    positions::Vector{SVector{3,Float64}}
    rng::Trng
end

function setup_montecarlo()
    MonteCarloSetup([rand(SVector{3,Float64}) for _ in 1:80], default_rng())
end

function MonteCarloSetup(mc::MonteCarloSetup, o=mc.positions)
    rng = deepcopy(mc.rng)
    MonteCarloSetup(copy(o), rng)
end


choose_random_species(mc::MonteCarloSetup) = rand(mc.rng)
