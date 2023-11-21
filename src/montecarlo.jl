struct MonteCarloSetup{Trng}
    positions::Vector{SVector{3,Float64}}
    rng::Trng
end

function setup_montecarlo()
    MonteCarloSetup([rand(SVector{3,Float64}) for _ in 1:80], default_rng())
end
