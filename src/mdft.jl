import Optim
using FFTW
import LinearAlgebra: norm

abstract type MDFTProblem <: Function end

struct IdealGasProblem <: MDFTProblem
    T::Float64    # in K
    P::Float64    # in bar
    ρ₀::Float64   # in mol/L
    externalV::Array{Float64,3} # energies in K
end
function IdealGasProblem(T::Float64, P::Float64, externalV::Array{Float64,3})
    IdealGasProblem(T, P, P/(0.0831446261815324*T), externalV)
end
function IdealGasProblem(gasname::String, T::Float64, P::Float64, externalV::Array{Float64,3})
    a, b = get(VDW_COEFF, gasname, (0.0, 0.0))
    a == 0.0 && error(lazy"""Gas $gasname does not appear in the database, please enter an explicit ideal density or omit the name to rely on the perfect gas model. Possible gases are: $(join(keys(VDW_COEFF, ", "))). Please open an issue if you want to add the Van der Waals coefficient corresponding to your gas.""")
    c = 0.0831446261815324*T # perfect gas constants in bar*L/K/mol
    # The following is the real root of (Px^2 + a)(x-b) = c*x^2 given by Wolfram Alpha
    x = cbrt(sqrt((18*a*b*P^2 - 9*a*c*P + 2*b^3*P^3 + 6*b^2*c*P^2 + 6*b*c^2*P + 2*c^3)^2 + 4*(3*a*P - (-b*P - c)^2)^3) + 18*a*b*P^2 - 9*a*c*P + 2*b^3*P^3 + 6*b^2*c*P^2 + 6*b*c^2*P + 2*c^3)/(3*cbrt(2)*P) - (cbrt(2)*(3*a*P - (-b*P - c)^2))/(3*P*cbrt(sqrt((18*a*b*P^2 - 9*a*c*P + 2*b^3*P^3 + 6*b^2*c*P^2 + 6*b*c^2*P + 2*c^3)^2 + 4*(3*a*P - (-b*P - c)^2)^3) + 18*a*b*P^2 - 9*a*c*P + 2*b^3*P^3 + 6*b^2*c*P^2 + 6*b*c^2*P + 2*c^3)) - (-b*P - c)/(3*P)
    IdealGasProblem(T, P, inv(x), externalV)
end


function (igp::IdealGasProblem)(F, flat_∂ψ, flat_ψ::Vector{Float64})
    try
    ψ = reshape(flat_ψ, size(igp.externalV))
    gradient = isnothing(flat_∂ψ) ? similar(ψ) : reshape(flat_∂ψ, size(igp.externalV))
    gradient .= exp.(ψ)
    ψ₀ = log(igp.ρ₀)
    value = sum(expr*v + igp.T*(expr*(r - ψ₀ - 1) + igp.ρ₀) for (r, v, expr) in zip(ψ, igp.externalV, gradient))
    if !isnothing(flat_∂ψ)
        @. gradient *= igp.externalV + igp.T*(ψ - ψ₀)
    end
    # isempty(flat_∂ψ) || @show extrema(flat_∂ψ)
    # @show value, extrema(ψ), norm(vec(gradient))
    value
    catch e
        Base.showerror(stderr, e)
        Base.show_backtrace(stderr, catch_backtrace())
        rethrow()
    end
end

function mdft(igp::IdealGasProblem)
    ψ₀ = log(igp.ρ₀)
    η = 1e-9
    ψ_init = (ψ₀ .- vec(igp.externalV)./igp.T)# .* rand(1 .+ (-η:η), length(igp.externalV))
    Optim.optimize(Optim.only_fg!(igp), ψ_init, Optim.LBFGS(), Optim.Options(iterations=100))
end

