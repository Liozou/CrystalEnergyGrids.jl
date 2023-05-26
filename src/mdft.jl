import Optim
using FFTW, Hankel
import BasicInterpolators: CubicSplineInterpolator, NoBoundaries, WeakBoundaries
using StaticArrays
import LinearAlgebra: norm, mul!
using Base.Threads

FFTW.set_num_threads(nthreads()÷2)

export IdealGasProblem, MonoAtomic

abstract type MDFTProblem <: Function end

struct IdealGasProblem <: MDFTProblem
    T::Float64    # temperature, in K
    P::Float64    # pressure, in bar
    ρ₀::Float64   # reference number density, in number/Å³ (obtained from VdW coefficients)
    externalV::Array{Float64,3} # external energy field in K
    δv::Float64   # volume of an elementary grid cube, in Å³
end
function IdealGasProblem(T::Float64, P::Float64, externalV::Array{Float64,3}, δv::Float64)
    IdealGasProblem(T, P, P/(0.0831446261815324*T), externalV, δv)
end
function IdealGasProblem(gasname::String, T::Float64, P::Float64, externalV::Array{Float64,3}, δv::Float64)
    a, b = get(VDW_COEFF, gasname, (0.0, 0.0))
    a == 0.0 && error(lazy"""Gas $gasname does not appear in the database, please enter an explicit ideal density or omit the name to rely on the perfect gas model. Possible gases are: $(join(keys(VDW_COEFF, ", "))). Please open an issue if you want to add the Van der Waals coefficient corresponding to your gas.""")
    c = 0.0831446261815324*T # perfect gas constants in bar*L/K/mol
    # The following is the real root of (Px^2 + a)(x-b) = c*x^2 given by Wolfram Alpha
    x = cbrt(sqrt((18*a*b*P^2 - 9*a*c*P + 2*b^3*P^3 + 6*b^2*c*P^2 + 6*b*c^2*P + 2*c^3)^2 + 4*(3*a*P - (-b*P - c)^2)^3) + 18*a*b*P^2 - 9*a*c*P + 2*b^3*P^3 + 6*b^2*c*P^2 + 6*b*c^2*P + 2*c^3)/(3*cbrt(2)*P) - (cbrt(2)*(3*a*P - (-b*P - c)^2))/(3*P*cbrt(sqrt((18*a*b*P^2 - 9*a*c*P + 2*b^3*P^3 + 6*b^2*c*P^2 + 6*b*c^2*P + 2*c^3)^2 + 4*(3*a*P - (-b*P - c)^2)^3) + 18*a*b*P^2 - 9*a*c*P + 2*b^3*P^3 + 6*b^2*c*P^2 + 6*b*c^2*P + 2*c^3)) - (-b*P - c)/(3*P)
    # inv(x) is ρ₀ in mol/L
    ρ₀ = inv(x) * 6.02214076e-4 # Nₐ×1e-27
    IdealGasProblem(T, P, ρ₀, externalV, δv)
end


function (igp::IdealGasProblem)(_, flat_∂ψ, flat_ψ::Vector{Float64})
    try
    ψ = reshape(flat_ψ, size(igp.externalV))
    ρ = isnothing(flat_∂ψ) ? similar(ψ) : reshape(flat_∂ψ, size(igp.externalV))
    ρ .= exp.(ψ)
    ψ₀ = log(igp.ρ₀)
    value = igp.δv*sum(ρr*v + igp.T*(ρr*(ψr - ψ₀ - 1) + igp.ρ₀) for (ψr, v, ρr) in zip(ψ, igp.externalV, ρ))
    if !isnothing(flat_∂ψ)
        @. ρ *= igp.δv*(igp.externalV + igp.T*(ψ - ψ₀))
    end
    # isempty(flat_∂ψ) || @show extrema(flat_∂ψ)
    @show value, extrema(ψ), norm(vec(ρ))
    value
    catch e
        Base.showerror(stderr, e)
        Base.show_backtrace(stderr, catch_backtrace())
        rethrow()
    end
end

function mdft(igp::IdealGasProblem)
    ψ₀ = log(igp.ρ₀)
    return (ψ₀ .- vec(igp.externalV)./igp.T)
    # η = 0.0
    # η = 1e-5
    # ψ_init = (ψ₀ .- vec(igp.externalV)./igp.T) .* rand(1 .+ (-η:η), length(igp.externalV))
    # ψ_init = zeros(length(igp.externalV))
    # Optim.optimize(Optim.only_fg!(igp), ψ_init, Optim.LBFGS(), Optim.Options(iterations=20000))
end


struct SemiTruncatedInterpolator
    f::CubicSplineInterpolator{Float64, NoBoundaries}
    R::Float64
end
(f::SemiTruncatedInterpolator)(x) = ((@assert x ≥ 0); x > f.R ? 0.0 : f.f(x))
function SemiTruncatedInterpolator(qdht::QDHT, ĉ)
    f = CubicSplineInterpolator(qdht.r, (qdht\ĉ) ./ ((2π)^(3/2)), NoBoundaries())
    SemiTruncatedInterpolator(f, qdht.R)
end


function compute_dcf1D(tcf1D, R, ρ)
    n = length(tcf1D)
    qdht = QDHT{0,2}(R, n)
    hr = CubicSplineInterpolator([0; LinRange(0, R, n+1)[1:end-1] .+ (qdht.R/(2*n))],
                                 [-1.0; tcf1D], WeakBoundaries())
    hk = ((2π)^(3/2)) .* (qdht * hr.(qdht.r))
    ck = hk ./ (1 .+ ρ .* hk)
    qdht, ck
    # return qdht \ ck
end

## MonoAtomic

struct MonoAtomic <: MDFTProblem
    igp::IdealGasProblem
    ĉ₂::Array{ComplexF64,3} # ĉ₂(k), Fourier transform of the direct correlation function c₂(r)
    plan::FFTW.rFFTWPlan{Float64, -1, false, 3, UnitRange{Int64}}
end

function MonoAtomic(gasname::String, T::Float64, P::Float64, externalV::Array{Float64,3}, qdht, ĉ₂vec::Vector{Float64}, _mat::AbstractMatrix{Float64})
    c₂r = SemiTruncatedInterpolator(qdht, ĉ₂vec)
    a1, a2, a3 = size(externalV)
    mat = SMatrix{3,3,Float64,9}(_mat)
    δv = det(mat)/(a1*a2*a3)
    invmat = inv(mat)
    buffer, ortho, safemin = prepare_periodic_distance_computations(mat)
    stepmatrix = similar(buffer)
    c₂ = Array{Float64}(undef, a1, a2, a3)
    for i1 in 1:a1, i2 in 1:a2, i3 in 1:a3
        stepmatrix .= mat[:,1].*((i1-1)/a1) .+ mat[:,2].*((i2-1)/a2) .+ mat[:,3].*((i3-1)/a3)
        mul!(buffer, invmat, stepmatrix)
        c₂[i1,i2,i3] = δv*c₂r(periodic_distance!(buffer, mat, ortho, safemin))
    end
    plan = plan_rfft(c₂)
    ĉ₂ = plan * c₂
    # rĉ₂ = real.(ĉ₂)
    # @show maximum(abs.(imag.(ĉ₂)))
    # @assert ĉ₂ ≈ rĉ₂
    # @assert plan \ rĉ₂ ≈ c₂
    igp = IdealGasProblem(gasname, T, P, externalV, δv)
    MonoAtomic(igp, ĉ₂, plan)
end

function (ma::MonoAtomic)(_, flat_∂ψ, flat_ψ::Vector{Float64})
    ψ = reshape(flat_ψ, size(ma.igp.externalV))
    ρ = isnothing(flat_∂ψ) ? similar(ψ) : reshape(flat_∂ψ, size(ma.igp.externalV))
    ρ .= ma.igp.ρ₀ .* ψ.^2
    convol = ρ .- ma.igp.ρ₀
    rfftΔρ = ma.plan * convol
    rfftΔρ .*= ma.ĉ₂
    FFTW.ldiv!(convol, ma.plan, rfftΔρ)
    logρ = max.(log.(ρ ./ ma.igp.ρ₀), -1.3407807929942596e154)
    # for (logr, v, ρr, CΔρ, ψr) in zip(logρ, ma.igp.externalV, ρ, convol, ψ)
    #     if isinf(ρr*v + ma.igp.T*(ρr*logr + (ma.igp.ρ₀-ρr) - CΔρ*(ρr-ma.igp.ρ₀)/2))
    #         @show logr, v, ρr, CΔρ
    #         break
    #     elseif abs(2*ma.igp.ρ₀*ma.igp.δv*ψr*(v + ma.igp.T*(logr - CΔρ))) > 1e100
    #         @show logr, v, ρr, CΔρ, ψr
    #         break
    #     end
    # end
    logρ₀ = log(ma.igp.ρ₀)
    value = ma.igp.δv*sum(ρr*v + ma.igp.T*(ρr*logr - ρr*logρ₀ + (ma.igp.ρ₀-ρr) - CΔρ*(ρr-ma.igp.ρ₀)/2)
                for (logr, v, ρr, CΔρ) in zip(logρ, ma.igp.externalV, ρ, convol))
    # value = ifelse(isnan(value), Inf, value)
    if !isnothing(flat_∂ψ) # gradient update
        @. ρ = $(2*ma.igp.ρ₀*ma.igp.δv)*ψ*(ma.igp.externalV + ma.igp.T*(logρ - logρ₀ - convol))
        # @. ρ = ifelse(abs(ρ) > 1.3407807929942596e100, 0.0, ρ)
        @show value, maximum(ρ), extrema(ρ), norm(vec(ρ))
    end
    value
end


function mdft(ma::MonoAtomic, ψ_init=exp.(.-vec(ma.igp.externalV)./(2*ma.igp.T)))
    # # ρ = ρ₀*ψ^2
    # Optim.optimize(Optim.only_fg!(ma), ψ_init, Optim.LBFGS(linesearch=Optim.LineSearches.Static()), Optim.Options(iterations=100))
    # ψ_init .*= [1-0.1/i for i in 1:length(ψ_init)]
    # ρ_init = ma.igp.ρ₀ .* ψ_init.^2
    # convol = ρ .- ma.igp.ρ₀
    # rfftΔρ = ma.plan * convol
    # rfftΔρ .*= ma.ĉ₂
    # FFTW.ldiv!(convol, ma.plan, rfftΔρ)
    Optim.optimize(Optim.only_fg!(ma), ψ_init, Optim.LBFGS(linesearch=Optim.BackTracking(order=2, ρ_hi=0.3)),
                   Optim.Options(iterations=1000, f_tol=1e-10))
end

function finaldensity(ma::MonoAtomic, opt)
    ψ = Optim.minimizer(opt)
    ρ = ma.igp.ρ₀ .* ψ.^2
    sum(ρ)*ma.igp.δv
end
