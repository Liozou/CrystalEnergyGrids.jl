import Optim
using FFTW, Hankel
import BasicInterpolators: CubicSplineInterpolator, NoBoundaries
import LinearAlgebra: norm

export IdealGasProblem, MonoAtomic

abstract type MDFTProblem <: Function end

struct IdealGasProblem <: MDFTProblem
    T::Float64    # in K
    P::Float64    # in bar
    ρ₀::Float64   # in number/Å³
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
    # inv(x) is ρ₀ in mol/L
    ρ₀ = inv(x) * 6.02214076e-4 # Nₐ×1e-27
    IdealGasProblem(T, P, ρ₀, externalV)
end


function (igp::IdealGasProblem)(_, flat_∂ψ, flat_ψ::Vector{Float64})
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
    return (ψ₀ .- vec(igp.externalV)./igp.T)
    # η = 1e-9
    # ψ_init = (ψ₀ .- vec(igp.externalV)./igp.T)# .* rand(1 .+ (-η:η), length(igp.externalV))
    # Optim.optimize(Optim.only_fg!(igp), ψ_init, Optim.LBFGS(), Optim.Options(iterations=100))
end

function compute_dcf1D(tcf1D, R, ρ)
    n = length(tcf1D)
    qdht = QDHT{0,1}(R, n)
    hr = CubicSplineInterpolator(LinRange(0, qdht.R, n), tcf1D)
    hk = qdht * hr.(qdht.r)
    ck = hk ./ (1 .+ ρ .* hk)
    qdht, ck
    # return qdht \ ck
end
    # K = qdht.K
    # hk = CubicSplineInterpolator(qdht.k, qdht * hr.(qdht.r), NoBoundaries())
    # ck3D = Array{Float64,3}(undef, 2m+1, 2m+1, 2m+1)
    # for i1 in -m:m, i2 in -m:m, i3 in -m:m
    #     k = hypot(i1/m*K, i2/m*K, i3/m*K)
    #     val = k < K ? hk(k) : 0.0
    #     ck3D[i1+m+1, i1+m+1, i3+m+1] = val/(1+ρ*val)
    # end
    # return ifft(ck3D)
# end

struct MonoAtomic0 <: MDFTProblem
    igp::IdealGasProblem
    # rmax::Float64
    # kmax::Float64
    ĉ₂::Array{Float64,3} # ĉ₂(k), Fourier transform of the direct correlation function c₂(r)
end

MonoAtomic = MonoAtomic0
function MonoAtomic0(gasname::String, T::Float64, P::Float64, externalV::Array{Float64,3}, qdht, ĉ₂vec::Vector{Float64})
    igp = IdealGasProblem(gasname, T, P, externalV)
    ĉ₂k = CubicSplineInterpolator([0.0; qdht.k], [onaxis(ĉ₂vec, qdht); ĉ₂vec])#, NoBoundaries())
    a1, a2, a3 = size(externalV)
    ks1 = rfftfreq(a1); ks2 = fftfreq(a2); ks3 = fftfreq(a3)
    a1 = length(ks1)
    ĉ₂ = Array{Float64,3}(undef, a1, a2, a3)
    for i1 in 1:a1, i2 in 1:a2, i3 in 1:a3
        k = hypot(ks1[i1], ks2[i2], ks3[i3])
        # val = k < K ? ĉ₂k(k) : 0.0
        ĉ₂[i1,i2,i3] = ĉ₂k(k)
    end
    MonoAtomic(igp, ĉ₂)
end

function (ma::MonoAtomic)(_, flat_∂ψ, flat_ψ::Vector{Float64})
    ψ = reshape(flat_ψ, size(igp.externalV))
    gradient = isnothing(flat_∂ψ) ? similar(ψ) : reshape(flat_∂ψ, size(igp.externalV))
    gradient .= exp.(ψ)
    ψ₀ = log(igp.ρ₀)
    value = sum(expr*v + igp.T*(expr*(r - ψ₀ - 1) + igp.ρ₀) for (r, v, expr) in zip(ψ, igp.externalV, gradient))
    if !isnothing(flat_∂ψ)
        @. gradient *= igp.externalV + igp.T*(ψ - ψ₀)
    end
    ψ = reshape(flat_ψ, size(igp.externalV))
    ρ .= exp.(ψ)
    return value
end

function mdft(ma::MonoAtomic)
    ψ₀ = log(ma.igp.ρ₀)
    ψ_init = (ψ₀ .- vec(ma.igp.externalV)./ma.igp.T)
    Optim.optimize(Optim.only_fg!(ma), ψ_init, Optim.LBFGS(), Optim.Options(iterations=100))
end
