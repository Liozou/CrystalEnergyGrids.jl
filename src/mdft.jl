import LinearAlgebra: norm, mul!
using Base.Threads
using Serialization

using StaticArrays
import Optim
using FFTW, Hankel
import BasicInterpolators: CubicSplineInterpolator, NoBoundaries, WeakBoundaries

FFTW.set_num_threads(nthreads()÷2)

export isotherm_hnc, IdealGasProblem, MonoAtomic

abstract type MDFTProblem <: Function end

struct IdealGasProblem <: MDFTProblem
    ρ₀::Float64   # reference number density, in number/Å³ (obtained from VdW coefficients)
    T::Float64    # temperature, in K
    P::Float64    # pressure, in bar
    externalV::Array{Float64,3} # external energy field in K
    mat::SMatrix{3,3,Float64,9} # unit cell
    δv::Float64   # volume of an elementary grid cube, in Å³, derived from externalV and mat
end
function IdealGasProblem(T::Float64, P::Float64, externalV::Array{Float64,3}, mat::SMatrix{3,3,Float64,9})
    IdealGasProblem(P/(0.0831446261815324*T), T, P, externalV, mat, det(mat)/length(externalV))
end
function IdealGasProblem(gasname::String, T::Float64, P::Float64, externalV::Array{Float64,3}, mat::SMatrix{3,3,Float64,9})
    a, b = get(VDW_COEFF, gasname, (0.0, 0.0))
    a == 0.0 && error(lazy"""Gas $gasname does not appear in the database, please enter an explicit ideal density or omit the name to rely on the perfect gas model. Possible gases are: $(join(keys(VDW_COEFF, ", "))). Please open an issue if you want to add the Van der Waals coefficient corresponding to your gas.""")
    c = 0.0831446261815324*T # perfect gas constants in bar*L/K/mol
    # The following is the real root of (Px^2 + a)(x-b) = c*x^2 given by Wolfram Alpha
    x = cbrt(sqrt((18*a*b*P^2 - 9*a*c*P + 2*b^3*P^3 + 6*b^2*c*P^2 + 6*b*c^2*P + 2*c^3)^2 + 4*(3*a*P - (-b*P - c)^2)^3) + 18*a*b*P^2 - 9*a*c*P + 2*b^3*P^3 + 6*b^2*c*P^2 + 6*b*c^2*P + 2*c^3)/(3*cbrt(2)*P) - (cbrt(2)*(3*a*P - (-b*P - c)^2))/(3*P*cbrt(sqrt((18*a*b*P^2 - 9*a*c*P + 2*b^3*P^3 + 6*b^2*c*P^2 + 6*b*c^2*P + 2*c^3)^2 + 4*(3*a*P - (-b*P - c)^2)^3) + 18*a*b*P^2 - 9*a*c*P + 2*b^3*P^3 + 6*b^2*c*P^2 + 6*b*c^2*P + 2*c^3)) - (-b*P - c)/(3*P)
    # inv(x) is ρ₀ in mol/L
    ρ₀ = inv(x) * 6.02214076e-4 # Nₐ×1e-27
    IdealGasProblem(ρ₀, T, P, externalV, mat, det(mat)/length(externalV))
end
function IdealGasProblem(ρ₀::Float64, T::Float64, P::Float64, externalV::Array{Float64,3}, mat::AbstractMatrix)
    IdealGasProblem(ρ₀, T, P, externalV, SMatrix{3,3,Float64,9}(mat), det(mat)/length(externalV))
end
function IdealGasProblem(gasname::String, T::Float64, P::Float64, externalV::Array{Float64,3}, mat::AbstractMatrix)
    IdealGasProblem(gasname, T, P, externalV, SMatrix{3,3,Float64,9}(mat))
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

function expand_correlation(c₂r, (a1, a2, a3)::NTuple{3,Int}, mat)
    δv = det(mat)/(a1*a2*a3)
    invmat = inv(mat)
    buffer, ortho, safemin = prepare_periodic_distance_computations(mat)
    stepmatrix = similar(buffer)
    c₂ = Array{Float64}(undef, a1, a2, a3)
    for i3 in 1:a3, i2 in 1:a2, i1 in 1:a1
        stepmatrix .= mat[:,1].*((i1-1)/a1) .+ mat[:,2].*((i2-1)/a2) .+ mat[:,3].*((i3-1)/a3)
        mul!(buffer, invmat, stepmatrix)
        c₂[i1,i2,i3] = δv*c₂r(periodic_distance!(buffer, mat, ortho, safemin))
    end
    c₂
end

struct MonoAtomic <: MDFTProblem
    igp::IdealGasProblem
    ĉ₂::Array{ComplexF64,3} # ĉ₂(k), Fourier transform of the direct correlation function c₂(r)
    plan::FFTW.rFFTWPlan{Float64, -1, false, 3, NTuple{3,Int}}
    c₂r::SemiTruncatedInterpolator
end

function MonoAtomic(gasname_or_ρ₀, T::Float64, P::Float64, externalV::Array{Float64,3}, c₂r::SemiTruncatedInterpolator, _mat::AbstractMatrix{Float64})
    mat = SMatrix{3,3,Float64,9}(_mat)
    c₂ = expand_correlation(c₂r, size(externalV), mat)
    plan = plan_rfft(c₂)
    ĉ₂ = plan * c₂
    # rĉ₂ = real.(ĉ₂)
    # @show maximum(abs.(imag.(ĉ₂)))
    # @assert ĉ₂ ≈ rĉ₂
    # @assert plan \ rĉ₂ ≈ c₂
    igp = IdealGasProblem(gasname_or_ρ₀, T, P, externalV, mat)
    MonoAtomic(igp, ĉ₂, plan, c₂r)
end

function MonoAtomic(gasname_or_ρ₀, T::Float64, P::Float64, externalV::Array{Float64,3}, qdht::QDHT{0,2,Float64}, ĉ₂vec::Vector{Float64}, _mat::AbstractMatrix{Float64})
    c₂r = SemiTruncatedInterpolator(qdht, ĉ₂vec)
    MonoAtomic(gasname_or_ρ₀, T, P, externalV, c₂r, _mat)
end

function (ma::MonoAtomic)(_, flat_∂ψ, flat_ψ::Vector{Float64})
    ψ = reshape(flat_ψ, size(ma.igp.externalV))
    ρ = isnothing(flat_∂ψ) ? similar(ψ) : reshape(flat_∂ψ, size(ma.igp.externalV))
    ρ .= ma.igp.ρ₀ .* ψ.^2
    convol = ρ .- ma.igp.ρ₀
    rfftΔρ = ma.plan * convol
    rfftΔρ .*= ma.ĉ₂
    FFTW.ldiv!(convol, ma.plan, rfftΔρ)
    logρmρ₀ = max.(log.(ρ ./ ma.igp.ρ₀), -1.3407807929942596e154)
    # for (logr, v, ρr, CΔρ, ψr) in zip(logρ, ma.igp.externalV, ρ, convol, ψ)
    #     if isinf(ρr*v + ma.igp.T*(ρr*logr + (ma.igp.ρ₀-ρr) - CΔρ*(ρr-ma.igp.ρ₀)/2))
    #         @show logr, v, ρr, CΔρ
    #         break
    #     elseif abs(2*ma.igp.ρ₀*ma.igp.δv*ψr*(v + ma.igp.T*(logr - CΔρ))) > 1e100
    #         @show logr, v, ρr, CΔρ, ψr
    #         break
    #     end
    # end
    value = ma.igp.δv*sum(ρr*v + ma.igp.T*(ρr*logrmr₀ + (ma.igp.ρ₀-ρr) - CΔρ*(ρr-ma.igp.ρ₀)/2)
                for (logrmr₀, v, ρr, CΔρ) in zip(logρmρ₀, ma.igp.externalV, ρ, convol))
    # value = ifelse(isnan(value), Inf, value)
    if !isnothing(flat_∂ψ) # gradient update
        @. ρ = $(2*ma.igp.ρ₀*ma.igp.δv)*ψ*(ma.igp.externalV + ma.igp.T*(logρmρ₀ - convol))
        # @. ρ = ifelse(abs(ρ) > 1.3407807929942596e100, 0.0, ρ)
        # @show value, maximum(ρ), extrema(ψ), norm(vec(ρ))
        # length(ma.igp.externalV) > 156248 && @show value, maximum(ρ)
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
    Optim.optimize(Optim.only_fg!(ma), ψ_init, Optim.LBFGS(linesearch=Optim.BackTracking(order=2, ρ_hi=0.1)),
                   Optim.Options(iterations=1000, f_tol=1e-10))
end


struct LinearMolecule <: MDFTProblem
    ρ₀::Float64   # reference number density, in number/Å³ (obtained from VdW coefficients)
    T::Float64    # temperature, in K
    P::Float64    # pressure, in bar
    externalV::Array{Float64,4} # external energy field in K (1st dimension: angles)
    lebedev_weights::Vector{Float64}
    mat::SMatrix{3,3,Float64,9} # unit cell
    δv::Float64   # volume of an elementary grid cube, in Å³, derived from externalV and mat
    ĉ₂::Array{ComplexF64,3} # ĉ₂(k), Fourier transform of the direct correlation function c₂(r)
    plan::FFTW.rFFTWPlan{Float64, -1, false, 3, NTuple{3,Int}}
    c₂r::SemiTruncatedInterpolator
end

function LinearMolecule(gasname_or_ρ₀, T::Float64, P::Float64, externalV::Array{Float64,4}, c₂r::SemiTruncatedInterpolator, _mat::AbstractMatrix{Float64})
    mat = SMatrix{3,3,Float64,9}(_mat)
    c₂ = expand_correlation(c₂r, size(externalV), mat)
    plan = plan_rfft(c₂)
    ĉ₂ = plan * c₂
    igp = IdealGasProblem(gasname_or_ρ₀, T, P, Array{Float64,3}(undef, size(externalV)[2:end]), mat)
    LinearMolecule(igp.ρ₀, igp.T, igp.P, externalV, lebedev_weights, mat, igp.δv, ĉ₂, plan, c₂r)
end
function LinearMolecule(gasname_or_ρ₀, T::Float64, P::Float64, externalV::Array{Float64,4}, qdht::QDHT{0,2,Float64}, ĉ₂vec::Vector{Float64}, _mat::AbstractMatrix{Float64})
    c₂r = SemiTruncatedInterpolator(qdht, ĉ₂vec)
    LinearMolecule(gasname_or_ρ₀, T, P, externalV, c₂r, _mat)
end

function (la::LinearMolecule)(_, flat_∂ψ, flat_ψ::Vector{Float64})
    a0, a1, a2, a3 = size(la.externalV)
    ρ₀π = la.ρ₀ / (4π)
    ψ = reshape(flat_ψ, a0, a1, a2, a3)
    ρ_average = Array{Float64}(undef, a1, a2, a3)
    ρ = isnothing(flat_∂ψ) ? similar(ψ) : reshape(flat_∂ψ, size(ma.igp.externalV))
    logρmρ₀ = similar(ψ)
    Fext_contrib = Vector{Float64}(undef, a3)
    Fid_contrib = Vector{Float64}(undef, a3)
    @threads for i3 in 1:a3
        fext_c = 0.0
        fid_c = 0.0
        for i2 in 1:a2, i1 in 1:a1
            tot = 0.0
            for i0 in 1:a0
                ψ2 = ψ[i0,i1,i2,i3]^2
                logrmr₀ = log(ψ2) # log(ρ / ρ₀)
                logρmρ₀[i0,i1,i2,i3] = logrmr₀
                wψ2 = la.lebedev_weights[i0]*ψ2
                tot += wψ2
                fext_c += wψ2*la.externalV[i0,i1,i2,i3]
                fid_c += la.lebedev_weights[i0]*(ρ₀π*ψ2*(logrmr₀ - 1) + ρ₀π)
            end
            ρ_average[i1,i2,i3] = tot*la.ρ₀
        end
        Fext_contrib[i3] = fext_c
        Fid_contrib[i3] = fid_c
    end
    Fext = la.ρ₀*sum(Fext_contrib)
    convol = ρ_average .- la.ρ₀
    rfftΔρ = la.plan * convol
    rfftΔρ .*= la.ĉ₂
    FFTW.ldiv!(convol, la.plan, rfftΔρ)

    # ρ .= ρ₀π .* ψ.^2
    # logρmρ₀ = max.(log.(ρ ./ ρ₀π), -1.3407807929942596e154)
    # Fid = la.T*sum(ρr*(logr - 1) + ρ₀π for (ρr, logr) in zip(ρ, logρmρ₀))
    Fexc = -sum(CΔρ*(ρ_ave-la.ρ₀) for (CΔρ, ρ_ave) in zip(convol, ρ_average))/2
    if !isnothing(flat_∂ψ) # gradient update
        @. ρ = $(2*ρ₀π*la.δv)*ψ
        for i0 in 1:a0
            @. ρ[i0,:,:,:] *= la.externalV[i0,:,:,:] + la.T*(logρmρ₀*la.lebedev_weights[i0] - convol)
        end
        # @. ρ = ifelse(abs(ρ) > 1.3407807929942596e100, 0.0, ρ)
        # @show value, maximum(ρ), extrema(ψ), norm(vec(ρ))
    end
    (Fid + Fext + Fexc)*la.δv
end


function mdft(la::LinearMolecule, ψ_init=exp.(.-vec(la.externalV)./(2*la.T)))
    Optim.optimize(Optim.only_fg!(pa), ψ_init, Optim.LBFGS(linesearch=Optim.BackTracking(order=2, ρ_hi=0.1)),
                   Optim.Options(iterations=1000, f_tol=1e-10))
end



function __inter(i, a, b)
    m = a÷b
    (1+(i-1)*m):(i == b ? a : i*m)
end

function iterative_mdft(ma::MonoAtomic)
    rψ_init = exp.(.-ma.igp.externalV./(2*ma.igp.T))
    a1, a2, a3 = size(ma.igp.externalV)
    b1, b2, b3 = (a1, a2, a3) .÷ 20
    ψ = Vector{Float64}(undef, b1*b2*b3)
    rψ = reshape(ψ, b1, b2, b3)
    @threads for i3 in 1:b3
        for i2 in 1:b2, i1 in 1:b1
            rψ[i1,i2,i3] = mean(@view(rψ_init[__inter(i1, a1, b1), __inter(i2, a2, b2), __inter(i3, a3, b3)]))
        end
    end
    opt = nothing
    for M in (20, 15, 11, 8, 5, 3, 1)
        @show M
        c1, c2, c3 = size(rψ)
        (b1, b2, b3) = (a1, a2, a3) .÷ M
        newψ = Vector{Float64}(undef, b1*b2*b3)
        newrψ = reshape(newψ, b1, b2, b3)
        @threads for i3 in 1:c3
            for i2 in 1:c2, i1 in 1:c1
                J1 = __inter(i1, b1, c1)
                J2 = __inter(i2, b2, c2)
                J3 = __inter(i3, b3, c3)
                for j3 in J3, j2 in J2, j1 in J1
                    newrψ[j1,j2,j3] = rψ[i1,i2,i3]
                end
            end
        end
        newexternalV = Array{Float64}(undef, b1, b2, b3)
        @threads for i3 in 1:b3
            for i2 in 1:b2, i1 in 1:b1
                newexternalV[i1,i2,i3] = mean(@view(ma.igp.externalV[__inter(i1, a1, b1), __inter(i2, a2, b2), __inter(i3, a3, b3)]))
            end
        end
        newma = MonoAtomic(ma.igp.ρ₀, ma.igp.T, ma.igp.P, newexternalV, ma.c₂r, ma.igp.mat)
        ψ, newψ = newψ, ψ
        opt = Optim.optimize(Optim.only_fg!(newma), ψ, Optim.LBFGS(linesearch=Optim.BackTracking(order=2, ρ_hi=0.1)),
                   Optim.Options(iterations=10000, f_tol=1e-10))
        display(opt)
        @show finaldensity(newma, opt)
    end
    opt
end



function finaldensity(ma::MonoAtomic, opt)
    if !Optim.converged(opt)
        @error "Optimizer failed to converge; proceeeding with partial result"
    end
    ψ = Optim.minimizer(opt)
    ρ = ma.igp.ρ₀ .* ψ.^2
    sum(ρ)*ma.igp.δv
end


function compute_store_hnc!(potential, ffname, T, ρ₀, qdht)
    store = "$ffname-$(hash(potential))-$T-$ρ₀-$(qdht.R)-$(qdht.N)"
    global scratchspace
    path = joinpath(scratchspace, store)
    if isfile(path)
        stored = deserialize(path)
        result = hnc(potential, qdht.R, T, ρ₀; qdht, γ0=copy(stored[2]))
        if result != stored
            @warn "Previously stored hnc computation stale for FF $ffname (T = $T, ρ₀ = $ρ₀); refreshing."
            serialize(path, result)
        end
        return result
    end
    ret = hnc(potential, qdht.R, T, ρ₀; qdht)
    serialize(path, ret)
    return ret
end

function isotherm_hnc(ff::ForceField, mol::AbstractSystem, egrid, temperature, pressures, mat;
                      molname=identify_molecule(atomic_symbol(mol)), qdht=QDHT{0,2}(100, 1000))
    potential = compute_average_self_potential(mol, ff, qdht.r)
    m = length(pressures)
    ck_hnc = Vector{Vector{Float64}}(undef, m)
    @threads for i in 1:m
        P = pressures[i]
        igp = IdealGasProblem(molname, temperature, P, egrid, mat)
        c, _ = compute_store_hnc!(potential, ff.name, temperature, igp.ρ₀, qdht)
        ck_hnc[i] = qdht * c
    end
    isotherm = Vector{Float64}(undef, m)
    for i in 1:m
        ck = ck_hnc[i]
        P = pressures[i]
        ma = MonoAtomic(molname, temperature, P, egrid, qdht, ck, mat)
        opt = mdft(ma)
        isotherm[i] = finaldensity(ma, opt)
        @show P, isotherm[i]
        display(opt)
    end
    isotherm
end
