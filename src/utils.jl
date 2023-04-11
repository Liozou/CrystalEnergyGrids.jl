# General utils

function wrap_atom(point, mat, invmat)
    abc = invmat * (point isa SVector ? point : SVector{3,Float64}(point))
    # return mat * (abc .- floor.(abc))
    x, y, z = mat * (abc .- floor.(abc))
    ε = 1e-14
    return SVector{3,Float64}(clamp(x, ε, mat[1,1]-ε), clamp(y, ε, mat[2,2]-ε), clamp(z, ε, mat[3,3]-ε))
end

function offsetpoint(point, mat, invmat, shift, size, dims)
    newpoint = wrap_atom(point, mat, invmat)
    @. (newpoint - shift)/size*dims + 1
end

function inverse_abc_offsetpoint(ipoint, invmat, shift, size, dims)
    abc = invmat*(@. (ipoint - 1)/dims*size + shift)
    abc .- floor.(abc)
end

# The following are copied from PeriodicGraphEmbeddings.jl

function cell_parameters(mat::AbstractMatrix)
    _a, _b, _c = eachcol(mat)
    a = norm(_a)
    b = norm(_b)
    c = norm(_c)
    α = acosd(_b'_c/(b*c))
    β = acosd(_c'_a/(c*a))
    γ = acosd(_a'_b/(a*b))
    (a, b, c), (α, β, γ)
end

function prepare_periodic_distance_computations(mat)
    (a, b, c), (α, β, γ) = cell_parameters(mat)
    ortho = all(x -> isapprox(Float16(x), 90; rtol=0.02), (α, β, γ))
    _a, _b, _c = eachcol(mat)
    safemin = min(Float64(dot(cross(_b, _c), _a)/(b*c)),
                  Float64(dot(cross(_c, _a), _b)/(a*c)),
                  Float64(dot(cross(_a, _b), _c)/(a*b)))/2
    # safemin is the half-distance between opposite planes of the unit cell
    return MVector{3,Float64}(undef), ortho, safemin
end

function periodic_distance!(buffer, mat, ortho, safemin)
    @simd for i in 1:3
        diff = buffer[i] + 0.5
        buffer[i] = diff - floor(diff) - 0.5
    end
    ref = norm(mat*buffer)
    (ortho || ref ≤ safemin) && return ref
    @inbounds for i in 1:3
        buffer[i] += 1
        newnorm = norm(mat*buffer)
        newnorm < ref && return newnorm # in a reduced lattice, there should be at most one
        buffer[i] -= 2
        newnorm = norm(mat*buffer)
        newnorm < ref && return newnorm
        buffer[i] += 1
    end
    return ref
end

function fraction_sites(egrid)
    bins = zeros(100)
    for val in egrid
        -10000 < val < 0 || continue
        bins[1+floor(Int, -val/100)] += 1
    end
    s = sum(bins)
    @show s, s/length(egrid)
    bins ./= s
    bins
end

function is_zaxis_linear(system)
    poss = position(system) ./ ANG_UNIT
    n = length(poss)
    xpos = poss[1][1]; ypos = poss[1][2]
    for (i, pos) in enumerate(poss)
        pos[1] ≈ xpos || return false
        pos[2] ≈ ypos || return false
    end
    xpos == 0 && ypos == 0 || @warn "Molecule linear along axis z should not be offset with respect to that axis (current offset: ($xpos, $ypos))."
    return true
end

function is_zaxis_linear_symmetric(system, islin)
    islin || return false
    poss = position(system) ./ ANG_UNIT
    n = length(poss)
    zpos = Vector{Float64}(undef, n)
    for (i, pos) in enumerate(poss)
        zpos[i] = pos[3]
    end
    I = sortperm(zpos)
    for i in 1:cld(n, 2)
        fwd = I[i]
        bwd = I[end-i+1]
        zpos[fwd] ≈ -zpos[bwd] || return false
        atomic_number(system, fwd) == atomic_number(system, bwd) || return false
    end
    return true
end

function zaxis_rotate(x2=rand(), x3=rand())
    s = sqrt(x3)
    v = SVector{3,Float64}(cospi(2*x2)*s, sinpi(2*x2)*s, sqrt(1-x3))
    2*v*v'-LinearAlgebra.I
end

function sphere_sample()
    v = randn(SVector{3,Float64})
    v ./ norm(v)
end

function add_big!(acc::BigFloat, x::BigFloat)
    @ccall Base.MPFR.libmpfr.mpfr_add(acc::Ref{BigFloat}, acc::Ref{BigFloat}, x::Ref{BigFloat}, Base.MPFR.ROUNDING_MODE[]::Base.MPFR.MPFRRoundingMode)::Int32
end
function div_big!(acc::BigFloat, x::BigFloat)
    @ccall Base.MPFR.libmpfr.mpfr_div(acc::Ref{BigFloat}, acc::Ref{BigFloat}, x::Ref{BigFloat}, Base.MPFR.ROUNDING_MODE[]::Base.MPFR.MPFRRoundingMode)::Int32
end
function muld_big!(acc::BigFloat, x::Cdouble)
    @ccall Base.MPFR.libmpfr.mpfr_mul_d(acc::Ref{BigFloat}, acc::Ref{BigFloat}, x::Cdouble, Base.MPFR.ROUNDING_MODE[]::Base.MPFR.MPFRRoundingMode)::Int32
end
function set_big!(acc::BigFloat, x::Cdouble)
    @ccall Base.MPFR.libmpfr.mpfr_set_d(acc::Ref{BigFloat}, x::Cdouble, Base.MPFR.ROUNDING_MODE[]::Base.MPFR.MPFRRoundingMode)::Int32
end
function exp_big!(acc::BigFloat, x::BigFloat)
    @ccall Base.MPFR.libmpfr.mpfr_exp(acc::Ref{BigFloat}, x::Ref{BigFloat}, Base.MPFR.ROUNDING_MODE[]::Base.MPFR.MPFRRoundingMode)::Int32
end


function meanBoltzmann(A, T)
    B = Array{Float64}(undef, size(A)[2:end]...)
    n = size(A, 1)
    # xs = Array{Float64}(undef, n)
    # τs = Array{BigFloat}(undef, n)
    Base.Threads.@threads for i in eachindex(IndexCartesian(), B)
        tmp = BigFloat()
        val = zero(BigFloat)
        tot = zero(BigFloat)
        for j in 1:n
            x = A[j,i]
            set_big!(tmp, -x/T)
            exp_big!(tmp, tmp)
            add_big!(tot, tmp)
            muld_big!(tmp, x)
            add_big!(val, tmp)
        end
        div_big!(val, tot)
        B[i] = Float64(val)
    end
    B
end
