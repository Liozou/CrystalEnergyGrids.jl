# General utils

function wrap_atom(point, mat, invmat)
    abc = invmat * (point isa SVector ? point : SVector{3,Float64}(point))
    abc = map(x -> ifelse(-eps() < x < eps(), 0.0, x), abc)
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

function periodic_distance(buffer, mat)
    _, ortho, safemin = prepare_periodic_distance_computations(mat)
    periodic_distance!(buffer, mat, ortho, safemin)
end

"""
    periodic_distance_with_ofs!(buffer, ofs, mat, ortho, safemin)

Equivalent to `periodic_distance!` but also sets the corresponding integer offset in `ofs`.

If `buffer == inv(mat)*(pA - pB)` where `pA` and `pB` are the positions of two atoms A and
B, then the resulting periodic distance will be equal to `norm(pA - (pB + mat*ofs))`, i.e.
the (aperiodic) distance between atom A in cell (0,0,0) and the image of B in cell `ofs`.
"""
function periodic_distance_with_ofs!(buffer, ofs, mat, ortho, safemin)
    @simd for i in 1:3
        diff = buffer[i] + 0.5
        buffer[i] = diff - begin ofs[i] = floor(Int, diff) end - 0.5
    end
    ref = norm(mat*buffer)
    (ortho || ref ≤ safemin) && return ref
    @inbounds for i in 1:3
        buffer[i] += 1
        newnorm = norm(mat*buffer)
        if newnorm < ref
            ofs[i] -= 1
            return newnorm # in a reduced lattice, there should be at most one
        end
        buffer[i] -= 2
        ofs[i] += 2
        newnorm = norm(mat*buffer)
        if newnorm < ref
            ofs[i] += 1
            return newnorm # in a reduced lattice, there should be at most one
        end
        buffer[i] += 1
    end
    return ref
end

# Root mean-square deviation i.e. √[(∑(xᵢ-yᵢ)²)/n]
function rmsd(a::AbstractVector{T}, b::AbstractVector{T}) where T
    ret = zero(T)
    n = length(a)
    @assert eachindex(a) == eachindex(b)
    for (x, y) in zip(a, b)
        # if !isfinite(x)
        #     isfinite(y) && return abs(x)
        #     continue
        # elseif !isfinite(y)
        #     isfinite(x) && return abs(y)
        #     continue
        # end
        ret += (x - y)^2
    end
    sqrt(ret / n)
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
    poss = position(system) / u"Å"
    xpos = poss[1][1]; ypos = poss[1][2]
    for pos in poss
        pos[1] ≈ xpos || return false
        pos[2] ≈ ypos || return false
    end
    xpos == 0 && ypos == 0 || @warn "Molecule linear along axis z should not be offset with respect to that axis (current offset: ($xpos, $ypos))."
    return true
end

function is_zaxis_linear_symmetric(system, islin)
    islin || return false
    poss = position(system) / u"Å"
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
        AtomsBase.atomic_number(system, fwd) == AtomsBase.atomic_number(system, bwd) || return false
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

# function add_big!(acc::BigFloat, x::BigFloat)
#     @ccall Base.MPFR.libmpfr.mpfr_add(acc::Ref{BigFloat}, acc::Ref{BigFloat}, x::Ref{BigFloat}, Base.MPFR.ROUNDING_MODE[]::Base.MPFR.MPFRRoundingMode)::Int32
# end
# function div_big!(acc::BigFloat, x::BigFloat)
#     @ccall Base.MPFR.libmpfr.mpfr_div(acc::Ref{BigFloat}, acc::Ref{BigFloat}, x::Ref{BigFloat}, Base.MPFR.ROUNDING_MODE[]::Base.MPFR.MPFRRoundingMode)::Int32
# end
# function muld_big!(acc::BigFloat, x::Cdouble)
#     @ccall Base.MPFR.libmpfr.mpfr_mul_d(acc::Ref{BigFloat}, acc::Ref{BigFloat}, x::Cdouble, Base.MPFR.ROUNDING_MODE[]::Base.MPFR.MPFRRoundingMode)::Int32
# end
# function set_big!(acc::BigFloat, x::Cdouble)
#     @ccall Base.MPFR.libmpfr.mpfr_set_d(acc::Ref{BigFloat}, x::Cdouble, Base.MPFR.ROUNDING_MODE[]::Base.MPFR.MPFRRoundingMode)::Int32
# end
# function exp_big!(acc::BigFloat, x::BigFloat)
#     @ccall Base.MPFR.libmpfr.mpfr_exp(acc::Ref{BigFloat}, x::Ref{BigFloat}, Base.MPFR.ROUNDING_MODE[]::Base.MPFR.MPFRRoundingMode)::Int32
# end


# function meanBoltzmannBig(A, T)
#     B = Array{Float64}(undef, size(A)[2:end]...)
#     n = size(A, 1)
#     # xs = Array{Float64}(undef, n)
#     # τs = Array{BigFloat}(undef, n)
#     Base.Threads.@threads for i in eachindex(IndexCartesian(), B)
#         tmp = BigFloat()
#         val = zero(BigFloat)
#         tot = zero(BigFloat)
#         for j in 1:n
#             x = A[j,i]
#             set_big!(tmp, -x/T)
#             exp_big!(tmp, tmp)
#             add_big!(tot, tmp)
#             muld_big!(tmp, x)
#             add_big!(val, tmp)
#         end
#         div_big!(val, tot)
#         B[i] = Float64(val)
#     end
#     B
# end

function meanBoltzmann(A, T)
    n = size(A, 1)
    if ndims(A) == 1
        facts = 0.0
        tot = 0.0
        M = minimum(A) - 30.0*T
        for j in 1:n
            y = A[j]
            fact = exp((M-y)/T)
            facts += fact
            tot += fact*y
        end
        return tot/facts
    end
    B = Array{Float64}(undef, size(A)[2:end]...)
    Base.Threads.@threads for i in eachindex(IndexCartesian(), B)
        factors = 0.0
        total = 0.0
        m = minimum(view(A, 1:n, i)) - 30.0*T # exp(30.0) ≈ 1e13
        for j in 1:n
            x = A[j,i]
            factor = exp((m-x)/T)
            factors += factor
            total += factor*x
        end
        B[i] = total / factors
    end
    return B
end


function downsize(X::Array{T,3}, n1, n2, n3) where {T}
    [mean(@view X[1+(i1-1)*n1:i1*n1, 1+(i2-1)*n2:i2*n2, 1+(i3-1)*n3:i3*n3])
     for i1 in 1:(size(X,1)÷n1), i2 in 1:(size(X,2)÷n2), i3 in 1:(size(X,3)÷n3)]
end
downsize(X::Array{T,3}, n) where {T} = downsize(X, n, n, n)


# Accurately Computing log(1 − exp(− |a|)) Assessed by the Rmpfr package by Martin Mächler
log1pexp(x::Float64) = x <= -36.7368005696771 ? exp(x) : x < 18.021826694558577 ? log1p(exp(x)) : x < 33.23111882352963 ? x + exp(-x) : x
# log1mexp(x) = x < log(2) ? log(-expm1(-x)) : log1p(-exp(-x))
# LogExpFunctions.jl for coefficient fine-tuning
logexpm1(x::Float64) = x < 18.021826694558577 ? log(expm1(x)) : x < 33.23111882352963 ? x - exp(-x) : x
logistic(x::Float64) = x < -744.4400719213812 ? 0.0 : x < 36.7368005696771 ? (e = exp(x); e/(1.0 + e)) : 1.0

const ALLATOMS = Set(["D", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", "Uue"])

function identify_atom(symb)
    s = titlecase(split(String(symb), x -> x === 'w' || !isletter(x); limit=2)[1]; strict=true)
    while s != "" && s ∉ ALLATOMS
        s = s[1:prevind(s, end)]
    end
    s
end

"""
    identify_molecule(atomsymbols)

Given `atomsymbols`, a list of `Symbol`s corresponding to atoms of a molecule, attempts to
retrieve the name of the molecule.

Examples:
```jldoctest
julia> CrystalEnergyGrids.identify_molecule([:Oa, :C_co2, :Ob])
"CO2"

julia> CrystalEnergyGris.identify_molecule([:Hw, :O5, :Hz])
"H2O"
```
"""
function identify_molecule(atomsymbols)
    symbs = unique!(sort(atomsymbols))
    atoms = identify_atom.(symbs)
    I = sortperm(atoms)
    n = length(symbs)
    vmap = Vector{Int}(undef, n)
    last_symbol = Symbol("")
    m = 0
    for i in 1:n
        j = I[i]
        if atoms[j] != last_symbol
            last_symbol = atoms[j]
            m += 1
        end
        vmap[j] = m
    end
    unique!(atoms)
    dict = Dict(s => vmap[i] for (i,s) in enumerate(symbs))
    nums = [0 for _ in 1:m]
    for s in atomsymbols
        nums[dict[s]] += 1
    end
    parts = String[]
    for i in 1:m
        push!(parts, atoms[i])
        if nums[i] > 1
            push!(parts, string(nums[i]))
        end
    end
    join(parts)
end
