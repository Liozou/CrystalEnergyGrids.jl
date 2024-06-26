# Types

const TÅ = typeof(1.0u"Å")
const TK = typeof(1.0u"K")
const Te_au = typeof(1.0u"e_au")


# General utils

function perpendicular_lengths(a, b, c)
    axb1 = a[2]*b[3] - a[3]*b[2]
    axb2 = a[3]*b[1] - a[1]*b[3]
    axb3 = a[1]*b[2] - a[2]*b[1]
    axb = SVector{3}((axb1, axb2, axb3))
    bxc1 = b[2]*c[3] - b[3]*c[2]
    bxc2 = b[3]*c[1] - b[1]*c[3]
    bxc3 = b[1]*c[2] - b[2]*c[1]
    bxc = SVector{3}((bxc1, bxc2, bxc3))
    cxa1 = c[2]*a[3] - c[3]*a[2]
    cxa2 = c[3]*a[1] - c[1]*a[3]
    cxa3 = c[1]*a[2] - c[2]*a[1]
    cxa = SVector{3}((cxa1, cxa2, cxa3))
    volume = abs(dot(a, bxc))
    cx = volume/norm(bxc)
    cy = volume/norm(cxa)
    cz = volume/norm(axb)
    SVector{3}((cx, cy, cz))
end
perpendicular_lengths(mat) = perpendicular_lengths(@view(mat[:,1]), @view(mat[:,2]), @view(mat[:,3]))

"""
    find_supercell(mat::AbstractMatrix, cutoff)
    find_supercell((a, b, c)::NTuple{3,<:AbstractVector}, cutoff)
    find_supercell((cx, cy, cz)::NTuple{3,AbstractFloat}, cutoff)
    find_supercell(syst::AbstractSystem{3}, cutoff)

Return the triplet `(ix, iy, iz)` such that the `ix × iy × iz` supercell of the input has
perpendicular widths each at least twice as long as the cutoff.

The input can be the unit cell matrix `mat`, its three axes `(a, b, c)` or its three
perpendicular widths `(cx, cy, cz)`.
"""
function find_supercell((cx, cy, cz)::NTuple{3,T}, cutoff::T) where T
    ceil(Int, 2cutoff/cx), ceil(Int, 2cutoff/cy), ceil(Int, 2*cutoff/cz)
end
function find_supercell((ux, uy, uz)::NTuple{3,T}, cutoff::Real) where {T<:Quantity}
    cx, _cx = promote(ux, cutoff)
    cy, _cy = promote(uy, cutoff)
    cz, _cz = promote(uz, cutoff)
    @assert _cx == _cy == _cz
    find_supercell((cx, cy, cz), _cx * unit(T))
end
function find_supercell((ux, uy, uz)::NTuple{3,T}, cutoff::Quantity) where {T<:Real}
    cx, _cx = promote(ux, cutoff)
    cy, _cy = promote(uy, cutoff)
    cz, _cz = promote(uz, cutoff)
    @assert _cx == _cy == _cz
    U = unit(cutoff)
    find_supercell((cx*U, cy*U, cz*U), _cx)
end
find_supercell(mat::AbstractMatrix, cutoff) = find_supercell((view(mat, :, 1), view(mat, :, 2), view(mat, :, 3)), cutoff)
find_supercell((a, b, c)::NTuple{3,<:AbstractVector}, cutoff) = find_supercell(perpendicular_lengths(a, b, c), cutoff)
find_supercell(syst::AbstractSystem{3}, cutoff) = find_supercell(bounding_box(syst), cutoff)
function find_supercell(x, cutoff)
    (x isa NTuple || length(x) != 3) && throw(MethodError(find_supercell, (x, cutoff)))
    find_supercell((x[1], x[2], x[3]), cutoff)
end

function stack3((a, b, c)::NTuple{3,<:AbstractVector{<:Real}})
    SMatrix{3,3,Float64,9}((a[1], a[2], a[3], b[1], b[2], b[3], c[1], c[2], c[3]))
end
function stack3((a, b, c)::NTuple{3,<:AbstractVector{<:(Quantity{U,Unitful.𝐋} where U)}})
    SMatrix{3,3,Float64,9}((NoUnits(a[1]/u"Å"), NoUnits(a[2]/u"Å"), NoUnits(a[3]/u"Å"),
                            NoUnits(b[1]/u"Å"), NoUnits(b[2]/u"Å"), NoUnits(b[3]/u"Å"),
                            NoUnits(c[1]/u"Å"), NoUnits(c[2]/u"Å"), NoUnits(c[3]/u"Å")))
end
function stack3(x)
    x isa NTuple && throw(MethodError(stack3, (x,)))
    stack3((x[1], x[2], x[3]))
end

@static if VERSION < v"1.9-"
    stack(x) = reduce(hcat, x)
end

function make_supercell(s::AbstractSystem{3}, (ΠA, ΠB, ΠC))
    AtomsBase.boundary_conditions(s) == SA[Periodic(), Periodic(), Periodic()] || error("Cannot make a supercell of a non 3D-periodic system.")
    Π = ΠA*ΠB*ΠC
    symbols = atomic_symbol(s)
    n = length(symbols)
    for i in 1:n
        spls = rsplit(String(symbols[i]), '_'; limit=2)
        if length(spls) == 2 && all(isdigit, spls[2])
            symbols[i] = Symbol(spls[1])
        end
    end
    refpositions = position(s)
    @assert length(refpositions) == n
    newpositions = Vector{SVector{3,TÅ}}(undef, n*Π)
    newsymbols = Vector{Symbol}(undef, n*Π)
    A, B, C = bounding_box(s)
    k = 0
    for iA in 0:(ΠA-1)
        vA = iA * A
        for iB in 0:(ΠB-1)
            vB = iB * B + vA
            for iC in 0:(ΠC-1)
                vC = iC * C + vB
                for (i, (p, s)) in enumerate(zip(refpositions, symbols))
                    newpositions[k+i] = p + vC
                    newsymbols[k+i] = Symbol(s, :_, k+i)
                end
                k += n
            end
        end
    end
    FastSystem(
        bounding_box(s) .* (ΠA, ΠB, ΠC),
        (Periodic(), Periodic(), Periodic()),
        newpositions,
        newsymbols,
        repeat(atomic_number(s), Π),
        repeat(atomic_mass(s), Π),
    )
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

function mat_from_parameters((a, b, c), (α, β, γ))
    cosα = cosd(α); cosβ = cosd(β); sinγ, cosγ = sincosd(γ)
    ω = sqrt(1 - cosα^2 - cosβ^2 - cosγ^2 + 2*cosα*cosβ*cosγ)
    SMatrix{3,3,eltype(a),9}(a, 0, 0, b*cosγ, b*sinγ, 0, c*cosβ, c*(cosα - cosβ*cosγ)/sinγ, c*ω/sinγ)
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


function norm2(u::S, v::T) where {S,T}
    r2 = zero(eltype(S))*zero(eltype(T)) # == zero(u)^2 == zero(v)^2
    rx = eachindex(u)
    @simd for i in rx
        r2 += (u[i] - v[i])^2
    end
    r2
end
function norm2(u::T) where {T}
    r2 = zero(eltype(T))^2
    @simd for x in u
        r2 += x^2
    end
    r2
end

"""
    periodic_distance2_fromcartesian!(buffer, mat, invmat, ortho, safemin2, buffer2)

Similar to `periodic_distance!` except that:
- `buffer` contains cartesian coordinates, not fractional ones.
- `safemin2` should be `safemin^2`
- `buffer2` should be a `MVector{3,Float64}`.
- the result is the squared distance.

After this function returns the squared distance, `buffer` contains the cartesian
coordinates corresponding to the closest image.
"""
function periodic_distance2_fromcartesian!(buffer, mat, invmat, ortho, safemin2, buffer2)
    mul!(buffer2, invmat, buffer)
    periodic_distance2!(buffer, mat, ortho, safemin2, buffer2)
end

"""
    periodic_distance2!(buffer, mat, ortho, safemin2, buffer2)

Similar to `periodic_distance!` except that:
- the input fractional coordinates are stored in `buffer2`
- `safemin2` should be `safemin^2`
- the result is the squared distance.

After this function returns the squared distance, `buffer` contains the cartesian
coordinates corresponding to the closest image.
"""
function periodic_distance2!(buffer, mat, ortho, safemin2, buffer2)
    @simd for i in 1:3
        diff = buffer2[i] + 0.5
        buffer2[i] = diff - floor(diff) - 0.5
    end
    mul!(buffer, mat, buffer2)
    ref2 = norm2(buffer)
    (ortho || ref2 ≤ safemin2) && return ref2
    @inbounds for i in 1:3
        buffer2[i] += 1
        mul!(buffer, mat, buffer2)
        newnorm2 = norm2(buffer)
        newnorm2 < ref2 && return newnorm2 # in a reduced lattice, there should be at most one
        buffer2[i] -= 2
        mul!(buffer, mat, buffer2)
        newnorm2 = norm2(buffer)
        newnorm2 < ref2 && return newnorm2
        buffer2[i] += 1
    end
    return ref2
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
        newnorm = norm(mat*buffer)
        if newnorm < ref
            ofs[i] += 1
            return newnorm # in a reduced lattice, there should be at most one
        end
        buffer[i] += 1
    end
    return ref
end

struct CellMatrix
    mat::SMatrix{3,3,TÅ,9}
    invmat::SMatrix{3,3,typeof(1.0u"Å^-1"),9}
end
CellMatrix(mat::AbstractMatrix) = CellMatrix(mat, inv(ustrip.(u"Å", mat))*u"Å^-1")
CellMatrix(mat::AbstractMatrix{Float64}) = CellMatrix(mat*u"Å")
CellMatrix(system::AbstractSystem{3}) = CellMatrix(SMatrix{3,3,TÅ,9}(stack(bounding_box(system))))
function CellMatrix(mat::AbstractMatrix{Float64}, invmat::AbstractMatrix{Float64})
    CellMatrix(mat*u"Å", invmat*u"Å^-1")
end
CellMatrix() = CellMatrix(SMatrix{3,3}([Inf 0.0 0.0; 0.0 Inf 0.0; 0.0 0.0 Inf]u"Å"))

function unsafe_periodic_distance2!(buffer, buffer2, cell::CellMatrix)
    mul!(buffer2, cell.invmat, buffer)
    @simd for i in 1:3
        diff = buffer2[i] + 0.5
        buffer2[i] = diff - floor(diff) - 0.5
    end
    mul!(buffer, cell.mat, buffer2)
    norm2(buffer)
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
    poss = map.(Base.Fix1(ustrip, u"Å"), position(system))
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
        AtomsBase.atomic_symbol(system, fwd) == AtomsBase.atomic_symbol(system, bwd) || return false
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

function meanBoltzmann(A, T, weights=nothing)
    n = size(A, 1)
    if ndims(A) == 1
        facts = 0.0
        tot = 0.0
        M = minimum(A) - 30.0*T
        for j in 1:n
            y = A[j]
            fact = exp((M-y)/T)*(weights isa Nothing ? true : weights[j])
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
            factor = exp((m-x)/T)*(weights isa Nothing ? true : weights[j])
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

const ATOMS_PER_NUM = [:H, :He, :Li, :Be, :B, :C, :N, :O, :F, :Ne, :Na, :Mg, :Al, :Si, :P, :S, :Cl, :Ar, :K, :Ca, :Sc, :Ti, :V, :Cr, :Mn, :Fe, :Co, :Ni, :Cu, :Zn, :Ga, :Ge, :As, :Se, :Br, :Kr, :Rb, :Sr, :Y, :Zr, :Nb, :Mo, :Tc, :Ru, :Rh, :Pd, :Ag, :Cd, :In, :Sn, :Sb, :Te, :I, :Xe, :Cs, :Ba, :La, :Ce, :Pr, :Nd, :Pm, :Sm, :Eu, :Gd, :Tb, :Dy, :Ho, :Er, :Tm, :Yb, :Lu, :Hf, :Ta, :W, :Re, :Os, :Ir, :Pt, :Au, :Hg, :Tl, :Pb, :Bi, :Po, :At, :Rn, :Fr, :Ra, :Ac, :Th, :Pa, :U, :Np, :Pu, :Am, :Cm, :Bk, :Cf, :Es, :Fm, :Md, :No, :Lr, :Rf, :Db, :Sg, :Bh, :Hs, :Mt, :Ds, :Rg, :Cn, :Nh, :Fl, :Mc, :Lv, :Ts, :Og, :Uue]
const ALLATOMS = Set(["D", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", "Uue"])

function identify_atom(symb)
    s = titlecase(split(String(symb), x -> x === 'w' || !isletter(x); limit=2)[1]; strict=true)
    while s != "" && s ∉ ALLATOMS
        s = s[1:prevind(s, end)]
    end
    s
end

"""
    identify_molecule(atomsymbols::Union{AbstractVector,Tuple})
    identify_molecule(system::AbstractSystem)

Given `atomsymbols`, a list of `Symbol`s corresponding to atoms of a molecule, or `system`
an `AbstractSystem` obeying `AtomsBase` API, attempts to retrieve the name of the molecule.

Examples:
```jldoctest
julia> CrystalEnergyGrids.identify_molecule([:Oa, :C_co2, :Ob])
"CO2"

julia> CrystalEnergyGris.identify_molecule([:Hw, :O5, :Hz])
"H2O"
```
"""
function identify_molecule(atomsymbols::Union{AbstractVector,Tuple})
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
identify_molecule(system) = identify_molecule(atomic_symbol(system))


function get_atom_name(atom)
    name = atom isa String ? atom : String(atom)
    @assert isascii(name)
    any(!isletter, name) || return name
    s = split(name, '_')
    if length(s) > 1
        if all(isnumeric, last(s))
            return join(@view(s[1:end-1]), '_')
        end
    else
        i = length(name)
        while isnumeric(name[i])
            i -= 1
        end
        i < length(name) && return name[1:i]
    end
    return name
end


Base.@assume_effects :foldable function typeof_psystem(::Val{N}) where N
    typeof(ParticleSystem(;
        xpositions=SVector{N,typeof(1.0u"Å")}[],
        ypositions=SVector{N,typeof(1.0u"Å")}[],
        unitcell=SMatrix{3,3,Float64,9}(LinearAlgebra.I)*30.0u"Å",
        cutoff=12.0u"Å",
        parallel=true,
        output=0.0u"K"
    ))
end


function block_average(l; num=nothing, size=nothing)
    n = length(l)
    sizeint = size isa AbstractFloat ? round(Int, n*size) : size
    limits = if num isa Integer
        round.(Int, LinRange(1, n+1, num+1))
    elseif sizeint isa Integer
        _limits = collect(1:sizeint:n+1)
        if last(_limits) != n+1
            push!(_limits, n+1)
        end
        _limits
    else
        error("Please give either num or size as keyword argument")
    end
    [mean(@view(l[limits[i]:(limits[i+1]-1)])) for i in 1:length(limits)-1]
end

# Multithreading

using Base.Threads

struct Counter3D
    A::Vector{Vector{Tuple{Int,Vector{Int}}}}
    lck::ReentrantLock
    counter2::Base.RefValue{Int}
    counter3::Base.RefValue{Int}
end
Counter3D() = Counter3D(Vector{Vector{Int}}[], ReentrantLock(), Ref(0), Ref(0))
@inbounds function Base.getindex(x::Counter3D, i, j, k)
    lock(x.lck)
    n = length(x.A)
    n < i && append!(x.A, Tuple{Int,Vector{Int}}[] for _ in 1:(i-n))
    I = x.A[i]
    m = length(I)
    if m < j
        c2 = x.counter2[]
        newc2 = c2 + j - m
        append!(I, (c, Int[]) for c in c2:(newc2-1))
        x.counter2[] = newc2
    end
    ret2, J = I[j]
    p = length(J)
    if p < k
        c3 = x.counter3[]
        newc3 = c3 + k - p
        append!(J, c3:(newc3-1))
        x.counter3[] = newc3
    end
    ret3 = J[k]
    unlock(x.lck)
    ret2, ret3
end

function stripspawn(@nospecialize(expr), dict::IdDict{Symbol,Symbol}, inspawncall=false)
    if expr isa Expr
        if expr.head === :(=) && Meta.isexpr(expr.args[2], :macrocall) && expr.args[2].args[1] === Symbol("@spawn")
            left = expr.args[1]
            newsymb = if left isa Symbol
                get!(() -> Symbol(left::Symbol, :_nospawn), dict, left)
            else
                stripspawn(left, dict, false)
            end
            Expr(:(=), newsymb, stripspawn(last(expr.args[2].args), dict, true))
        elseif expr.head === :$ && inspawncall
            stripspawn(expr.args[1], dict, inspawncall)
        else
            x = Expr(expr.head)
            append!(x.args, stripspawn(arg, dict, inspawncall) for arg in expr.args)
            x
        end
    elseif expr isa Symbol
        get(dict, expr, expr)
    else
        expr
    end
end

macro spawnif(dospawn, expr)
    stripped = stripspawn(expr, IdDict{Symbol,Symbol}())
    esc(quote
        if $dospawn
            $expr
        else
            $stripped
        end
    end)
end


mutable struct SubAtomic{T}
    @atomic x::T
end

"""
    LoadBalancer{T}

Channel-like abstraction to balance function calls across a fixed number of tasks.
See [`LoadBalancer{T}(f, n::Integer)`](@ref) to create an instance.
"""
struct LoadBalancer{T}
    channel::Channel{T}
    tasks::Vector{Task}
    busy::SubAtomic{Int}
    event::Event
end

using Serialization

"""
    LoadBalancer{T}(f, n::Integer=nthreads())

Given a function `f`, create `n` tasks that wait on an input `x::T` to execute `f(x, i)`
where `i` is between 1 and `n` and uniquely identifies each task.

With `lb = LoadBalancer(f, n)`, use `put!(lb, x)` to send `x` to one of the tasks: this
call is not blocking and will always return immediately, but the `f(x, i)` call will not
occur until one of the `n` tasks is free.

Use `wait(lb)` to wait until all inputs have been processed.

## Example

```julia
lb = LoadBalancer{Tuple{Int,String}}() do (k,s), _
    some_complicated_function(k, s)
end

for i in 1:1000
    put!(lb, (i, "")) # setup works
end

do_something_else()

wait(lb)
```

See also [`@loadbalance`](@ref)
"""
function LoadBalancer{T}(f, n::Integer=nthreads()) where T
    busy::SubAtomic{Int} = SubAtomic{Int}(0)
    event::Event = Event(true)
    channel::Channel{T} = Channel{T}(n)
    tasks = [(@spawn begin
        _busy = $busy
        while true
            x = take!($channel)
            $f(x, $i)
            nbusy = @atomic :acquire_release _busy.x -= 1
            nbusy == 0 && notify($event)
        end
    end) for i in 1:n]
    foreach(errormonitor, tasks)
    LoadBalancer{T}(channel, tasks, busy, event)
end
function Base.put!(lb::LoadBalancer, x)
    @atomic :acquire_release lb.busy.x += 1
    put!(lb.channel, x)
end
function Base.wait(lb::LoadBalancer)
    while (@atomic :acquire lb.busy.x) > 0
        wait(lb.event)
    end
end
function Base.close(lb::LoadBalancer)
    if (@atomic :acquire lb.busy.x) != 0
        @error "Forcibly closing the LoadBalancer may yield errors from the still-working tasks"
    end
    close(lb.channel)
    val = -max(0, @atomic :monotonic lb.busy.x)
    @atomic lb.busy.x = val # make lb.busy zero or negative to signal the close
    notify(lb.event)
end

function loadbalance_body(expr::Expr; balance=0.5)
    Meta.isexpr(expr, :for) || error("@loadbalance requires a `for` loop expression")
    @assert Meta.isexpr(expr.args[1], :(=)) # TODO: handle multi-loops
    @assert length(expr.args[1].args) == 2
    lhs, rhs = expr.args[1].args
    body = expr.args[2]
    quote
        balance = $(esc(balance))
        0 ≤ balance ≤ 1 || error("The number between @loadbalance and the `for` loop must be between 0 and 1.")
        balance /= Base.Threads.nthreads()
        rhs = $(esc(rhs))
        itsize = Base.IteratorSize(rhs)
        len = (itsize isa Base.HasLength || itsize isa Base.HasShape) ? length(rhs) : nothing
        el = Base.IteratorEltype(rhs) isa Base.HasEltype ? eltype(rhs) : Any
        lb = if rhs isa AbstractArray
            i = firstindex(rhs)
            stop = lastindex(rhs)
            lb1 = LoadBalancer{UnitRange{Int}}() do rnge, $(esc(:taskid))
                for $(esc(lhs)) in view(rhs, rnge)
                    $(esc(body))
                end
            end
            count = 0
            while i ≤ stop
                count += 1
                j = i + round(Int, balance*(stop - i))
                put!(lb1.channel, i:j)
                i = j + 1
            end
            i == stop && put!(lb1, i:i)
            @atomic :acquire_release lb1.busy.x += count
            lb1
        else
            lb0 = LoadBalancer{el}() do $(esc(lhs)), $(esc(:taskid))
                $(esc(body))
            end
            # The following sequence is equivalent to doing `put!(lb0, x)` for all x in rhs,
            # but it avoids doing as many @atomic operations
            count = 0
            for x in rhs
                count += 1
                put!(lb0.channel, x)
            end
            @atomic :acquire_release lb0.busy.x += count
            lb0
        end
        wait(lb)
    end
end

"""
    @loadbalance [0.5] for ...

Similar to `Base.Threads.@threads` but using [`LoadBalancer`](@ref) underneath. This allows
using the special variable `taskid` inside the body of the `for` loop, where `taskid` is a
number between 1 and `Base.Threads.nthreads()` that uniquely identifies the current thread.
As a consequence, `taskid` may be used as an index in an array where each thread can write
concurrently without risking race condition.

Scheduling can be tweaked by putting a number between 0 and 1 before the `for` loop. A
value of 0 means that each element of the loop will be sent individually to a thread, while
higher values mean that groups of contiguous elements will be sent together.
This scheduling is only possible when the `for` loop iterates over an `AbstractArray`, such
as a range or an array for example.
"""
macro loadbalance(balance, expr)
    loadbalance_body(expr; balance)
end
macro loadbalance(expr)
    loadbalance_body(expr)
end
