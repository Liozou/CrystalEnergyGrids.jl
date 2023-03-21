using LinearAlgebra: norm, det
using StaticArrays
using OffsetArrays
using AtomsBase

export initialize_ewald

function cell_lengths(mat::AbstractMatrix)
    a, b, c = eachcol(mat)
    return norm(a), norm(b), norm(c)
end
function cell_angles(mat::AbstractMatrix)
    _a, _b, _c = eachcol(mat)
    a, b, c = cell_lengths(mat)
    α = acosd(_b'_c/(b*c))
    β = acosd(_c'_a/(c*a))
    γ = acosd(_a'_b/(a*b))
    return (a, b, c), (α, β, γ)
end

nint(x) = ifelse(x>=0.0, x+0.5, x-0.5)
complexofpoint(x) = begin (s,c) = sincos(x); Complex{Float64}(c, s) end

"""
    initialize_ewald(syst::AbstractSystem{3}, precision=1e-6)

Given `syst`, which contains a fixed system, return an object `x` to feed to
`compute_ewald(x, mol)` to compute the Fourier contribution in the Ewald summation for the
interaction between fixed system `syst` and molecule `mol`.
"""
function initialize_ewald(syst::AbstractSystem{3}, precision=1e-6)
    @assert all(==(Periodic()), boundary_conditions(syst))

    cutoff_coulomb = 12.0
    ε = cutoff_coulomb*min(0.5, abs(precision))
    tol = sqrt(abs(log(ε)))
    α = sqrt(abs(log(ε*tol)))/cutoff_coulomb
    tol1 = sqrt(-log(ε*4.0*(tol*α)^2))

    mat = SMatrix{3,3,Float64,9}((stack(bounding_box(syst)) ./ ANG_UNIT))
    len_a, len_b, len_c = cell_lengths(mat)
    volume = abs(det(mat))
    __α = α*tol1/π
    kx = floor(Int, nint(0.25 + __α*len_a))
    ky = floor(Int, nint(0.25 + __α*len_b))
    kz = floor(Int, nint(0.25 + __α*len_c))

    recip_cutoff2 = (1.05*max(kx, ky, kz))^2
    num_kvecs = 0
    for i in 0:kx, j in -ky:ky, k in -kz:kz
        r2_a = i^2 + j^2 + k^2
        num_kvecs += (r2_a != 0) & (r2_a < recip_cutoff2)
    end

    invmat::SMatrix{3,3,Float64,9} = inv(mat)
    il_ax, il_bx, il_cx, il_ay, il_by, il_cy, il_az, il_bz, il_cz = invmat

    volume_factor = COULOMBIC_CONVERSION_FACTOR*π/volume
    α_factor = -0.25/α^2

    kvecs = Vector{SVector{3,Float64}}(undef, num_kvecs)
    kfactors = Vector{Float64}(undef, num_kvecs)

    idx_b = 0
    for i in 0:kx
        rk2x = i*il_ax
        rk2y = i*il_bx
        rk2z = i*il_cx
        for j in -ky:ky
            rk1x = rk2x + j*il_ay
            rk1y = rk2y + j*il_by
            rk1z = rk2z + j*il_cy
            for k in -kz:kz
                r2_b = i^2 + j^2 + k^2
                ((r2_b >= recip_cutoff2) | (r2_b == 0)) && continue
                idx_b += 1
                rkx = 2π*(rk1x + k*il_az)
                rky = 2π*(rk1y + k*il_bz)
                rkz = 2π*(rk1z + k*il_cz)
                kvecs[idx_b] = SVector{3,Float64}(rkx, rky, rkz)
                rksqr = rkx^2 + rky^2 + rkz^2
                kfactors[idx_b] = volume_factor*(2.0+2.0*(i!=0))*exp(α_factor*rksqr)/rksqr
            end
        end
    end

    poss = position(syst)
    numsites = length(poss)

    UIon = COULOMBIC_CONVERSION_FACTOR*α/sqrt(π) - sum(kfactors)

    Eikx = OffsetArray(Vector{Complex{Float64}}(undef, (kx+1)*numsites), 0:(kx+1)*numsites-1)
    Eiky = OffsetArray(Vector{Complex{Float64}}(undef, (2ky+1)*numsites), -ky*numsites:(ky+1)*numsites-1)
    Eikz = OffsetArray(Vector{Complex{Float64}}(undef, (2kz+1)*numsites), -kz*numsites:(kz+1)*numsites-1)

    for i in 0:(numsites-1)
        Eikx[i] = 1.0
        Eiky[i] = 1.0
        Eikz[i] = 1.0
        i₊ = numsites+i
        i₋ = -numsites+i
        px, py, pz = invmat * (poss[i+1] ./ (ANG_UNIT/(2π)))
        Eikx[i₊] = complexofpoint(px)
        tmp = Eiky[i₊] = complexofpoint(py)
        Eiky[i₋] = conj(tmp)
        tmp = Eikz[i₊] = complexofpoint(pz)
        Eikz[i₋] = conj(tmp)
    end

    for j in 2:kx, i in 0:(numsites-1)
        Eikx[j*numsites+i] = Eikx[(j-1)*numsites+i]*Eikx[numsites+i]
    end
    for j in 2:ky, i in 0:(numsites-1)
        tmp = Eiky[j*numsites+i] = Eiky[(j-1)*numsites+i]*Eiky[numsites+i]
        Eiky[-j*numsites+i] = conj(tmp)
    end
    for j in 2:kz, i in 0:(numsites-1)
        tmp = Eikz[j*numsites+i] = Eikz[(j-1)*numsites+i]*Eikz[numsites+i]
        Eikz[-j*numsites+i] = conj(tmp)
    end

    charges::Vector{Float64} = syst[:,:atomic_charge] ./ CHARGE_UNIT

    Eikr_xy = Vector{Complex{Float64}}(undef, numsites)
    Eikr = Vector{Complex{Float64}}(undef, numsites)
    StoreRigidChargeFramework = zeros(Complex{Float64}, num_kvecs)
    UChargeChargeFrameworkRigid = 0.0
    idx = 0
    for jx in 0:kx, jy in -ky:ky
        for i in 0:(numsites-1)
            Eikr_xy[i+1] = Eikx[jx*numsites+i]*Eiky[jy*numsites+i]
        end
        for jz in -kz:kz
            r2 = jx^2 + jy^2 + jz^2
            ((r2 == 0) | (r2 >= recip_cutoff2)) && continue
            idx += 1
            # rk = kvecs[idx] # only useful with dipoles?
            for i in 1:numsites
                Eikr[i] = Eikr_xy[i]*Eikz[jz*numsites+i-1]
                StoreRigidChargeFramework[idx] += charges[i] * Eikr[i]
            end
            UChargeChargeFrameworkRigid += kfactors[idx]*abs2(StoreRigidChargeFramework[idx])
        end
    end

    UChargeChargeFrameworkRigid += UIon*sum(charges)^2

    return α, kx, ky, kz, recip_cutoff2, volume_factor, Eikx, Eiky, Eikz, kvecs, kfactors, UIon, StoreRigidChargeFramework, UChargeChargeFrameworkRigid
end





