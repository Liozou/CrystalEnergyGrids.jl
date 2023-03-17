using LinearAlgebra: norm, det
using StaticArrays
using OffsetArrays

const ELECTRONIC_CHARGE_UNIT = 1.60217733e-19
const ELECTRIC_CONSTANT = 8.8541878176e-12
const DielectricConstantOfTheMedium = 1.0
const COULOMBIC_CONVERSION_FACTOR = ELECTRONIC_CHARGE_UNIT^2/(4π*ELECTRIC_CONSTANT*ANGSTROM*ENERGY_CONVERSION_FACTOR*DielectricConstantOfTheMedium)

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

function initialize_ewald(_mat::AbstractMatrix, numsites, precision=1e-6)
    mat = SMatrix{3,3,Float64,9}(_mat)
    cutoff_coulomb = 12.0
    ε = cutoff_coulomb*min(0.5, abs(precision))
    tol = sqrt(abs(log(ε)))
    α = sqrt(abs(log(ε*tol)))/cutoff_coulomb
    tol1 = sqrt(-log(ε*4.0*(tol*α)^2))

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

    il_ax, il_bx, il_cx, il_ay, il_by, il_cy, il_az, il_bz, il_cz = inv(mat)

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

    UIon = 0.0
    Eikr = Vector{Complex{Float64}}(undef, numsites)
    Eikr_xy = Vector{Complex{Float64}}(undef, numsites)
    Eikx = OffsetArray(Vector{Complex{Float64}}(undef, (2kx+1)*numsites), -kx*numsites:(kx+1)*numsites-1)
    Eiky = OffsetArray(Vector{Complex{Float64}}(undef, (2ky+1)*numsites), -ky*numsites:(ky+1)*numsites-1)
    Eikz = OffsetArray(Vector{Complex{Float64}}(undef, (2kz+1)*numsites), -kz*numsites:(kz+1)*numsites-1)

    Eikx[-numsites] = Eiky[-numsites] = Eikz[-numsites] = 1.0
    Eikx[0] = Eiky[0] = Eikz[0] = 1.0
    Eikx[numsites] = Eiky[numsites] = Eikz[numsites] = 1.0

    for j in 2:kx
        Eikx[j*numsites] = 1.0
    end
    for j in 2:ky
        Eiky[-j*numsites] = Eiky[j*numsites] = 1.0
    end
    for j in 2:kz
        Eikz[-j*numsites] = Eikz[j*numsites] = 1.0
    end

    idx_c = 0
    for i in 0:kx, j in -ky:ky
        idx_i = i*numsites
        idx_j = j*numsites
        Eikr_xy[1] = Eikx[idx_i]*Eiky[idx_j]
        _r2_c = i^2 + j^2
        for k in -kz:kz
            r2_c = _r2_c + k^2
            ((r2_c >= recip_cutoff2) | (r2_c == 0)) && continue
            idx_k = k*numsites
            idx_c += 1
            UIon += kfactors[idx_c]*abs2(Eikr_xy[1]*Eikz[idx_k])
        end
    end

    UIon = COULOMBIC_CONVERSION_FACTOR*α/sqrt(π) - UIon



end




