using LinearAlgebra: norm, det
using StaticArrays
using OffsetArrays
using AtomsBase
using SpecialFunctions: erf

export initialize_ewald, compute_ewald

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

complexofpoint(x) = begin (s,c) = sincos(x); ComplexF64(c, s) end


function setup_Eik(systems, numsites, kx, ky, kz, invmat, (superA, superB, superC))
    Eikx = Vector{ComplexF64}(undef, (kx+1)*numsites)
    Eiky = OffsetArray(Vector{ComplexF64}(undef, (2ky+1)*numsites), 1-ky*numsites:(ky+1)*numsites)
    Eikz = OffsetArray(Vector{ComplexF64}(undef, (2kz+1)*numsites), 1-kz*numsites:(kz+1)*numsites)

    i = 0
    for syst in systems
        for sA in 0:(superA-1), sB in 0:(superB-1), sC in 0:(superC-1)
            for j in 1:length(syst)
                i += 1
                Eikx[i] = 1.0
                Eiky[i] = 1.0
                Eikz[i] = 1.0
                i₊ = numsites+i
                i₋ = -numsites+i
                px, py, pz = invmat * (position(syst, j) ./ ANG_UNIT)
                px = 2π*(px + sA/superA)
                py = 2π*(py + sB/superB)
                pz = 2π*(pz + sC/superC)
                Eikx[i₊] = complexofpoint(px)
                tmp = Eiky[i₊] = complexofpoint(py)
                Eiky[i₋] = conj(tmp)
                tmp = Eikz[i₊] = complexofpoint(pz)
                Eikz[i₋] = conj(tmp)
            end
        end
    end

    for j in 2:kx, i in 1:numsites
        Eikx[j*numsites+i] = Eikx[(j-1)*numsites+i]*Eikx[numsites+i]
    end
    for j in 2:ky, i in 1:numsites
        tmp = Eiky[j*numsites+i] = Eiky[(j-1)*numsites+i]*Eiky[numsites+i]
        Eiky[-j*numsites+i] = conj(tmp)
    end
    for j in 2:kz, i in 1:numsites
        tmp = Eikz[j*numsites+i] = Eikz[(j-1)*numsites+i]*Eikz[numsites+i]
        Eikz[-j*numsites+i] = conj(tmp)
    end

    return Eikx, Eiky, Eikz
end

function ewald_main_loop(allcharges, numsites, kx, ky, kz, recip_cutoff2, num_kvecs, Eikx, Eiky, Eikz, Π)
    idx = 0
    sums = zeros(ComplexF64, num_kvecs)
    Eikr_xy = Vector{ComplexF64}(undef, numsites)
    @inbounds for jx in 0:kx, jy in -ky:ky
        for i in 1:numsites
            Eikr_xy[i] = Eikx[jx*numsites+i]*Eiky[jy*numsites+i]
        end
        for jz in -kz:kz
            r2 = jx^2 + jy^2 + jz^2
            ((r2 == 0) | (r2 >= recip_cutoff2)) && continue
            idx += 1
            ofs = jz*numsites
            i = 0
            for charges in allcharges
                m = length(charges)
                for _ in 1:Π
                    for c in 1:m
                        Eikr = Eikr_xy[i+c]*Eikz[ofs+i+c]
                        sums[idx] += charges[c] * Eikr
                    end
                    i += m
                end
            end
        end
    end
    sums
end

"""
    initialize_ewald(syst::AbstractSystem{3}, precision=1e-6)

Given `syst`, which contains a fixed system, return an object `x` to feed to
`compute_ewald(x, mol)` to compute the Fourier contribution in the Ewald summation for the
interaction between fixed system `syst` and molecule `mol`.
"""
function initialize_ewald(syst::AbstractSystem{3}, supercell=(1,1,1), precision=1e-6)
    @assert all(==(Periodic()), boundary_conditions(syst))

    cutoff_coulomb = 12.0
    ε = cutoff_coulomb*min(0.5, abs(precision))
    tol = sqrt(abs(log(ε)))
    α = sqrt(abs(log(ε*tol)))/cutoff_coulomb
    tol1 = sqrt(-log(ε*4.0*(tol*α)^2))

    mat = SMatrix{3,3,Float64,9}((stack(bounding_box(syst) .* supercell) ./ ANG_UNIT))
    len_a, len_b, len_c = cell_lengths(mat)
    volume = abs(det(mat))
    __α = α*tol1/π
    kx = nint(0.25 + __α*len_a)
    ky = nint(0.25 + __α*len_b)
    kz = nint(0.25 + __α*len_c)

    recip_cutoff2 = (1.05*max(kx, ky, kz))^2
    num_kvecs = 0
    for i in 0:kx, j in -ky:ky, k in -kz:kz
        r2_a = i^2 + j^2 + k^2
        num_kvecs += (r2_a != 0) & (r2_a < recip_cutoff2)
    end

    invmat::SMatrix{3,3,Float64,9} = inv(mat)
    il_ax, il_ay, il_az, il_bx, il_by, il_bz, il_cx, il_cy, il_cz = invmat
    volume_factor = COULOMBIC_CONVERSION_FACTOR*π/volume
    α_factor = -0.25/α^2

    # kvecs = Vector{SVector{3,Float64}}(undef, num_kvecs)
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
                # kvecs[idx_b] = SVector{3,Float64}(rkx, rky, rkz)
                rksqr = rkx^2 + rky^2 + rkz^2
                kfactors[idx_b] = volume_factor*(2.0*(1+(i!=0)))*exp(α_factor*rksqr)/rksqr
            end
        end
    end

    Π = prod(supercell)
    numsites = length(syst)*Π

    UIon = COULOMBIC_CONVERSION_FACTOR*α/sqrt(π) - sum(kfactors)
    charges::Vector{Float64} = syst[:,:atomic_charge] ./ CHARGE_UNIT
    energy_framework_self = sum(abs2, charges)*(COULOMBIC_CONVERSION_FACTOR/sqrt(π))*α

    Eikx, Eiky, Eikz = setup_Eik((syst,), numsites, kx, ky, kz, invmat, supercell)
    StoreRigidChargeFramework = ewald_main_loop((charges,), numsites, kx, ky, kz, recip_cutoff2, num_kvecs, Eikx, Eiky, Eikz, Π)
    # UChargeChargeFrameworkRigid = 0.0
    # for idx in 1:num_kvecs
    #     UChargeChargeFrameworkRigid += kfactors[idx]*abs2(StoreRigidChargeFramework[idx])
    # end
    # UChargeChargeFrameworkRigid += UIon*sum(charges)^2
    net_charges_framework = sum(charges)

    return α, kx, ky, kz, invmat, recip_cutoff2, num_kvecs, volume_factor, energy_framework_self, kfactors, UIon, StoreRigidChargeFramework, net_charges_framework
end



function compute_ewald((α, kx, ky, kz, invmat, recip_cutoff2, num_kvecs, volume_factor, energy_framework_self, kfactors, UIon, StoreRigidChargeFramework, net_charges_framework),
                       systems)
    numsites = sum(length, systems)
    allcharges::Vector{Vector{Float64}} = [syst[:,:atomic_charge] ./ CHARGE_UNIT for syst in systems]
    chargefactor = (COULOMBIC_CONVERSION_FACTOR/sqrt(π))*α
    energy_adsorbate_self = sum(sum(abs2, charges)*chargefactor for charges in allcharges)
    net_charges = sum(sum(charges) for charges in allcharges)

    Eikx, Eiky, Eikz = setup_Eik(systems, numsites, kx, ky, kz, invmat, (1,1,1))
    StoreTotalChargeAdsorbates = ewald_main_loop(allcharges, numsites, kx, ky, kz, recip_cutoff2, num_kvecs, Eikx, Eiky, Eikz, 1)
    energy_framework_adsorbate = 0.0
    energy_adsorbate_adsorbate = 0.0
    for idx in 1:num_kvecs
        temp = kfactors[idx]
        _re_f, _im_f = reim(StoreRigidChargeFramework[idx])
        _re_a, _im_a = reim(StoreTotalChargeAdsorbates[idx])
        energy_framework_adsorbate += temp*(_re_f*_re_a + _im_f*_im_a)
        energy_adsorbate_adsorbate += temp*abs2(StoreTotalChargeAdsorbates[idx])
    end

    UHostAdsorbateChargeChargeFourier = 2*(energy_framework_adsorbate + UIon*net_charges_framework*net_charges)

    m = length(systems)
    energy_adsorbate_excluded = 0.0
    for i in 1:m
        syst = systems[i]
        charges = allcharges[i]
        n = length(syst)
        for atomA in 1:n
            chargeA = charges[atomA]
            posA = position(syst, atomA) ./ ANG_UNIT
            for atomB in (atomA+1):n
                chargeB = charges[atomB]
                posB = position(syst, atomB) ./ ANG_UNIT
                r = norm(posB .- posA)
                energy_adsorbate_excluded += erf(α*r)*chargeA*chargeB/r
            end
        end
    end
    energy_adsorbate_excluded *= COULOMBIC_CONVERSION_FACTOR

    UAdsorbateAdsorbateChargeChargeFourier = energy_adsorbate_adsorbate - energy_adsorbate_self - energy_adsorbate_excluded + UIon*net_charges^2

    UHostAdsorbateChargeChargeFourier*ENERGY_TO_KELVIN, UAdsorbateAdsorbateChargeChargeFourier*ENERGY_TO_KELVIN
end
