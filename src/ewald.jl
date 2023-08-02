using LinearAlgebra: norm, det
using StaticArrays
using OffsetArrays
using AtomsBase
using SpecialFunctions: erf, erfc

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


struct EwaldKspace
    ks::NTuple{3,Int}
    num_kvecs::Int
    kindices::Vector{Tuple{Int,Int,UnitRange{Int},Int}}
end

struct EwaldFramework
    kspace::EwaldKspace
    α::Float64
    invmat::SMatrix{3,3,Float64,9}
    volume_factor::Float64
    energy_framework_self::Float64
    kfactors::Vector{Float64}
    UIon::Float64
    StoreRigidChargeFramework::Vector{ComplexF64}
    net_charges_framework::Float64
    ε::Float64
end


function make_line_pos!(Eikt, start, numk, eikt)
    Eikt[start] = 1
    Eikt[start+1] = eikt
    for k in 2:numk
        Eikt[start+k] = Eikt[start+k-1]*eikt
    end
    nothing
end

function make_line_neg!(Eikt, start, numk, eikt)
    ceikt = conj(eikt)
    mid = start + numk
    Eikt[mid-1] = ceikt
    for k in 2:numk
        Eikt[mid-k] = Eikt[mid-k+1]*ceikt
    end
    make_line_pos!(Eikt, mid, numk, eikt)
end


"""
    setup_Eik!((Eikx, Eiky, Eikz), (kx, ky, kz), systems, invmat, (ΠA, ΠB, ΠC))

Setup the computation of the e^(ik...) necessary for Ewald summation, for the given
`systems`, with `numsites` total number of point charges (including those in the supercell)
along a `kx × ky × kz` K-space grid.
`invmat` is the inverse of the real-space unit cell. That unit cell is multiplied by
`ΠA × ΠB × ΠC` to obtain the supercell used for the computation.

Return `(Eikx, Eiky, Eikz)` where, for `m` in `1:numsites`
- `Eikx[(m-1)*(kx+1) + 1 + k] == exp(2ikπ * fracpos_x[m])` with `k ∈ 0:kx`
- `Eiky[(m-1)*(2ky+1) + ky + 1 + k] == exp(±2ikπ * fracpos_y[m])` with `k ∈ -ky:ky`
- `Eikz[(m-1)*(2ky+1) + ky + 1 + k] == exp(±2ikπ * fracpos_z[m])` with `k ∈ -kz:kz`

with `fracpos_x[m], fracpos_y[m], fracpos_z[m] == invmat*position[m]`
"""
function setup_Eik(systems, (kx, ky, kz), invmat, (ΠA, ΠB, ΠC))
    ΠBC = ΠB*ΠC
    ΠABC = ΠA*ΠBC
    numsites = ΠABC * sum(length, systems)

    Eikx = Vector{ComplexF64}(undef, (kx+1)*numsites)
    Eiky = Vector{ComplexF64}(undef, (2ky+1)*numsites)
    Eikz = Vector{ComplexF64}(undef, (2kz+1)*numsites)

    ofs = -1
    @inbounds for syst in systems
        for j in 1:length(syst)
            px, py, pz = invmat * (NoUnits.(position(syst, j)/u"Å"))
            for πs in CartesianIndices((0:(ΠA-1), 0:(ΠB-1), 0:(ΠC-1)))
                πa, πb, πc = Tuple(πs)
                πofs = πa*ΠBC + πb*ΠC + πc

                ix = 1 + ((j+ofs)*ΠABC + πofs)*(kx+1)
                eikx = cispi(2*(px + πa/ΠA))
                make_line_pos!(Eikx, ix, kx, eikx)

                eiky = cispi(2*(py + πb/ΠB))
                iy = 1 + ((j+ofs)*ΠABC + πofs)*(2ky+1)
                make_line_neg!(Eiky, iy, ky, eiky)

                eikz = cispi(2*(pz + πc/ΠC))
                iz = 1 + ((j+ofs)*ΠABC + πofs)*(2kz+1)
                make_line_neg!(Eikz, iz, kz, eikz)
            end
        end
        ofs += length(syst)
    end

    (Eikx, Eiky, Eikz)
end


function ewald_main_loop(allcharges, kspace::EwaldKspace, Eiks)
    sums = zeros(ComplexF64, kspace.num_kvecs)

    Eikx, Eiky, Eikz = Eiks
    kx, ky, kz = kspace.ks
    kxp = kx + 1
    tkyp = 2ky + 1
    tkzp = 2kz + 1

    ix = 0 # reference index of the current site for Eikx
    iy = ky+1 # reference index of the current site for Eiky
    iz = kz+1 # reference index of the current site for Eikz

    @inbounds for charges in allcharges, c in charges
        for (jy, jz, jxrange, rangeidx) in kspace.kindices
            Eik_xy = c*Eiky[iy+jy]*Eikz[iz+jz]
            n = length(jxrange)
            ofs = first(jxrange) + ix
            @simd for I in 1:n
                sums[rangeidx+I] += Eikx[ofs+I]*Eik_xy
            end
        end
        ix += kxp
        iy += tkyp
        iz += tkzp
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

    mat = stack3(bounding_box(syst) .* supercell)
    len_a, len_b, len_c = cell_lengths(mat)
    volume = abs(det(mat))
    __α = α*tol1/π
    kx = nint(0.25 + __α*len_a)
    ky = nint(0.25 + __α*len_b)
    kz = nint(0.25 + __α*len_c)

    recip_cutoff2 = (1.05*max(kx, ky, kz))^2
    num_kvecs = 0
    kindices = Tuple{Int,Int,UnitRange{Int},Int}[]
    for j in -ky:ky, k in -kz:kz
        start = -1
        for i in 0:kx
            r2_a = i^2 + j^2 + k^2
            if (r2_a != 0) & (r2_a < recip_cutoff2)
                num_kvecs += 1
                if start == -1
                    start = i
                end
            elseif start != -1
                _, _, lastrange, lastrangeidx = isempty(kindices) ? (0,0,1:0,0) : last(kindices)
                push!(kindices, (j, k, start:(i-1), lastrangeidx + length(lastrange)))
                start = -1
                break
            end
        end
        if start != -1
            _, _, lastrange, lastrangeidx = isempty(kindices) ? (0,0,1:0,0) : last(kindices)
            push!(kindices, (j, k, start:kx, lastrangeidx + length(lastrange)))
        end
    end

    invmat::SMatrix{3,3,Float64,9} = inv(mat)
    il_ax, il_ay, il_az, il_bx, il_by, il_bz, il_cx, il_cy, il_cz = invmat
    volume_factor = COULOMBIC_CONVERSION_FACTOR*π/volume
    α_factor = -0.25/α^2

    # kvecs = Vector{SVector{3,Float64}}(undef, num_kvecs)
    kfactors = Vector{Float64}(undef, num_kvecs)

    @threads for (j, k, irange, rangeidx) in kindices
        rk0x = j*il_ay + k*il_az
        rk0y = j*il_by + k*il_bz
        rk0z = j*il_cy + k*il_cz
        for (I, i) in enumerate(irange)
            rkx = 2π*(rk0x + i*il_ax)
            rky = 2π*(rk0y + i*il_bx)
            rkz = 2π*(rk0z + i*il_cx)
            # kvecs[idx_b] = SVector{3,Float64}(rkx, rky, rkz)
            rksqr = rkx^2 + rky^2 + rkz^2
            kfactors[rangeidx+I] = volume_factor*(2.0*(1+(i!=0)))*exp(α_factor*rksqr)/rksqr
        end
    end

    UIon = COULOMBIC_CONVERSION_FACTOR*α/sqrt(π) - sum(kfactors)
    charges::Vector{Float64} = Float64.(syst[:,:atomic_charge]/u"e_au")
    energy_framework_self = sum(abs2, charges)*(COULOMBIC_CONVERSION_FACTOR/sqrt(π))*α

    kspace = EwaldKspace((kx, ky, kz), num_kvecs, kindices)
    Eiks = setup_Eik((syst,), kspace.ks, invmat, supercell)
    StoreRigidChargeFramework = ewald_main_loop((repeat(charges, prod(supercell)),), kspace, Eiks)
    # UChargeChargeFrameworkRigid = 0.0
    # for idx in 1:num_kvecs
    #     UChargeChargeFrameworkRigid += kfactors[idx]*abs2(StoreRigidChargeFramework[idx])
    # end
    # UChargeChargeFrameworkRigid += UIon*sum(charges)^2
    net_charges_framework = sum(charges)

    return EwaldFramework(kspace, α, invmat, volume_factor, energy_framework_self,
                          kfactors, UIon, StoreRigidChargeFramework, net_charges_framework,
                          ε)
end


struct EwaldContext
    eframework::EwaldFramework
    Eiks::NTuple{3,Vector{ComplexF64}}
    allcharges::Vector{Vector{Float64}}
    static_contribution::Float64
    energy_net_charges::Float64
end

function EwaldContext(eframework::EwaldFramework, systems)
    allcharges::Vector{Vector{Float64}} = [NoUnits.(syst[:,:atomic_charge]/u"e_au") for syst in systems]
    chargefactor = (COULOMBIC_CONVERSION_FACTOR/sqrt(π))*eframework.α
    energy_adsorbate_self = sum(sum(abs2, charges)*chargefactor for charges in allcharges)
    net_charges = sum(sum(charges) for charges in allcharges)
    static_contribution = energy_adsorbate_self +  eframework.UIon*net_charges^2
    energy_net_charges = eframework.UIon*eframework.net_charges_framework*net_charges

    Eiks = setup_Eik(systems, eframework.kspace.ks, eframework.invmat, (1,1,1))
    EwaldContext(eframework, Eiks, allcharges, static_contribution, energy_net_charges)
end


function compute_ewald(ctx::EwaldContext, systems)
    newcharges = ewald_main_loop(ctx.allcharges, ctx.eframework.kspace, ctx.Eiks)
    energy_framework_adsorbate = 0.0
    energy_adsorbate_adsorbate = 0.0
    for idx in 1:ctx.eframework.kspace.num_kvecs
        temp = ctx.eframework.kfactors[idx]
        _re_f, _im_f = reim(ctx.eframework.StoreRigidChargeFramework[idx])
        _re_a, _im_a = reim(newcharges[idx])
        energy_framework_adsorbate += temp*(_re_f*_re_a + _im_f*_im_a)
        energy_adsorbate_adsorbate += temp*abs2(newcharges[idx])
    end

    UHostAdsorbateChargeChargeFourier = 2*(energy_framework_adsorbate + ctx.energy_net_charges)

    m = length(systems)
    energy_adsorbate_excluded = 0.0
    for i in 1:m
        syst = systems[i]
        charges = ctx.allcharges[i]
        n = length(syst)
        for atomA in 1:n
            chargeA = charges[atomA]
            posA = Float64.(position(syst, atomA)/u"Å")
            for atomB in (atomA+1):n
                chargeB = charges[atomB]
                posB = Float64.(position(syst, atomB)/u"Å")
                r = norm(posB .- posA)
                energy_adsorbate_excluded += erf(ctx.eframework.α*r)*chargeA*chargeB/r
            end
        end
    end
    energy_adsorbate_excluded *= COULOMBIC_CONVERSION_FACTOR

    UAdsorbateAdsorbateChargeChargeFourier = energy_adsorbate_adsorbate - energy_adsorbate_excluded - ctx.static_contribution

    UHostAdsorbateChargeChargeFourier*ENERGY_TO_KELVIN, UAdsorbateAdsorbateChargeChargeFourier*ENERGY_TO_KELVIN
end

function compute_ewald(eframework::EwaldFramework, systems)
    compute_ewald(EwaldContext(eframework, systems), systems)
end

function derivatives_ewald(ewald::EwaldFramework, charge::Float64, r2::Float64)
    r = sqrt(r2)
    r3 = r2*r
    r5 = r3*r2
    α = ewald.α
    r2α2 = r2*α^2
    er2α2 = 2α*r*exp(-r2α2)/sqrt(π)
    erfαr = erfc(α*r)
    value = charge*erfαr/r
    ∂1 = -charge*(er2α2 + erfαr)/r3
    ∂2 = charge*(er2α2*(3+2*r2α2) + 3*erfαr)/r5
    ∂3 = charge*(-er2α2*(15 + 10*r2α2 + 4*r2α2^2) - 15*erfαr)/(r5*r2)
    (value, ∂1, ∂2, ∂3)
end
