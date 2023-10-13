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

"""
    EwaldKspace

Record of the k-vectors used for a computation, characterized by the unit cell, the coulomb
cutoff (12.0 Å in this package) and the required precision (default to 1e-6).
"""
struct EwaldKspace
    ks::NTuple{3,Int}
    num_kvecs::Int
    kindices::Vector{Tuple{Int,Int,UnitRange{Int},Int}}
end
Base.:(==)(e1::EwaldKspace, e2::EwaldKspace) = e1.ks == e2.ks && e1.num_kvecs == e2.num_kvecs && e1.kindices == e2.kindices

"""
    EwaldFramework

Pre-computation record for the Ewald summation scheme on a fixed framework.
"""
struct EwaldFramework
    kspace::EwaldKspace
    α::Float64
    mat::SMatrix{3,3,Float64,9}
    invmat::SMatrix{3,3,Float64,9}
    kfactors::Vector{Float64}
    UIon::Float64
    StoreRigidChargeFramework::Vector{ComplexF64}
    net_charges_framework::Float64
    precision::Float64
end
function EwaldFramework(mat::SMatrix{3,3,Float64,9})
    EwaldFramework(EwaldKspace((0,0,0), 0, Tuple{Int,Int,UnitRange{Int},Int}[]),
                   0.0, mat, inv(mat), Float64[], 0.0, ComplexF64[], 0.0, 0.0)
end

function Base.:(==)(e1::EwaldFramework, e2::EwaldFramework)
    e1.kspace == e2.kspace && e1.α == e2.α &&
                              e1.mat == e2.mat &&
                              e1.invmat == e2.invmat &&
                              e1.kfactors == e2.kfactors &&
                              e1.UIon == e2.UIon &&
                              e1.StoreRigidChargeFramework == e2.StoreRigidChargeFramework &&
                              e1.net_charges_framework == e2.net_charges_framework &&
                              e1.precision == e2.precision
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
    numsites = ΠABC * sum(length, systems; init=0)

    kxp = kx + 1
    tkyp = ky + ky + 1
    tkzp = kz + kz + 1

    Eikx = Vector{ComplexF64}(undef, (kx+1)*numsites)
    Eiky = Vector{ComplexF64}(undef, (ky+ky+1)*numsites)
    Eikz = Vector{ComplexF64}(undef, (kz+kz+1)*numsites)

    ofs = -1
    @inbounds for syst in systems
        for j in 1:length(syst)
            px, py, pz = invmat * (NoUnits.(position(syst, j)/u"Å"))
            for πs in CartesianIndices((0:(ΠC-1), 0:(ΠB-1), 0:(ΠA-1)))
                πc, πb, πa = Tuple(πs)
                jofs = (j + ofs)*ΠABC + πa*ΠBC + πb*ΠC + πc

                ix = 1 + jofs*kxp
                eikx = cispi(2*(px + πa/ΠA))
                make_line_pos!(Eikx, ix, kx, eikx)

                eiky = cispi(2*(py + πb/ΠB))
                iy = 1 + jofs*tkyp
                make_line_neg!(Eiky, iy, ky, eiky)

                eikz = cispi(2*(pz + πc/ΠC))
                iz = 1 + jofs*tkzp
                make_line_neg!(Eikz, iz, kz, eikz)
            end
        end
        ofs += length(syst)
    end

    (Eikx, Eiky, Eikz)
end


function ewald_main_loop!(sums::T, allcharges, kspace::EwaldKspace, Eiks, skipcontribution) where T
    fill!(sums, zero(ComplexF64))
    Eikx, Eiky, Eikz = Eiks
    kx, ky, kz = kspace.ks

    kxp = kx + 1
    tkyp = ky + ky + 1
    tkzp = kz + kz + 1

    ix = -kxp # reference index of the current site for Eikx
    iy = -ky # reference index of the current site for Eiky
    iz = -kz # reference index of the current site for Eikz

    counter = 0

    @inbounds for (idxsystem, charges) in enumerate(allcharges)
        counter += 1
        counter == skipcontribution && continue
        for c in charges
            ix += kxp
            iy += tkyp
            iz += tkzp
            # if !iszero(c) # TODO: remove zero charges when constructing EwaldContext
                for (jy, jz, jxrange, rangeidx) in kspace.kindices
                    Eik_yz = c*Eiky[iy+jy]*Eikz[iz+jz]
                    n = length(jxrange)
                    ofs = first(jxrange) + ix
                    if T === Vector{ComplexF64}
                        @simd for I in 1:n
                            sums[rangeidx+I] += Eikx[ofs+I]*Eik_yz
                        end
                    else
                        @simd for I in 1:n
                            sums[rangeidx+I,idxsystem] += Eikx[ofs+I]*Eik_yz
                        end
                    end
                end
            # end
        end
    end
    if T !== Vector{ComplexF64}
        sums[:,end] .= sum(sums; dims=2)
    end
    nothing
end


"""
    initialize_ewald(syst::AbstractSystem{3}, supercell=find_supercell(syst, 12.0u"Å"), precision=1e-6)

Given `syst`, which contains a fixed system, return an `EwaldFramework` object `x` to feed
to [`compute_ewald(x, molecules)`](@ref compute_ewald) to compute the Fourier contribution
in the Ewald summation for the interaction between fixed system `syst` and `molecules`.
"""
function initialize_ewald(syst::AbstractSystem{3}, supercell=find_supercell(syst, 12.0u"Å"), precision=1e-6)
    @assert all(==(Periodic()), boundary_conditions(syst))

    cutoff_coulomb = 12.0
    ε = cutoff_coulomb*min(0.5, abs(precision))
    tol = sqrt(abs(log(ε)))
    α = sqrt(abs(log(ε*tol)))/cutoff_coulomb
    tol1 = sqrt(-log(ε*4.0*(tol*α)^2))

    mat = stack3(bounding_box(syst) .* supercell)
    len_a, len_b, len_c = cell_lengths(mat)
    __α = α*tol1/π
    kx = nint(0.25 + __α*len_a)
    ky = nint(0.25 + __α*len_b)
    kz = nint(0.25 + __α*len_c)

    recip_cutoff2 = (1.05*max(kx, ky, kz))^2
    num_kvecs = 0
    kindices = Tuple{Int,Int,UnitRange{Int},Int}[]
    for j in -ky:ky, k in -kz:kz
        started = false
        for i in 0:kx
            r2_a = i^2 + j^2 + k^2
            if (r2_a != 0) & (r2_a < recip_cutoff2)
                num_kvecs += 1
                if !started
                    started = true
                end
            elseif started
                _, _, lastrange, lastrangeidx = isempty(kindices) ? (0,0,1:0,0) : last(kindices)
                push!(kindices, (j, k, (j==k==0):(i-1), lastrangeidx + length(lastrange)))
                started = false
                break
            end
        end
        if started
            _, _, lastrange, lastrangeidx = isempty(kindices) ? (0,0,1:0,0) : last(kindices)
            push!(kindices, (j, k, (j==k==0):kx, lastrangeidx + length(lastrange)))
        end
    end

    invmat::SMatrix{3,3,Float64,9} = inv(mat)
    il_ax, il_ay, il_az, il_bx, il_by, il_bz, il_cx, il_cy, il_cz = invmat
    volume_factor = COULOMBIC_CONVERSION_FACTOR*2π/(det(mat))
    @assert volume_factor > 0
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
            kfactors[rangeidx+I] = volume_factor*(1+(i!=0))*exp(α_factor*rksqr)/rksqr
        end
    end

    UIon = COULOMBIC_CONVERSION_FACTOR*α/sqrt(π) - sum(kfactors; init=0.0)
    charges::Vector{Float64} = Float64.(syst[:,:atomic_charge]/u"e_au")

    kspace = EwaldKspace((kx, ky, kz), num_kvecs, kindices)
    Eiks = setup_Eik((syst,), kspace.ks, invmat, supercell)
    StoreRigidChargeFramework = Vector{ComplexF64}(undef, num_kvecs)
    Π = prod(supercell)
    repeatedcharges = [c for c in charges for _ in 1:Π]
    ewald_main_loop!(StoreRigidChargeFramework, repeatedcharges, kspace, Eiks, 0)
    # UChargeChargeFrameworkRigid = 0.0
    # for idx in 1:num_kvecs
    #     UChargeChargeFrameworkRigid += kfactors[idx]*abs2(StoreRigidChargeFramework[idx])
    # end
    net_charges_framework = sum(charges; init=0.0)*Π
    # UChargeChargeFrameworkRigid += UIon*net_charges_framework^2

    return EwaldFramework(kspace, α, mat, invmat, kfactors, UIon,
                          StoreRigidChargeFramework, net_charges_framework, precision)
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


"""
    EwaldContext

Context to be fed to [`compute_ewald`](@ref) to compute the Fourier contribution to the
energy of a system composed of a rigid framework and any number of rigid molecules.
"""
struct EwaldContext
    eframework::EwaldFramework
    Eiks::NTuple{3,Vector{ComplexF64}}
    allcharges::Vector{Vector{Float64}}
    offsets::Vector{Int}
    # offsets[i] is the number of charges before the start of the system at allcharges[i] minus 1
    static_contribution::Float64
    energy_net_charges::Float64
end

function Base.:(==)(e1::EwaldContext, e2::EwaldContext)
    e1.eframework == e2.eframework && e1.Eiks == e2.Eiks &&
                                      e1.allcharges == e2.allcharges &&
                                      e1.offsets == e2.offsets &&
                                      e1.static_contribution == e2.static_contribution &&
                                      e1.energy_net_charges && e2.energy_net_charges
end

# function EwaldContext(invmat::SMatrix{3,3,Float64,9})
#     eframework = EwaldFramework(invmat)
#     EwaldContext(eframework, ntuple(Returns(ComplexF64[]), 3), Vector{Float64}[], Int[], 0.0, 0.0)
# end

function move_one_system!(Eiks, ctx::EwaldContext, positions, ofs)
    Eikx, Eiky, Eikz = Eiks
    kx, ky, kz = ctx.eframework.kspace.ks
    kxp = kx + 1; tkyp = ky + ky + 1; tkzp = kz + kz + 1
    invmat = ctx.eframework.invmat
    isnumber = eltype(positions) <: AbstractVector{<:AbstractFloat}
    n = length(positions)
    @inbounds for j in 1:n
        px, py, pz = invmat * (isnumber ? positions[j] : NoUnits.(positions[j]/u"Å"))
        jofs = ofs + j
        make_line_pos!(Eikx, 1 + jofs*kxp, kx, cispi(2*px))
        make_line_neg!(Eiky, 1 + jofs*tkyp, ky, cispi(2*py))
        make_line_neg!(Eikz, 1 + jofs*tkzp, kz, cispi(2*pz))
    end
end
move_one_system!(ctx::EwaldContext, k::Int, positions) = move_one_system!(ctx.Eiks, ctx, positions, ctx.offsets[k])

"""
    EwaldContext(eframework::EwaldFramework, systems)

Build an [`EwaldContext`](@ref) for a fixed framework and any number of rigid molecules in
`systems`.
`eframework` can be obtained from [`initialize_ewald`](@ref).
"""
function EwaldContext(eframework::EwaldFramework, systems)
    iszero(eframework.α) && return EwaldContext(eframework, ntuple(Returns(ComplexF64[]), 3), Vector{Float64}[], Int[], 0.0, 0.0)
    allcharges::Vector{Vector{Float64}} = [NoUnits.(syst[:,:atomic_charge]/u"e_au") for syst in systems]

    chargefactor = (COULOMBIC_CONVERSION_FACTOR/sqrt(π))*eframework.α
    energy_adsorbate_self = sum(sum(abs2, charges; init=0.0)*chargefactor for charges in allcharges; init=0.0)
    net_charges = sum(sum(charges; init=0.0) for charges in allcharges; init=0.0)
    if abs(net_charges + eframework.net_charges_framework) > 1e-5
        @warn "Framework charge of $(eframework.net_charges_framework) and additional charge of $net_charges do not make a neutral system"
    end
    buffer, ortho, safemin = prepare_periodic_distance_computations(eframework.mat)
    buffer2 = MVector{3,Float64}(undef)
    m = length(systems)
    energy_adsorbate_excluded = 0.0
    for i in 1:m
        syst = systems[i]
        charges = allcharges[i]
        n = length(syst)
        for atomA in 1:n
            chargeA = charges[atomA]
            posA = Float64.(position(syst, atomA)/u"Å")
            for atomB in (atomA+1):n
                chargeB = charges[atomB]
                buffer2 .= Float64.(position(syst, atomB)/u"Å") .- posA
                mul!(buffer, eframework.invmat, buffer2)
                r = periodic_distance!(buffer, eframework.mat, ortho, safemin)
                energy_adsorbate_excluded += erf(eframework.α*r)*chargeA*chargeB/r
            end
        end
    end
    static_contribution = eframework.UIon*net_charges^2 - energy_adsorbate_self -
                          energy_adsorbate_excluded*COULOMBIC_CONVERSION_FACTOR

    energy_net_charges = eframework.UIon*eframework.net_charges_framework*net_charges
    offsets = Vector{Int}(undef, m)
    m > 0 && (offsets[1] = -1)
    @inbounds for i in 1:(m-1)
        offsets[i+1] = offsets[i] + length(systems[i])
    end
    Eiks = setup_Eik(systems, eframework.kspace.ks, eframework.invmat, (1,1,1))
    return EwaldContext(eframework, Eiks, allcharges, offsets, static_contribution, energy_net_charges)
end

"""
    compute_ewald(eframework::EwaldFramework, systems, skipcontribution=0)
    compute_ewald(ctx::EwaldContext, skipcontribution=0)

Compute the Fourier contribution to the Coulomb part of the interaction energy between
a framework and a set of rigid molecules.

If `skipcontribution` is set to `k`, the contributions of the charges of the `k`-th species
are not taken into account.
"""
function compute_ewald(ctx::EwaldContext, skipcontribution=0)
    iszero(ctx.eframework.α) && return 0.0u"K"
    newcharges::Vector{ComplexF64} = get!(task_local_storage(), :buffer_vector) do
        Vector{ComplexF64}(undef, ctx.eframework.kspace.num_kvecs)
    end
    resize!(newcharges, ctx.eframework.kspace.num_kvecs)
    ewald_main_loop!(newcharges, ctx.allcharges, ctx.eframework.kspace, ctx.Eiks, skipcontribution)
    framework_adsorbate = 0.0
    adsorbate_adsorbate = 0.0
    for i_kvec in 1:ctx.eframework.kspace.num_kvecs
        temp = ctx.eframework.kfactors[i_kvec]
        _re_f, _im_f = reim(ctx.eframework.StoreRigidChargeFramework[i_kvec])
        _re_a, _im_a = reim(newcharges[i_kvec])
        framework_adsorbate += temp*(_re_f*_re_a + _im_f*_im_a)
        adsorbate_adsorbate += temp*(_re_a*_re_a + _im_a*_im_a)
    end

    UHostAdsorbateChargeChargeFourier = 2*(framework_adsorbate + ctx.energy_net_charges)

    UAdsorbateAdsorbateChargeChargeFourier = adsorbate_adsorbate + ctx.static_contribution

    return (UHostAdsorbateChargeChargeFourier + UAdsorbateAdsorbateChargeChargeFourier)*ENERGY_TO_KELVIN
end

function compute_ewald(eframework::EwaldFramework, systems, skipcontribution=0)
    compute_ewald(EwaldContext(eframework, systems), skipcontribution)
end


struct IncrementalEwaldContext
    ctx::EwaldContext
    sums::Matrix{ComplexF64} # need to change to ElasticArray if number of charges can change
    tmpEiks::NTuple{3,Vector{ComplexF64}}
    tmpsums::Vector{ComplexF64}
    last::Base.RefValue{Int} # last inquired change
end
function IncrementalEwaldContext(ctx::EwaldContext)
    n = length(ctx.allcharges)
    m = ctx.eframework.kspace.num_kvecs
    kx, ky, kz = ctx.eframework.kspace.ks
    kxp = kx+1; tkyp = ky+ky+1; tkzp = kz+kz+1
    Eikx = Vector{ComplexF64}(undef, kxp)
    Eiky = Vector{ComplexF64}(undef, tkyp)
    Eikz = Vector{ComplexF64}(undef, tkzp)
    tmpsums = Vector{ComplexF64}(undef, m)
    IncrementalEwaldContext(ctx, Matrix{ComplexF64}(undef, m, n+1), (Eikx, Eiky, Eikz), tmpsums, Ref(-1))
end

"""
    compute_ewald(ewald::IncrementalEwaldContext)

See [`compute_ewald(ctx::EwaldContext)`](@ref).

!!! warning
    When used with an `IncrementalEwaldContext`, `compute_ewald` modifies a hidden state of
    its input. This means that it cannot be called from multiple threads in parallel on the
    same input.
"""
function compute_ewald(ewald::IncrementalEwaldContext)
    iszero(ewald.ctx.eframework.α) && return 0.0u"K"
    ewald_main_loop!(ewald.sums, ewald.ctx.allcharges, ewald.ctx.eframework.kspace, ewald.ctx.Eiks, 0)
    framework_adsorbate = 0.0
    adsorbate_adsorbate = 0.0
    for idxkvec in 1:ewald.ctx.eframework.kspace.num_kvecs
        adsorbate_contribution = ewald.sums[idxkvec,end]
        _re_a, _im_a = reim(adsorbate_contribution)
        temp = ewald.ctx.eframework.kfactors[idxkvec]
        _re_f, _im_f = reim(ewald.ctx.eframework.StoreRigidChargeFramework[idxkvec])
        framework_adsorbate += temp*(_re_f*_re_a + _im_f*_im_a)
        adsorbate_adsorbate += temp*(_re_a*_re_a + _im_a*_im_a)
    end

    UHostAdsorbateChargeChargeFourier = 2*(framework_adsorbate + ewald.ctx.energy_net_charges)

    UAdsorbateAdsorbateChargeChargeFourier = adsorbate_adsorbate + ewald.ctx.static_contribution

    ewald.last[] = 0 # signal for single_contribution_ewald

    (UHostAdsorbateChargeChargeFourier + UAdsorbateAdsorbateChargeChargeFourier)*ENERGY_TO_KELVIN
end


"""
    single_contribution_ewald(ewald::IncrementalEwaldContext, k, positions)

Compute the contribution of species number `k` at the given `positions` (one position per
atom) to the reciprocal Ewald sum, so that
`single_contribution_ewald(ewald, k, poss2) - single_contribution_ewald(ewald, k, poss1)`
is the reciprocal Ewald sum energy difference between species `k` at positions `poss2` and
`poss1`.

If `positions == nothing`, use the position for the species currently stored in `ewald`
(note that this particular computation is substantially faster).

!!! info
    This function can only be called if `compute_ewald(ewald)` has been called prior.
    Otherwise it will error.
    This is necessary to ensure that the `ewald.sums` are correctly set.

!!! warning
    This function is not thread-safe if `positions !== nothing`. This means that it must
    not be called from multiple threads in parallel on the same `ewald` with
    `positions !== nothing`.
"""
function single_contribution_ewald(ewald::IncrementalEwaldContext, k, positions)
    iszero(ewald.ctx.eframework.α) && return 0.0u"K"
    ewald.last[] == -1 && error("Please call `compute_ewald(ewald)` before calling `single_contribution_ewald(ewald, ...)`")
    if positions isa Nothing
        contribution = @view ewald.sums[:,k]
    else
        kx, ky, kz = ewald.ctx.eframework.kspace.ks
        kxp = kx+1; tkyp = ky+ky+1; tkzp = kz+kz+1
        Eiks = ewald.tmpEiks
        Eikx, Eiky, Eikz = Eiks
        m = length(positions)
        resize!(Eikx, m*kxp); resize!(Eiky, m*tkyp); resize!(Eikz, m*tkzp)
        move_one_system!(Eiks, ewald.ctx, positions, -1) # this only modifies Eiks, not ewald
        contribution = ewald.tmpsums
        contribution .= zero(ComplexF64)
        charges = ewald.ctx.allcharges[k]
        kindices = ewald.ctx.eframework.kspace.kindices
        ix = 0
        iy = ky+1
        iz = kz+1
        for c in charges
            for (jy, jz, jxrange, rangeidx) in kindices
                Eik_xy = c*Eiky[iy+jy]*Eikz[iz+jz]
                n = length(jxrange)
                ofs = first(jxrange) + ix
                @simd for I in 1:n
                    contribution[rangeidx+I] += Eikx[ofs+I]*Eik_xy
                end
            end
            ix += kxp
            iy += tkyp
            iz += tkzp
        end
        ewald.last[] = k # signal for update_ewald_context!
    end

    rest_single = 0.0
    single_single = 0.0
    n = length(ewald.ctx.allcharges)
    @inbounds for idxkvec in 1:ewald.ctx.eframework.kspace.num_kvecs
        rest = ewald.ctx.eframework.StoreRigidChargeFramework[idxkvec]
        rest += ewald.sums[idxkvec,end] - ewald.sums[idxkvec,k]
        _re_f, _im_f = reim(rest)
        _re_a, _im_a = reim(contribution[idxkvec])
        temp = ewald.ctx.eframework.kfactors[idxkvec]
        rest_single += temp*(_re_f*_re_a + _im_f*_im_a)
        single_single += temp*(_re_a*_re_a + _im_a*_im_a)
    end

    (2*rest_single + single_single)*ENERGY_TO_KELVIN
end

"""
    update_ewald_context!(ewald::IncrementalEwaldContext)

Following a call to [`single_contribution_ewald(ewald, k, positions)`](@ref), update the
internal state of `ewald` so that the species `k` is now at the given `positions`.

Note that this modifies the underlying [`EwaldContext`](@ref) `ewald.ctx`.
"""
function update_ewald_context!(ewald::IncrementalEwaldContext)
    k = ewald.last[]
    k ≤ 0 && return # last moved molecule bears no charge TODO: check
    ewald.sums[:,end] .+= ewald.tmpsums .- @view(ewald.sums[:,k])
    ewald.sums[:,k] .= ewald.tmpsums
    Eikx, Eiky, Eikz = ewald.ctx.Eiks
    newEikx, newEiky, newEikz = ewald.tmpEiks
    kx, ky, kz = ewald.ctx.eframework.kspace.ks
    chargepos = 1 + ewald.ctx.offsets[k]
    copyto!(Eikx, 1 + chargepos*(kx+1), newEikx, 1, length(newEikx))
    copyto!(Eiky, 1 + chargepos*(ky+ky+1), newEiky, 1, length(newEiky))
    copyto!(Eikz, 1 + chargepos*(kz+kz+1), newEikz, 1, length(newEikz))
    ewald.last[] = 0
    nothing
end
