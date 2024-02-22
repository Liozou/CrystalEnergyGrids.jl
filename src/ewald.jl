using LinearAlgebra: norm, det
using StaticArrays
using AtomsBase
using SpecialFunctions: erf, erfc
using ElasticArrays: ElasticMatrix

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
    α::typeof(1.0u"Å^-1")
    mat::SMatrix{3,3,TÅ,9}
    invmat::SMatrix{3,3,typeof(1.0u"Å^-1"),9}
    kfactors::Vector{Float64}
    UIon::typeof(1.0u"K/e_au^2")
    StoreRigidChargeFramework::Vector{ComplexF64}
    net_charges_framework::Te_au
    precision::Float64
end
function EwaldFramework(mat::AbstractMatrix)
    cell = CellMatrix(mat)
    EwaldFramework(EwaldKspace((0,0,0), 0, Tuple{Int,Int,UnitRange{Int},Int}[]),
                   0.0u"Å^-1", cell.mat, cell.invmat, Float64[], 0.0u"K/e_au^2", ComplexF64[], 0.0u"e_au", 0.0)
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


function make_line_pos!(Eikt, j, eikt)
    k = size(Eikt, 1)
    Eikt[1,j] = 1
    Eikt[2,j] = eikt
    @inbounds for i in 3:k
        Eikt[i,j] = Eikt[i-1,j]*eikt
    end
    nothing
end

function make_line_neg!(Eikt, j, k, eikt)
    ceikt = conj(eikt)
    Eikt[k,j] = ceikt
    @inbounds for i in k-1:-1:1
        Eikt[i,j] = Eikt[i+1,j]*ceikt
    end
    make_line_pos!(@view(Eikt[k+1:end,j]), 1, eikt)
end


"""
    setup_Eik(systems, (kx, ky, kz), invmat, (ΠA, ΠB, ΠC), numspecies)

Setup the computation of the e^(ik...) necessary for Ewald summation, for the given
`systems`, with `numsites` total number of point charges (including those in the supercell)
along a `kx × ky × kz` K-space grid.
`invmat` is the inverse of the real-space unit cell. That unit cell is multiplied by
`ΠA × ΠB × ΠC` to obtain the supercell used for the computation.

Return `(Eikx, Eiky, Eikz)` where, for `m` in `1:numsites`
- `Eikx[(m-1)*(kx+1) + 1 + k] == exp(2ikπ * fracpos_x[m])` with `k ∈ 0:kx`
- `Eiky[(m-1)*(2ky+1) + ky + 1 + k] == exp(±2ikπ * fracpos_y[m])` with `k ∈ -ky:ky`
- `Eikz[(m-1)*(2ky+1) + ky + 1 + k] == exp(±2ikπ * fracpos_z[m])` with `k ∈ -kz:kz`

- `Eikx[k,m] == exp(2ikπ * fracpos_x[m])` with `k ∈ 0:kx`
- `Eiky[k,m] == exp(±2ikπ * fracpos_y[m])` with `k ∈ -ky:ky`
- `Eikz[k,m] == exp(±2ikπ * fracpos_z[m])` with `k ∈ -kz:kz`

with `fracpos_x[m], fracpos_y[m], fracpos_z[m] == invmat*position[m]`
"""
function setup_Eik(systems, (kx, ky, kz), invmat, (ΠA, ΠB, ΠC), numspecies)
    ΠBC = ΠB*ΠC
    ΠABC = ΠA*ΠBC
    numsites = ΠABC * sum(sum(length(kind[i_species]) for i_species in 1:num; init=0)::Int for (num, kind) in zip(numspecies, systems); init=0)::Int

    kxp = kx + 1
    tkyp = ky + ky + 1
    tkzp = kz + kz + 1

    Eikx = ElasticMatrix{ComplexF64}(undef, kxp, numsites)
    Eiky = ElasticMatrix{ComplexF64}(undef, tkyp, numsites)
    Eikz = ElasticMatrix{ComplexF64}(undef, tkzp, numsites)

    ofs = -1
    @inbounds for (i_kind, kind) in enumerate(systems), i_syst in 1:numspecies[i_kind]
        syst = kind[i_syst]
        for j in 1:length(syst)
            px, py, pz = invmat * position(syst, j)
            for πs in CartesianIndices((0:(ΠC-1), 0:(ΠB-1), 0:(ΠA-1)))
                πc, πb, πa = Tuple(πs)
                jofs = 1 + (j + ofs)*ΠABC + πa*ΠBC + πb*ΠC + πc

                eikx = cispi(2*(px + πa/ΠA))
                make_line_pos!(Eikx, jofs, eikx)

                eiky = cispi(2*(py + πb/ΠB))
                make_line_neg!(Eiky, jofs, ky, eiky)

                eikz = cispi(2*(pz + πc/ΠC))
                make_line_neg!(Eikz, jofs, kz, eikz)
            end
        end
        ofs += length(syst)
    end

    (Eikx, Eiky, Eikz)
end


function ewald_main_loop!(sums::T, ctx, kspace::EwaldKspace, skipcontribution) where T
    fill!(sums, zero(ComplexF64))
    Eikx, Eiky, Eikz = ctx.Eiks
    _, ky, kz = kspace.ks

    atoms = ctx.atoms
    speciesof = ctx.speciesof
    charges = ctx.charges

    for (chargecounter, (ij,k)) in enumerate(atoms)
        skipcontribution == ij && continue
        c = ustrip(u"e_au", charges[speciesof[ij]][k])
        ijp1 = ij + 1
        # if !iszero(c) # TODO: remove zero charges when constructing EwaldContext
            for (jy, jz, jxrange, rangeidx) in kspace.kindices
                Eik_yz = c*Eiky[ky+1+jy,chargecounter]*Eikz[kz+1+jz,chargecounter]
                n = length(jxrange)
                ofs = first(jxrange)
                if T === Vector{ComplexF64}
                    @simd for I in 1:n
                        sums[rangeidx+I] += Eikx[ofs+I,chargecounter]*Eik_yz
                    end
                else
                    @simd for I in 1:n
                        sums[rangeidx+I,ijp1] += Eikx[ofs+I,chargecounter]*Eik_yz
                    end
                end
            end
        # end
    end
    if T !== Vector{ComplexF64}
        sums[:,1] .= sum(sums; dims=2)
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

    cutoff_coulomb = 12.0u"Å"
    # From a dimensional point of view, the definition of tol is weird... We follow the
    # code of RASPA.
    ε = ustrip(u"Å", cutoff_coulomb*min(0.5, abs(precision)))
    tol = sqrt(abs(log(ε)))
    α = sqrt(abs(log(ε*tol)))/cutoff_coulomb
    tol1 = sqrt(-log(ε*4.0*(tol*α*u"Å")^2))

    mat = stack3(bounding_box(syst) .* supercell)
    len_a, len_b, len_c = cell_lengths(mat)
    __α = α*u"Å"*tol1/π
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
    volume_factor = COULOMBIC_CONVERSION_FACTOR*2π/(det(mat)*u"Å^3")
    @assert volume_factor > 0*u"K/e_au^2/Å^2"
    α_factor = -0.25/(α*u"Å")^2

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
            kfactors[rangeidx+I] = volume_factor*(1+(i!=0))*exp(α_factor*rksqr)/rksqr/u"K/e_au^2/Å^2"
        end
    end

    UIon = COULOMBIC_CONVERSION_FACTOR*α/sqrt(π) - sum(kfactors; init=0.0)*u"K/e_au^2"

    Π = prod(supercell)
    charges::Vector{Te_au} = repeat(syst[:,:atomic_charge]; inner=Π)

    kspace = EwaldKspace((kx, ky, kz), num_kvecs, kindices)
    Eiks = setup_Eik(((syst,),), kspace.ks, invmat*u"Å^-1", supercell, (1,))
    StoreRigidChargeFramework = Vector{ComplexF64}(undef, num_kvecs)
    atoms = [(1,k) for k in 1:length(charges)]
    pseudoctx = (; atoms, Eiks, charges=[charges], speciesof=[1])
    ewald_main_loop!(StoreRigidChargeFramework, pseudoctx, kspace, 0)
    # UChargeChargeFrameworkRigid = 0.0
    # for idx in 1:num_kvecs
    #     UChargeChargeFrameworkRigid += kfactors[idx]*abs2(StoreRigidChargeFramework[idx])
    # end
    net_charges_framework = sum(charges; init=0.0u"e_au")
    # UChargeChargeFrameworkRigid += UIon*net_charges_framework^2

    return EwaldFramework(kspace, α, mat*u"Å", invmat*u"Å^-1", kfactors, UIon,
                          StoreRigidChargeFramework, net_charges_framework, precision)
end

"""
    initialize_ewald(mat::AbstractMatrix, supercell=find_supercell(mat, 12.0u"Å"), precision=1e-6)

Return an `EwaldFramework` object representing the empty system of dimensions given by its
unit cell matrix `mat` and optionally a supercell.

See also [initialize_ewald(syst::AbstractSystem{3}, supercell=find_supercell(syst, 12.0u"Å"), precision=1e-6)](@ref)
"""
function initialize_ewald(mat::AbstractMatrix, supercell=find_supercell(mat, 12.0u"Å"), precision=1e-6)
    col1, col2, col3 = eltype(mat) <: Quantity ? broadcast.(Base.Fix1(ustrip, u"Å"), eachcol(mat)) : eachcol(mat)
    smat = SA[SVector{3,Float64}(col1)*u"Å", SVector{3,Float64}(col2)*u"Å", SVector{3,Float64}(col3)*u"Å"]
    emptysystem = RASPASystem(smat, SVector{3,TÅ}[], Symbol[], Int[], typeof(1.0u"u")[], Te_au[], false)
    initialize_ewald(emptysystem, supercell, precision)
end


function derivatives_ewald(ewald::EwaldFramework, charge::Float64, r2::Float64)
    r = sqrt(r2)
    r3 = r2*r
    r5 = r3*r2
    α = ewald.α*u"Å"
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
    Eiks::NTuple{3,ElasticMatrix{ComplexF64,Vector{ComplexF64}}}
    # allcharges::Vector{Vector{Float64}}
    atoms::Vector{Tuple{Int,Int}}
    # atoms[l] is (ij,k) where `ij` is the index of the specific molecule and `k` is the
    # atom number within the molecule
    speciesof::Vector{Int} # speciesof[ij] is the kind i of the corresponding species
    indexof::Vector{Vector{Int}}
    # indexof[ij][k] is the index `l` such that atoms[l] == (ij, k)
    numspecies::Vector{Int} # numspecies[i] is the number of ij such that speciesof[ij]==i
    static_contribution::Base.RefValue{TK}
    charges::Vector{Vector{Te_au}} # charges[i] is the list of atom charges of species i
    energies::Vector{TK} # energies[i] is the negative energy contribution to
    # static_contribution of each species of kind i
    energy_net_charges::TK
end

function Base.:(==)(e1::EwaldContext, e2::EwaldContext)
    e1.eframework == e2.eframework && e1.Eiks == e2.Eiks &&
                                      e1.atoms == e2.atoms &&
                                      e1.charges == e2.charges &&
                                      e1.speciesof == e2.speciesof &&
                                      # no need to compare indexof because of atoms
                                      # no need to compare numspecies because of speciesof
                                      e1.static_contribution[] ≈ e2.static_contribution[] &&
                                      e1.energy_net_charges ≈ e2.energy_net_charges
end

function move_one_system!(Eiks, ctx::EwaldContext, ij::Union{Nothing,Int}, positions)
    Eikx, Eiky, Eikz = Eiks
    _, ky, kz = ctx.eframework.kspace.ks
    invmat = ctx.eframework.invmat
    isnumber = eltype(positions) <: AbstractVector{<:AbstractFloat}
    n = length(positions)
    indexofij = ij isa Nothing ? nothing : ctx.indexof[ij]
    @inbounds for k in 1:n
        px, py, pz = invmat * (isnumber ? positions[k]*u"Å" : positions[k])
        ofs = ij isa Nothing ? k : indexofij[k]
        make_line_pos!(Eikx, ofs, cispi(2*px))
        make_line_neg!(Eiky, ofs, ky, cispi(2*py))
        make_line_neg!(Eikz, ofs, kz, cispi(2*pz))
    end
end

"""
    move_one_system!(ctx::EwaldContext, ij::Int, positions)

Move one molecule given by its index `ij` in the [`EwaldContext`](@ref) to the given atom
`positions`.
"""
move_one_system!(ctx::EwaldContext, ij::Int, positions) = move_one_system!(ctx.Eiks, ctx, ij, positions)

"""
    add_one_system!(ctx::EwaldContext, i::Int, positions)

Add to the [`EwaldContext`](@ref) one system defined by its kind `i` and its atom
`positions`.

Return the index `ij` of the newly introduced molecule.
"""
function add_one_system!(ctx::EwaldContext, i::Int, positions)
    ctx.static_contribution[] -= ctx.energies[i]
    push!(ctx.speciesof, i)
    ctx.numspecies[i] += 1
    ij = length(ctx.speciesof)
    l0 = length(ctx.atoms)
    n = length(positions)
    push!(ctx.indexof, collect((l0+1):(l0+n)))
    oldn = length(ctx.atoms)
    resize!(ctx.atoms, oldn + n)
    for k in 1:n
        ctx.atoms[l0+k] = (ij,k)
    end
    Eikx, Eiky, Eikz = ctx.Eiks
    sx = size(Eikx, 1); sy = size(Eiky, 1); sz = size(Eikz, 1)
    resize!(Eikx, sx, oldn + n); resize!(Eiky, sy, oldn + n); resize!(Eikz, sz, oldn + n)
    move_one_system!(ctx.Eiks, ctx, ij, positions)
    ij
end

"""
    remove_one_system!(ctx::EwaldContext, ij::Int)

Remove from the [`EwaldContext`](@ref) one system identified by its index `ij`, and return
the index `oldij` of the system which is now referred to as `ij`. It may return
`oldij == ij` to signify that there is no change to the current indices.

In other words, if `remove_one_system!(ctx, 7)` returns `12`, then the system which used to
be called `12` is now called `7`, and there is no system which is now called `12`.

## Implementation note
The value returned by `remove_one_system!` is always the highest index referring to a
system in `ctx`.
"""
function remove_one_system!(ctx::EwaldContext, ij::Int)
    isempty(ctx.charges) && return 0

    i = ctx.speciesof[ij]
    ctx.static_contribution[] += ctx.energies[i]
    ctx.numspecies[i] -= 1
    prev_atoms = sort!(ctx.indexof[ij])
    oldij = length(ctx.speciesof)
    if oldij != ij
        ctx.speciesof[ij] = ctx.speciesof[oldij]
        for l in (ctx.indexof[ij] = ctx.indexof[oldij])
            _oldij, _k = ctx.atoms[l]
            @assert _oldij == oldij
            ctx.atoms[l] = (ij, _k)
        end
    end
    pop!(ctx.speciesof)
    pop!(ctx.indexof)

    #=
    Put removed atoms at the end. To do so, exchange atom l with atom first_l, where l goes
    decreasing in the indices of ctx.atoms, and first_l goes increasing through prev_atoms.
    If l == last_l, skip since it means that the removed atom of index last_l is among the
    last indices of atoms.
    =#
    n = length(prev_atoms)
    oldn = length(ctx.atoms)
    lastidx = n
    firstidx = 1
    Eikx, Eiky, Eikz = ctx.Eiks
    for k in 1:n
        l = oldn + 1 - k
        if l == prev_atoms[lastidx]
            lastidx -= 1
        else
            first_l = prev_atoms[firstidx]
            (atomij, atomk) = ctx.atoms[l]
            ctx.indexof[atomij][atomk] = first_l
            ctx.atoms[first_l] = (atomij, atomk)
            Eikx[:,first_l] .= Eikx[:,l]
            Eiky[:,first_l] .= Eiky[:,l]
            Eikz[:,first_l] .= Eikz[:,l]
            firstidx += 1
        end
    end

    resize!(ctx.atoms, oldn - n)
    sx = size(Eikx, 1); sy = size(Eiky, 1); sz = size(Eikz, 1)
    resize!(Eikx, sx, oldn - n); resize!(Eiky, sy, oldn - n); resize!(Eikz, sz, oldn - n)
    return oldij
end

"""
    EwaldContext(eframework::EwaldFramework, systems)

Build an [`EwaldContext`](@ref) for a fixed framework and any number of rigid molecules in
`systems`.
`eframework` can be obtained from [`initialize_ewald`](@ref).
"""
function EwaldContext(eframework::EwaldFramework, systems, emptysystems=())
    iszero(eframework.α) && return EwaldContext(eframework, ntuple(Returns(ElasticMatrix{ComplexF64}(undef, 0, 0)), 3), Tuple{Int,Int}[], Int[], Vector{Int}[], Int[], Ref(0.0u"K"), Vector{Vector{Te_au}}[], TK[], 0.0u"K")
    allcharges::Vector{Vector{Te_au}} = [first(kind)[:,:atomic_charge] for kind in systems]
    numspecies::Vector{Int} = [length(kind) for kind in systems]
    for i_empty in emptysystems
        @assert numspecies[i_empty] == 1
        numspecies[i_empty] = 0
    end

    chargefactor = COULOMBIC_CONVERSION_FACTOR/sqrt(π)*eframework.α
    energies = [sum(abs2, charges; init=0.0u"e_au^2")*chargefactor for charges in allcharges]
    energy_adsorbate_self = sum(eas*num for (num, eas) in zip(numspecies, energies); init=0.0u"K")
    net_charges = sum.(allcharges; init=0.0u"e_au")
    total_net_charges = sum(charge*num for (num, charge) in zip(numspecies, net_charges); init=0.0u"e_au")
    if abs(total_net_charges + eframework.net_charges_framework) > 1e-5u"e_au" && PRINT_CHARGE_WARNING[]
        @warn lazy"Framework charge of $(eframework.net_charges_framework) and additional charge of $total_net_charges do not make a neutral system."
    end

    buffer2, ortho, safemin = prepare_periodic_distance_computations(eframework.mat)
    safemin2 = safemin^2
    buffer = MVector{3,TÅ}(undef)
    # m = length(systems)
    energy_adsorbate_excluded = 0.0u"K"
    speciesof = zeros(Int, sum(numspecies))
    atoms = Tuple{Int,Int}[]
    indexof = Vector{Vector{Int}}(undef, length(speciesof))
    systcounter = 0
    for (i, (num, kind)) in enumerate(zip(numspecies, systems))
        num == 0 && continue
        for syst in kind
            systcounter += 1
            speciesof[systcounter] = i
            n = length(syst)
            indexof[systcounter] = collect((length(atoms)+1):(length(atoms)+n))
            append!(atoms, (systcounter, k) for k in 1:n)
        end
    end
    for (i, (num, kind, charges)) in enumerate(zip(numspecies, systems, allcharges))
        syst = first(kind)
        this_energy = 0.0u"e_au^2/Å"
        n = length(syst)
        for atomA in 1:n
            chargeA = charges[atomA]
            posA = position(syst, atomA)
            for atomB in (atomA+1):n
                chargeB = charges[atomB]
                buffer .= position(syst, atomB) .- posA
                r = sqrt(periodic_distance2_fromcartesian!(buffer, eframework.mat, eframework.invmat, ortho, safemin2, buffer2))
                this_energy += erf(eframework.α*r)*chargeA*chargeB/r
            end
        end
        energies[i] += this_energy*COULOMBIC_CONVERSION_FACTOR
        energy_adsorbate_excluded += num*this_energy*COULOMBIC_CONVERSION_FACTOR
    end
    @assert !any(iszero, speciesof)

    static_contribution = eframework.UIon*total_net_charges^2 - energy_adsorbate_self - energy_adsorbate_excluded

    energy_net_charges = eframework.UIon*eframework.net_charges_framework*total_net_charges
    # offsets = Vector{Int}(undef, m)
    # m > 0 && (offsets[1] = 0)
    # @inbounds for i in 1:(m-1)
    #     offsets[i+1] = offsets[i] + length(systems[i])
    # end
    Eiks = setup_Eik(systems, eframework.kspace.ks, eframework.invmat, (1,1,1), numspecies)
    return EwaldContext(eframework, Eiks, atoms, speciesof, indexof, numspecies, Ref(static_contribution), allcharges, energies, energy_net_charges)
end

isdefined_ewald(ctx::EwaldContext) = !iszero(ctx.eframework.α)

"""
    compute_ewald(eframework::EwaldFramework, systems, skipcontribution=0)
    compute_ewald(ctx::EwaldContext, skipcontribution=0)

Compute the Fourier contribution to the Coulomb part of the interaction energy between
a framework and a set of rigid molecules.

If `skipcontribution` is set to `ij`, the contributions of the charges of the `ij`-th species
are not taken into account.
"""
function compute_ewald(ctx::EwaldContext, skipcontribution=0)
    isdefined_ewald(ctx) || return 0.0u"K"
    newcharges::Vector{ComplexF64} = get!(task_local_storage(), :CEG_buffer_ewald_context) do
        Vector{ComplexF64}(undef, ctx.eframework.kspace.num_kvecs)
    end
    resize!(newcharges, ctx.eframework.kspace.num_kvecs)
    ewald_main_loop!(newcharges, ctx, ctx.eframework.kspace, skipcontribution)
    framework_adsorbate = 0.0
    adsorbate_adsorbate = 0.0
    for i_kvec in 1:ctx.eframework.kspace.num_kvecs
        temp = ctx.eframework.kfactors[i_kvec]
        _re_f, _im_f = reim(ctx.eframework.StoreRigidChargeFramework[i_kvec])
        _re_a, _im_a = reim(newcharges[i_kvec])
        framework_adsorbate += temp*(_re_f*_re_a + _im_f*_im_a)
        adsorbate_adsorbate += temp*(_re_a*_re_a + _im_a*_im_a)
    end

    UHostAdsorbateChargeChargeFourier = 2*(framework_adsorbate*u"K" + ctx.energy_net_charges)

    UAdsorbateAdsorbateChargeChargeFourier = adsorbate_adsorbate*u"K" + ctx.static_contribution[]

    return UHostAdsorbateChargeChargeFourier + UAdsorbateAdsorbateChargeChargeFourier
end

function compute_ewald(eframework::EwaldFramework, systems, skipcontribution=0)
    compute_ewald(EwaldContext(eframework, systems), skipcontribution)
end


struct IncrementalEwaldContext
    ctx::EwaldContext
    sums::ElasticMatrix{ComplexF64,Vector{ComplexF64}}
    tmpEiks::NTuple{3,ElasticMatrix{ComplexF64,Vector{ComplexF64}}}
    tmpsums::Vector{ComplexF64}
    last::Base.RefValue{Int} # last inquired change
end
function IncrementalEwaldContext(ctx::EwaldContext)
    n = sum(ctx.numspecies)
    m = ctx.eframework.kspace.num_kvecs
    kx, ky, kz = ctx.eframework.kspace.ks
    kxp = kx+1; tkyp = ky+ky+1; tkzp = kz+kz+1
    Eikx = ElasticMatrix{ComplexF64}(undef, kxp, 1)
    Eiky = ElasticMatrix{ComplexF64}(undef, tkyp, 1)
    Eikz = ElasticMatrix{ComplexF64}(undef, tkzp, 1)
    tmpsums = Vector{ComplexF64}(undef, m)
    IncrementalEwaldContext(ctx, ElasticMatrix{ComplexF64}(undef, m, (m!=0)*(n+1)), (Eikx, Eiky, Eikz), tmpsums, Ref(-isdefined_ewald(ctx)))
end

isdefined_ewald(ewald::IncrementalEwaldContext) = isdefined_ewald(ewald.ctx)

"""
    compute_ewald(ewald::IncrementalEwaldContext)

See [`compute_ewald(ctx::EwaldContext)`](@ref).

!!! warning
    When used with an `IncrementalEwaldContext`, `compute_ewald` modifies a hidden state of
    its input. This means that it cannot be called from multiple threads in parallel on the
    same input.
"""
function compute_ewald(ewald::IncrementalEwaldContext)
    isdefined_ewald(ewald) || return 0.0u"K"
    ewald_main_loop!(ewald.sums, ewald.ctx, ewald.ctx.eframework.kspace, 0)
    framework_adsorbate = 0.0
    adsorbate_adsorbate = 0.0
    for idxkvec in 1:ewald.ctx.eframework.kspace.num_kvecs
        adsorbate_contribution = ewald.sums[idxkvec,1]
        _re_a, _im_a = reim(adsorbate_contribution)
        temp = ewald.ctx.eframework.kfactors[idxkvec]
        _re_f, _im_f = reim(ewald.ctx.eframework.StoreRigidChargeFramework[idxkvec])
        framework_adsorbate += temp*(_re_f*_re_a + _im_f*_im_a)
        adsorbate_adsorbate += temp*(_re_a*_re_a + _im_a*_im_a)
    end

    UHostAdsorbateChargeChargeFourier = 2*(framework_adsorbate*u"K" + ewald.ctx.energy_net_charges)

    UAdsorbateAdsorbateChargeChargeFourier = adsorbate_adsorbate*u"K" + ewald.ctx.static_contribution[]

    ewald.last[] = 0 # signal for single_contribution_ewald

    UHostAdsorbateChargeChargeFourier + UAdsorbateAdsorbateChargeChargeFourier
end


function update_sums!(sums, ewald::IncrementalEwaldContext, ij::Int, (Eikx, Eiky, Eikz), signal)
    sums .= zero(ComplexF64)
    charges = ewald.ctx.charges[ij < 0 ? -ij : ewald.ctx.speciesof[ij]]
    kindices = ewald.ctx.eframework.kspace.kindices
    _, ky, kz = ewald.ctx.eframework.kspace.ks
    for (chargecounter, c) in enumerate(charges)
        for (jy, jz, jxrange, rangeidx) in kindices
            Eik_xy = ustrip(u"e_au", c)*Eiky[jy+ky+1,chargecounter]*Eikz[jz+kz+1,chargecounter]
            n = length(jxrange)
            ofs = first(jxrange)
            @simd for I in 1:n
                sums[rangeidx+I] += Eikx[ofs+I,chargecounter]*Eik_xy
            end
        end
    end
    if signal
        ewald.last[] = ij # signal for update_ewald_context!
    end
    nothing
end

"""
    single_contribution_ewald(ewald::IncrementalEwaldContext, ij, positions; tmpEiks::NTuple{3,<:AbstractMatrix{ComplexF64}}=ewald.tmpEiks, tmpsums::AbstractVector{ComplexF64}=ewald.tmpsums)

Compute the contribution of species number `ij` at the given `positions` (one position per
atom) to the reciprocal Ewald sum, so that
`single_contribution_ewald(ewald, ij, poss2) - single_contribution_ewald(ewald, ij, poss1)`
is the reciprocal Ewald sum energy difference between species `ij` at positions `poss2` and
`poss1`.

If `positions == nothing`, use the position for the species currently stored in `ewald`
(note that this particular computation is substantially faster).

If `ij < 0`, compute the contribution of a new species of kind `-ij` at the given
`positions` (which cannot be `nothing`) to the reciprocal Ewald sum.

!!! info
    This function can only be called if `compute_ewald(ewald)` has been called prior.
    Otherwise it will error.
    This is necessary to ensure that the `ewald.sums` are correctly set.

!!! warning
    This function is not thread-safe if `positions !== nothing`. This means that it must
    not be called from multiple threads in parallel on the same `ewald` with
    `positions !== nothing`.

    !!! note
        It is possible to make it thread-safe even with `positions !== nothing`: to do so,
        you must pass thread-local buffers as keyword arguments `tmpEiks` and `tmpsums`.
        Note that in that case, this function cannot be used with `update_ewald_context!`,
        see the following danger.

!!! danger
    Specifying a value to either keyword arguments `tmpEiks` or `tmpsums` will prevent the
    state-local change required to perform [`update_ewald_context!`](@ref) afterwards.
    Calling [`update_ewald_context`](@ref) will thus not perform any actual update, but may
    not error either.
"""
function single_contribution_ewald(ewald::IncrementalEwaldContext, ij, positions;
                                   tmpEiks::NTuple{3,<:AbstractMatrix{ComplexF64}}=ewald.tmpEiks,
                                   tmpsums::AbstractVector{ComplexF64}=ewald.tmpsums
                                  )
    isdefined_ewald(ewald) || return 0.0u"K"
    ewald.last[] == -1 && error("Please call `compute_ewald(ewald)` before calling `single_contribution_ewald(ewald, ...)`")
    if positions isa Nothing
        contribution = @view ewald.sums[:,ij+1]
    else
        kx, ky, kz = ewald.ctx.eframework.kspace.ks
        kxp = kx+1; tkyp = ky+ky+1; tkzp = kz+kz+1
        Eikx, Eiky, Eikz = tmpEiks
        m = length(positions)
        resize!(Eikx, kxp, m); resize!(Eiky, tkyp, m); resize!(Eikz, tkzp, m)
        move_one_system!(tmpEiks, ewald.ctx, nothing, positions) # this only modifies tmpEiks, not ewald
        contribution = tmpsums
        update_sums!(contribution, ewald, ij, tmpEiks, (ij > 0) & (tmpEiks===ewald.tmpEiks) & (tmpsums===ewald.tmpsums))
    end

    rest_single = 0.0
    single_single = 0.0
    @inbounds for idxkvec in 1:ewald.ctx.eframework.kspace.num_kvecs
        rest = ewald.ctx.eframework.StoreRigidChargeFramework[idxkvec]
        if ij < 0
            rest += ewald.sums[idxkvec,1]
        else
            rest += ewald.sums[idxkvec,1] - ewald.sums[idxkvec,ij+1]
        end
        _re_f, _im_f = reim(rest)
        _re_a, _im_a = reim(contribution[idxkvec])
        temp = ewald.ctx.eframework.kfactors[idxkvec]
        rest_single += temp*(_re_f*_re_a + _im_f*_im_a)
        single_single += temp*(_re_a*_re_a + _im_a*_im_a)
    end

    (2*rest_single + single_single)*u"K"
end

"""
    update_ewald_context!(ewald::IncrementalEwaldContext)

Following a call to [`single_contribution_ewald(ewald, ij, positions)`](@ref), update the
internal state of `ewald` so that the species `ij` is now at the given `positions`.

Note that this modifies the underlying [`EwaldContext`](@ref) `ewald.ctx`.
"""
function update_ewald_context!(ewald::IncrementalEwaldContext)
    ij = ewald.last[]
    ij ≤ 0 && return # last moved molecule bears no charge TODO: check
    ewald.sums[:,1] .+= ewald.tmpsums .- @view(ewald.sums[:,ij+1])
    ewald.sums[:,ij+1] .= ewald.tmpsums
    Eikx, Eiky, Eikz = ewald.ctx.Eiks
    newEikx, newEiky, newEikz = ewald.tmpEiks
    chargepos = ewald.ctx.indexof[ij]
    Eikx[:,chargepos] .= newEikx
    Eiky[:,chargepos] .= newEiky
    Eikz[:,chargepos] .= newEikz
    ewald.last[] = 0
    nothing
end

function add_one_system!(ewald::IncrementalEwaldContext, i::Int, positions)
    isempty(ewald.ctx.charges) && return 0
    ij = add_one_system!(ewald.ctx, i, positions)
    initialized = ewald.last[] != -1
    a, b = size(ewald.sums)
    resize!(ewald.sums, a, b+1)
    if initialized
        Eikx, Eiky, Eikz = ewald.tmpEiks
        m = length(positions)
        kx, ky, kz = ewald.ctx.eframework.kspace.ks
        kxp = kx+1; tkyp = ky+ky+1; tkzp = kz+kz+1
        resize!(Eikx, kxp, m); resize!(Eiky, tkyp, m); resize!(Eikz, tkzp, m)
        newsums = @view(ewald.sums[:,b+1])
        update_sums!(newsums, ewald, ij, (Eikx, Eiky, Eikz), false)
        @views ewald.sums[:,1] .+= newsums
        ewald.last[] = 0 # system is still initialized but tmpEiks were modified
    end
    return ij
end

function remove_one_system!(ewald::IncrementalEwaldContext, ij::Int)
    isempty(ewald.ctx.charges) && return 0
    oldij = remove_one_system!(ewald.ctx, ij)
    initialized = ewald.last[] != -1
    a, b = size(ewald.sums)
    if initialized
        ijp1 = ij + 1
        @views ewald.sums[:,1] .-= ewald.sums[:,ijp1]
        if ijp1 != b
            ewald.sums[:,ijp1] .= @view ewald.sums[:,b]
        end
        ewald.last[] = 0 # system is still initialized
    end
    resize!(ewald.sums, a, b-1)
    return oldij
end
