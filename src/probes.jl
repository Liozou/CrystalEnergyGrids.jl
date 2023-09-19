"""
    ProbeSystem(framework::AbstractSystem{3}, forcefield::ForceField, [atom::Symbol])

System composed of a `framework` and an optional probe `atom` with a given `forcefield`.
Such a system is used to compute the VdW grid of the atom, or a Coulomb grid if no atom is
provided.

!!! note
    If the system has perpendicular lengths lower than 12.0 Å, the smallest supercell above
    that threshold will be used instead.
"""
struct ProbeSystem
    positions::Vector{SVector{3,Float64}} # cartesian coordinates
    mat::SMatrix{3,3,Float64,9}
    invmat::SMatrix{3,3,Float64,9}
    forcefield::ForceField
    atomkinds::Vector{Int}
    charges::Vector{Float64} # only used for Coulomb grid (otherwise empty)
    probe::Int # only used for VdW grid (otherwise 0)
end
function ProbeSystem(framework::AbstractSystem{3}, forcefield::ForceField, atom::Symbol=Symbol(""))
    n = length(framework)
    _atomkinds = [forcefield.sdict[Symbol(get_atom_name(atomic_symbol(framework, i)))] for i in 1:n]
    nx, ny, nz = find_supercell(framework, 12.0u"Å")
    numsupercell = nx*ny*nz
    vx, vy, vz = bounding_box(framework)
    mat = stack3((nx .* vx, ny .* vy, nz .* vz))
    invmat = inv(mat)
    atomkinds, positions = if numsupercell == 1
        _atomkinds, [NoUnits.(position(framework, i) ./ u"Å") for i in 1:n]
    else
        wx = NoUnits.(vx/u"Å"); wy = NoUnits.(vy/u"Å"); wz = NoUnits.(vz/u"Å")
        newpositions = Vector{SVector{3,Float64}}(undef, numsupercell*n)
        for i in 1:n
            newpositions[i] = NoUnits.(position(framework, i) ./ u"Å") # only call this once
        end
        for iz in 0:(nz-1)
            izn = iz*n
            nnz = n*nz
            stepz = iz*wz
            for iy in 0:(ny-1)
                nnynz = ny*nnz
                iyz = iy*nnz + izn
                stepyz = iy*wy + stepz
                for ix in (iz==iy==0):(nx-1)
                    ixyz = ix*nnynz + iyz
                    stepxyz = ix*wx + stepyz
                    for i in 1:n
                        newpositions[ixyz + i] = newpositions[i] + stepxyz
                    end
                end
            end
        end
        repeat(_atomkinds, numsupercell), newpositions
    end
    if atom === Symbol("")
        charges = repeat([NoUnits(x[:atomic_charge] / u"e_au") for x in framework], numsupercell)
        ProbeSystem(positions, mat, invmat, forcefield, atomkinds, charges, 0)
    else
        probe = forcefield.sdict[Symbol(get_atom_name(atom))]
        ProbeSystem(positions, mat, invmat, forcefield, atomkinds, Float64[], probe)
    end
end

function Base.show(io::IO, ::MIME"text/plain", x::ProbeSystem)
    print(io, "ProbeSystem with atom ")
    print(io, findfirst(==(x.probe), x.forcefield.sdict))
    nothing
end

function compute_derivatives_vdw(s::ProbeSystem, pos::SVector{3,typeof(1.0u"Å")})
    buffer, ortho, safemin = prepare_periodic_distance_computations(s.mat)
    safemin2 = safemin^2
    buffer2 = MVector{3,Float64}(undef)
    cutoff2 = NoUnits(s.forcefield.cutoff^2/u"Å^2")
    value = 0.0
    ∂1 = zero(MVector{3,Float64})
    ∂2 = zero(MVector{3,Float64})
    ∂3 = 0.0
    for i in 1:length(s.positions)
        buffer .= NoUnits.(pos./u"Å") .- s.positions[i]
        d2 = periodic_distance2_fromcartesian!(buffer, s.mat, s.invmat, ortho, safemin2, buffer2)
        d2 ≥ cutoff2 && continue
        v, p1, p2, p3 = derivatives_nocutoff(s.forcefield, s.atomkinds[i], s.probe, d2)
        value += v
        ∂1 .+= p1 .* buffer
        d13 = buffer[1]*buffer[3]
        ∂2 .+= p2 .* SVector{3,Float64}((buffer[1]*buffer[2], d13, buffer[2]*buffer[3]))
        ∂3 += p3 * d13 * buffer[2]
    end
    (value, SVector{3,Float64}(∂1), SVector{3,Float64}(∂2), ∂3)
end

function compute_derivatives_ewald(s::ProbeSystem, ewald, pos::SVector{3,typeof(1.0u"Å")})
    buffer, ortho, safemin = prepare_periodic_distance_computations(s.mat)
    safemin2 = safemin^2
    buffer2 = MVector{3,Float64}(undef)
    cutoff2 = NoUnits(s.forcefield.cutoff^2/u"Å^2")
    value = 0.0
    ∂1 = zero(MVector{3,Float64})
    ∂2 = zero(MVector{3,Float64})
    ∂3 = 0.0
    smallest_d2 = Inf
    for i in 1:length(s.positions)
        buffer .= NoUnits.(pos./u"Å") .- s.positions[i]
        d2 = periodic_distance2_fromcartesian!(buffer, s.mat, s.invmat, ortho, safemin2, buffer2)
        d2 ≥ cutoff2 && continue
        smallest_d2 = min(smallest_d2, d2)
        v, p1, p2, p3 = derivatives_ewald(ewald, s.charges[i], d2)
        value += v
        ∂1 .+= p1 .* buffer
        d13 = buffer[1]*buffer[3]
        ∂2 .+= p2 .* SVector{3,Float64}((buffer[1]*buffer[2], d13, buffer[2]*buffer[3]))
        ∂3 += p3 * d13 * buffer[2]
    end
    ((smallest_d2 < 1.0 ? Inf : value), SVector{3,Float64}(∂1), SVector{3,Float64}(∂2), ∂3)
end
