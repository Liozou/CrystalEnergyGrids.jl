struct ProbeSystem{T<:AbstractSystem{3}}
    system::T
    mat::SMatrix{3,3,Float64,9}
    invmat::SMatrix{3,3,Float64,9}
    forcefield::Forcefield
    atomkinds::Vector{Int}
    probe::Int
end
function ProbeSystem(system::AbstractSystem3, forcefield::Forcefield, atom::Symbol)
    atomkinds = [get_atom_name(atomic_symbol(system, i)) for i in 1:length(system)]
    mat = stack3(bounding_box(system))
    invmat = inv(mat)
    probe = get_atom_name(atom)
    ProbeSystem(system, mat, invmat, forcefield, atomkinds, probe)
end

function compute_derivatives_vdw(s::ProbeSystem, pos::SVector{3,Float64})
    buffer, ortho, safemin = prepare_periodic_distance_computations(s.mat)
    buffer2 = MVector{3,Float64}(undef)
    value = 0.0
    ∂1 = zero(MVector{3,Float64})
    ∂2 = zero(MVector{3,Float64})
    ∂3 = 0.0
    n = length(s.framework)
    for i in 1:n
        buffer .= position(s.framework, i) .- pos
        d = periodic_distance_cartesian!(buffer, mat, invmat, ortho, safemin, buffer2)
        d ≥ s.forcefield.cutoff && continue
        v, p1, p2, p3 = derivatives_nocutoff(s.forcefield, s.atomkinds[i], s.probe, d)
        value += v
        ∂1 .+= p1 .* buffer
        d13 = buffer[1]*buffer[3]
        ∂2 .+= p2 .* SVector{3,Float64}((buffer[1]*buffer[2], d13, d13))
        ∂3 += p3 * d13 * buffer[2]
    end
    (value, ∂1, ∂2, ∂3)
end