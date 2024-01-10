# forwards calls to the underlying system, except for the positions which are fixed
using Unitful

"""
    ChangePositionSystem{N,T<:AbstractSystem{N}} <: AbstractSystem{N}

`AbstractSystem` consisting of a small layer above a reference `AbstractSystem` with
modified positions.

See [ChangePositionSystem(s::AbstractSystem{N}, poss)](@ref) and
[ChangePositionSystem(s::AbstractSystem{N}; translate, rotate, bead)](@ref)
to construct an instance.
"""
struct ChangePositionSystem{N,T<:AbstractSystem{N}} <: AbstractSystem{N}
    system::T
    positions::Vector{SVector{N,TÅ}}
end

"""
    ChangePositionSystem(s::AbstractSystem{N}, poss) where N

Construct a [`ChangePositionSystem`](@ref) consisting of the reference system `s` with
atomic positions given by `poss`.

`poss` must be an iterator over atomic coordinates (with their unit), such as a
`Vector{SVector{N, typeof(1.0u"Å")}}` for instance.

See also [ChangePositionSystem(s::AbstractSystem{N}; translate, rotate, bead)](@ref).
"""
ChangePositionSystem(s::AbstractSystem{N}, poss) where {N} = ChangePositionSystem{N,typeof(s)}(s, [SVector{N,TÅ}(p) for p in poss])
ChangePositionSystem(s::ChangePositionSystem, poss) = ChangePositionSystem(s.system, poss)
function ChangePositionSystem(s::ChangePositionSystem, poss::Vector{SVector{N,TÅ}}) where N
    ChangePositionSystem(s.system, poss) # to avoid ambiguity
end

"""
ChangePositionSystem(s::AbstractSystem{N}; translate=nothing, rotate=nothing, bead=0) where N

Construct a [`ChangePositionSystem`](@ref) consisting of the reference system `s` after
rotation and translation.

If `rotate` is set to a rotation matrix, each atom will have its position multiplied by it.
If `bead` is set to the number of an atom of the system, the rotation will occur relative
to that atom (i.e., that atom will not be moved).

If `translate` is set to a position, each atom will have its position translated by that
amount.

If `bead == 0`, the translation (if any) occurs after the rotation (if any). Otherwise, the
order does not matter since the rotation is relative to the position of an atom that will
be translated with the others.

See also [ChangePositionSystem(s::AbstractSystem{N}, poss)](@ref).
"""
function ChangePositionSystem(s::AbstractSystem{N}; translate=nothing, rotate=nothing, bead=0) where N
    n = length(s)
    poss = Vector{SVector{N,TÅ}}(undef, n)
    _t = translate isa Nothing ? zero(SVector{N,TÅ}) : SVector{N,TÅ}(translate)
    for i in 1:n
        p = position(s, i)::SVector{N,TÅ}
        if !(rotate isa Nothing)
            if bead == 0
                p = SVector{N,TÅ}(rotate * p)
            else
                p = p[bead] .+ SVector{N,TÅ}(rotate * (p .- p[bead]))
            end
        end
        if !(translate isa Nothing)
            p = p .+ _t
        end
        poss[i] = p
    end
    ChangePositionSystem(s, poss)
end

for f in (:boundary_conditions,
          :bounding_box,
          :chemical_formula,
          :isinfinite,
          :n_dimensions,
          :periodicity,
          :atomkeys,
          :length,
          :size,
         )
    @eval AtomsBase.$f(m::ChangePositionSystem) = AtomsBase.$f(m.system)
end
for f in (:atomic_mass,
          :atomic_number,
          :atomic_symbol,
          :velocity,
         )
    @eval AtomsBase.$f(m::ChangePositionSystem) = AtomsBase.$f(m.system)
    @eval AtomsBase.$f(m::ChangePositionSystem, i) = AtomsBase.$f(m.system, i)
end
AtomsBase.hasatomkey(m::ChangePositionSystem, x::Symbol) = AtomsBase.hasatomkey(m.system, x)

AtomsBase.position(m::ChangePositionSystem) = m.positions
AtomsBase.position(m::ChangePositionSystem, i) = m.positions[i]

function Base.getindex(sys::ChangePositionSystem, i, s::Symbol)
    if s === :position
        sys.positions[i]
    else
        sys.system[i,s]
    end
end

Base.haskey(v::ChangePositionSystem, x::Symbol) = x === :position || haskey(v.system, x)
function Base.get(v::ChangePositionSystem, x::Symbol, default)
    x === :position ? v.positions : get(v.system, x, default)
end
function Base.getindex(v::ChangePositionSystem, x::Symbol)
    x === :position ? v.positions : getindex(v.system, x)
end
function Base.keys(v::ChangePositionSystem)
    if haskey(v.system, :position)
        keys(v.system)
    else
        (:position, keys(v.system)...)
    end
end


struct ChangePositionAtom{N,T<:AbstractSystem{N}}
    cps::ChangePositionSystem{N,T}
    idx::Int
end
AtomsBase.position(cpa::ChangePositionAtom)      = AtomsBase.position(cpa.cps, cpa.idx)

for f in (:atomic_mass,
    :atomic_number,
    :atomic_symbol,
    :velocity,
)
@eval AtomsBase.$f(m::ChangePositionAtom) = AtomsBase.$f(m.cps, m.idx)
end
AtomsBase.n_dimensions(cpa::ChangePositionAtom)  = AtomsBase.n_dimensions(cpa.cps)
AtomsBase.element(cpa::ChangePositionAtom)       = AtomsBase.element(AtomsBase.atomic_number(cpa))

Base.show(io::IO, at::ChangePositionAtom) = AtomsBase.show_atom(io, at)
Base.show(io::IO, mime::MIME"text/plain", at::ChangePositionAtom) = AtomsBase.show_atom(io, mime, at)

Base.getindex(v::ChangePositionAtom, x::Symbol) = getindex(v.cps, v.idx, x)
Base.haskey(v::ChangePositionAtom, x::Symbol)   = hasatomkey(v.cps, x)
function Base.get(v::ChangePositionAtom, x::Symbol, default)
    AtomsBase.hasatomkey(v.cps, x) ? v[x] : default
end
Base.keys(v::ChangePositionAtom) = AtomsBase.atomkeys(v.cps)
Base.pairs(at::ChangePositionAtom) = (k => at[k] for k in keys(at))


AtomsBase.species_type(::ChangePositionSystem{N,T}) where {N,T} = ChangePositionAtom{N,T}
Base.getindex(sys::ChangePositionSystem, i::Integer)  = ChangePositionAtom(sys, i)
