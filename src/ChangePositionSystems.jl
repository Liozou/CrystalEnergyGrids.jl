# forwards calls to the underlying system, except for the positions which are fixed
using Unitful

struct ChangePositionSystem{N,T<:AbstractSystem{N}} <: AbstractSystem{N}
    system::T
    positions::Vector{SVector{N,typeof(1.0u"Å")}}
end
ChangePositionSystem(s::AbstractSystem{N}, poss) where {N} = ChangePositionSystem{N,typeof(s)}(s, [SVector{N,typeof(1.0u"Å")}(p) for p in poss])
ChangePositionSystem(s::ChangePositionSystem, poss) = ChangePositionSystem(s.system, poss)

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
    @eval AtomsBase.$f(m::ChangePositionSystem, i) = atomsBase.$f(m.system, i)
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

struct ChangePositionAtom{N,T<:AbstractSystem{N}}
    cps::ChangePositionSystem{N,T}
    idx::Int
end
AtomsBase.position(cpa::ChangePositionAtom)      = AtomsBase.position(cpa.cps, cpa.idx)
AtomsBase.atomic_mass(cpa::ChangePositionAtom)   = AtomsBase.atomic_mass(cpa.cps, cpa.idx)
AtomsBase.atomic_symbol(cpa::ChangePositionAtom) = AtomsBase.atomic_symbol(cpa.cps, cpa.idx)
AtomsBase.atomic_number(cpa::ChangePositionAtom) = AtomsBase.atomic_number(cpa.cps, cpa.idx)
AtomsBase.n_dimensions(cpa::ChangePositionAtom)  = AtomsBase.n_dimensions(cpa.cps)
AtomsBase.element(cpa::ChangePositionAtom)       = AtomsBase.element(AtomsBase.atomic_number(cpa))

Base.show(io::IO, at::ChangePositionAtom) = AtomsBase.show_atom(io, at)
Base.show(io::IO, mime::MIME"text/plain", at::ChangePositionAtom) = AtomsBase.show_atom(io, mime, at)

Base.getindex(v::ChangePositionSystem, x::Symbol) = getindex(v.cps, v.idx, x)
Base.haskey(v::ChangePositionSystem, x::Symbol)   = hasatomkey(v.cps, x)
function Base.get(v::ChangePositionSystem, x::Symbol, default)
    AtomsBase.hasatomkey(v.cps, x) ? v[x] : default
end
Base.keys(v::ChangePositionSystem) = AtomsBase.atomkeys(v.cps)
Base.pairs(at::ChangePositionSystem) = (k => at[k] for k in keys(at))


AtomsBase.species_type(::ChangePositionSystem{N,T}) where {N,T} = ChangePositionAtom{N,T}
Base.getindex(sys::ChangePositionSystem, i::Integer)  = ChangePositionAtom(sys, i)
