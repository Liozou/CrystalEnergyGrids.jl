# Interface to RASPA2

using DoubleArrayTries: DoubleArrayTrie, CommonPrefixSearch, lookup
using AtomsIO
using Unitful, UnitfulAtomic
using Printf

const RASPADIR = Ref(joinpath(homedir(), "RASPA2", "simulations", "share", "raspa"))

"""
    setdir_RASPA!(path)

Sets the `path` to the directory containing the forcefields, frameworks, grids, molecules
and structures. Default is $(RASPADIR[]). Access the value via [`getdir_RASPA`](@ref).
"""
function setdir_RASPA!(pos)
    RASPADIR[] = pos
end

"""
    getdir_RASPA()

Gets the `path` to the directory containing the forcefields, frameworks, grids, molecules
and structures. Default is $(RASPADIR[]). Set the value via [`setdir_RASPA!`](@ref).
"""
function getdir_RASPA()
    ret::String = RASPADIR[]
    isdir(ret) || error(lazy"Could not find raspa directory at the given path $ret. Please set the correct path through the `setdir_RASPA!` function.")
    ret
end

"""
    PseudoAtomInfo

Information stored in pseudo_atoms.def for a given species
"""
struct PseudoAtomInfo
    type::Symbol
    printas::Symbol
    symbol::Symbol
    oxidation::Float64
    mass::Float64
    charge::Float64
    polarization::Float64
    Bfactor::Float64
    radius::Float64
    connectivity::Float64
    anisotropic::Float64
    anisotropy_absolute::Bool
    tinker_type::Int
end

"""
    PseudoAtomListing

Internal representation of the pseudo_atoms.def file.

Can be queried with the `getindex` syntax:

```julia
julia> pal = CrystalEnergyGrids.parse_pseudoatoms_RASPA("/path/to/pseudo_atoms.def");

julia> pal[:Rb]
CrystalEnergyGrids.PseudoAtomInfo(:Rb, :Rb, :Rb, 0.0, 85.4678, 0.9094, 0.0, 1.0, 1.0, 0.0, 0.0, true, 0)
```

By default, numbers trailing after the last underscore are removed: `Oz_3` is interpreted
as `Oz`. To change this behaviour, query with `CrystalEnergyGrids.get_strict` instead of
`getindex`.
"""
struct PseudoAtomListing
    dat::Union{DoubleArrayTrie,Nothing}
    priority::Vector{Int}
    exact::Dict{String,Int}
    info::Vector{PseudoAtomInfo}
end

"""
    get_strict(pal::PseudoAtomListing, atom)

Identical to calling `pal[atom]` except that the name `atom` is kept as is, with no
heuristic transformation. Specifically, while `pal[:He_2]` will be computed as `pal[:He]`,
`CrystalEnergyGrids.get_strict(pal, :He_2)` will look for an atom exactly called `He_2`.
"""
function get_strict(pal::PseudoAtomListing, atom)
    name = atom isa String ? atom : String(atom)
    curr_priority = get(pal.exact, name, 0)
    if pal.dat isa DoubleArrayTrie
        for (id, _) in CommonPrefixSearch(pal.dat, name)
            j = pal.priority[id]
            if j > curr_priority
                curr_priority = j
            end
        end
    end
    curr_priority == 0 && throw(ArgumentError("Atom $atom not found in pseudo_atoms.def"))
    return pal.info[curr_priority]
end
function Base.getindex(pal::PseudoAtomListing, atom)
    name = get_atom_name(atom)
    get_strict(pal, name)
end

"""
    parse_pseudoatoms_RASPA(path)

Return a [`PseudoAtomListing`](@ref) extracted from the `pseudo_atoms.def` file at `path`.
"""
function parse_pseudoatoms_RASPA(file)
    lines = filter!(x -> !isempty(x) && x[1] != '#', readlines(file))
    n = parse(Int, popfirst!(lines))
    n == length(lines) || error(lazy"Found $(length(lines)) non-empty lines but $n declared")
    names = String[]
    namespos = Int[]
    exact = Dict{String,Int}()
    info = Vector{PseudoAtomInfo}(undef, n)
    for (i,_l) in enumerate(lines)
        l = split(_l)
        length(l) == 14 || error(lazy"malformed line \"$l\" does not contain 14 fields.")
        is_open_ended = l[1][end] == '_'
        name = is_open_ended ? l[1][1:end-1] : l[1]
        if is_open_ended
            push!(names, name)
            push!(namespos, i)
        else
            exact[name] = i
        end
        info[i] = PseudoAtomInfo(Symbol(name),
                                 l[2] == "yes" ? Symbol(l[3]) : Symbol(""), Symbol(l[4]),
                                 parse(Float64, l[5]), parse(Float64, l[6]),
                                 parse(Float64, l[7]), parse(Float64, l[8]),
                                 parse(Float64, l[9]), parse(Float64, l[10]),
                                 parse(Float64, l[11]), parse(Float64, l[12]),
                                 l[13] == "absolute", parse(Int, l[14]))
    end

    isempty(names) && return PseudoAtomListing(nothing, Int[], exact, info)

    dat = DoubleArrayTrie(copy(names))
    rev_priority = [lookup(dat, name) for name in names]
    priority = Vector{Int}(undef, n)
    for (i,j) in enumerate(rev_priority)
        priority[j] = namespos[i]
    end

    return PseudoAtomListing(dat, priority, exact, info)
end


"""
    RASPASystem <: AbstractSystem{3}

`AtomsBase`-compliant system extracted from a RASPA input.

Use `load_framework_RASPA`[@ref] or `load_molecule_RASPA`[@ref] to create one.
"""
struct RASPASystem <: AbstractSystem{3}
    bounding_box::SVector{3,SVector{3,TÅ}}
    position::Vector{SVector{3,TÅ}}
    atomic_symbol::Vector{Symbol}
    atomic_number::Vector{Int}
    atomic_mass::Vector{typeof(1.0u"u")}
    atomic_charge::Vector{Te_au}
    ismolecule::Bool
end
function AtomsBase.boundary_conditions(sys::RASPASystem)
    SVector{3,BoundaryCondition}(if sys.ismolecule || any(Base.Fix1(any, isinf), bounding_box(sys))
        (DirichletZero(), DirichletZero(), DirichletZero())
    else
        (Periodic(), Periodic(), Periodic())
    end)
end
AtomsBase.bounding_box(sys::RASPASystem)     = sys.bounding_box
Base.length(sys::RASPASystem)                = length(sys.position)
Base.size(sys::RASPASystem)                  = size(sys.position)
AtomsBase.species_type(::RASPASystem)        = AtomView{RASPASystem}
Base.getindex(sys::RASPASystem, i::Integer)  = AtomView(sys, i)
AtomsBase.position(s::RASPASystem)           = s.position
AtomsBase.atomic_symbol(s::RASPASystem)      = s.atomic_symbol
AtomsBase.atomic_number(s::RASPASystem)      = s.atomic_number
AtomsBase.atomic_mass(s::RASPASystem)        = s.atomic_mass
AtomsBase.velocity(::RASPASystem)            = missing
AtomsBase.atomic_symbol(sys::RASPASystem, index) = atomic_symbol(sys[index])
function Base.getindex(system::RASPASystem, x::Symbol)
    x === :bounding_box && return bounding_box(system)
    x === :boundary_conditions && return boundary_conditions(system)
    throw(KeyError("Key $x not found"))
end
Base.keys(::RASPASystem) = (:bounding_box, :boundary_conditions)
Base.haskey(::RASPASystem, x::Symbol) = (x === :bounding_box) | (x === :boundary_conditions)
AtomsBase.atomkeys(::RASPASystem) = (:position, :atomic_symbol, :atomic_number, :atomic_mass, :atomic_charge)
AtomsBase.hasatomkey(system::RASPASystem, x::Symbol)      = x in atomkeys(system)
Base.getindex(system::RASPASystem, i::Integer, x::Symbol) = getfield(system, x)[i]
Base.getindex(system::RASPASystem, ::Colon, x::Symbol)    = getfield(system, x)

"""
    load_framework_RASPA(name::AbstractString, forcefield::AbstractString)

Load the framework called `name` with the given `forcefield` from the RASPA directory.
Return a [`RASPASystem`](@ref).
"""
function load_framework_RASPA(name::AbstractString, forcefield::AbstractString)
    raspa::String = getdir_RASPA()
    cif = joinpath(raspa, "structures", "cif", splitext(name)[2] == ".cif" ? name : name*".cif")
    system = load_system(AtomsIO.ChemfilesParser(), cif)
    pseudoatoms = parse_pseudoatoms_RASPA(joinpath(raspa, "forcefield", forcefield, "pseudo_atoms.def"))
    mass  = Vector{typeof(1.0u"u")}(undef, length(system))
    charges = Vector{Te_au}(undef, length(system))
    for (i, atom) in enumerate(system)
        pseudo::PseudoAtomInfo = pseudoatoms[AtomsBase.atomic_symbol(atom)]
        mass[i]    = pseudo.mass*u"u"
        charges[i] = pseudo.charge*u"e_au"
    end
    RASPASystem(AtomsBase.bounding_box(system), AtomsBase.position(system), AtomsBase.atomic_symbol(system), AtomsBase.atomic_number(system), mass, charges, false)
end

function load_framework_RASPA(input::AbstractMatrix{T}, ::AbstractString) where T
    mat = T <: AbstractFloat ? input*u"Å" : uconvert.(u"Å", input)
    bbox = SVector{3}([[mat[1,1], mat[2,1], mat[3,1]],
                       [mat[1,2], mat[2,2], mat[3,2]],
                       [mat[1,3], mat[2,3], mat[3,3]]])
    RASPASystem(bbox, SVector{3,TÅ}[], Symbol[], Int[], typeof(1.0u"u")[], Te_au[], false)
end
