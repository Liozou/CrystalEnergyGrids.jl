# interface to RASPA2

using DoubleArrayTries, AtomsIO
using Printf

export set_raspadir!, get_raspadir
export PseudoAtomListing, parse_pseudoatoms
export RaspaSystem, raspa_setup

const RASPADIR = Ref(joinpath(ENV["HOME"], "RASPA2", "simulations", "share", "raspa"))

"""
    set_raspadir!(path)

Sets the `path` to the directory containing the forcefields, frameworks, grids, molecules
and structures. Default is $(RASPADIR[]). Access the value via `get_raspadir`[@ref].
"""
function set_raspadir!(pos)
    RASPADIR[] = pos
end
"""
    get_raspadir()

Gets the `path` to the directory containing the forcefields, frameworks, grids, molecules
and structures. Default is $(RASPADIR[]). Set the value via `set_raspadir!`[@ref].
"""
get_raspadir() = RASPADIR[]

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
julia> pal = CrystalEnergyGrids.parse_pseudoatoms("/path/to/pseudo_atoms.def");

julia> pal[:Rb]
CrystalEnergyGrids.PseudoAtomInfo(:Rb, :Rb, :Rb, 0.0, 85.4678, 0.9094, 0.0, 1.0, 1.0, 0.0, 0.0, true, 0)
```

By default, numbers trailing after the last underscore are removed: `Oz_3` is interpreted
as `Oz`. To chage this behaviour, query with `CrystalEnergyGrids.get_strict` instead of
`getindex`.
"""
struct PseudoAtomListing
    dat::Union{DoubleArrayTrie,Nothing}
    priority::Vector{Int}
    exact::Dict{String,Int}
    info::Vector{PseudoAtomInfo}
end

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
    name = atom isa String ? atom : String(atom)
    s = split(name, '_')
    if length(s) > 1 && all(isnumeric, last(s))
        name = join(@view(s[1:end-1]), '_')
    end
    get_strict(pal, name)
end

"""
    parse_pseudoatoms(path)

Return a `PseudoAtomListing` extracted from the `pseudo_atoms.def` file at `path`.
"""
function parse_pseudoatoms(file)
    lines = filter!(x -> !isempty(x) && x[1] != '#', readlines(file))
    n = parse(Int, popfirst!(lines))
    n == length(lines) || error(lazy"Found $(length(lines)) non-empty lines but $n declared")
    names = String[]
    namespos = Int[]
    exact = Dict{String,Int}()
    info = Vector{PseudoAtomInfo}(undef, n)
    for (i,_l) in enumerate(lines)
        l = split(_l)
        length(l) == 14 || error("malformed line \"$l\" does not contain 14 fields.")
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

function parse_molecule(file)
    lines = filter!(x -> !isempty(x) && x[1] != '#', readlines(file))
    num_atoms = parse(Int, lines[4])
    num_groups = parse(Int, lines[5])
    @assert num_groups == 1
    @assert num_atoms == parse(Int, lines[7])
    @assert num_atoms == 1 || strip(lines[6]) == "rigid"
    positions = Vector{SVector{3,typeof(ANG_UNIT)}}(undef, num_atoms)
    symbols = Vector{Symbol}(undef, num_atoms)
    for i in 1:num_atoms
        l = split(lines[7+i])
        symbols[i] = Symbol(l[2])
        positions[i] = if num_atoms == 1
            zero(SVector{3,Float64})*ANG_UNIT
        else
            SVector{3,Float64}(parse(Float64, l[3]), parse(Float64, l[4]), parse(Float64, l[5]))*ANG_UNIT
        end
    end
    (symbols, positions)
end

"""
    RaspaSystem <: AbstractSystem{3}

`AtomsBase`-compliant system extracted from a RASPA input.

Use `load_RASPA_framework`[@ref] or `load_RASPA_molecule`[@ref] to create one.
"""
struct RaspaSystem <: AbstractSystem{3}
    bounding_box::SVector{3,SVector{3,typeof(ANG_UNIT)}}
    position::Vector{SVector{3,typeof(ANG_UNIT)}}
    atomic_symbol::Vector{Symbol}
    atomic_number::Vector{Int}
    atomic_mass::Vector{typeof(ATOMMASS_UNIT)}
    atomic_charge::Vector{typeof(CHARGE_UNIT)}
end
AtomsBase.bounding_box(sys::RaspaSystem)     = sys.bounding_box
AtomsBase.boundary_conditions(sys::RaspaSystem) = SVector{3,BoundaryCondition}(ntuple(Returns(ifelse(isinf(bounding_box(sys)[1][1]), DirichletZero(), Periodic())), 3))
Base.length(sys::RaspaSystem)                = length(sys.position)
Base.size(sys::RaspaSystem)                  = size(sys.position)
AtomsBase.species_type(::RaspaSystem)        = AtomView{RaspaSystem}
Base.getindex(sys::RaspaSystem, i::Integer)  = AtomView(sys, i)
AtomsBase.position(s::RaspaSystem)           = s.position
AtomsBase.atomic_symbol(s::RaspaSystem)      = s.atomic_symbol
AtomsBase.atomic_number(s::RaspaSystem)      = s.atomic_number
AtomsBase.atomic_mass(s::RaspaSystem)        = s.atomic_mass
AtomsBase.velocity(::RaspaSystem)            = missing
AtomsBase.atomic_symbol(sys::RaspaSystem, index) = atomic_symbol(sys[index])
function Base.getindex(system::RaspaSystem, x::Symbol)
    x === :bounding_box && return bounding_box(system)
    x === :boundary_conditions && return boundary_conditions(system)
    throw(KeyError("Key $x not found"))
end
Base.keys(::RaspaSystem)          = (:bounding_box, :boundary_conditions)
AtomsBase.atomkeys(::RaspaSystem) = (:position, :atomic_symbol, :atomic_number, :atomic_mass, :atomic_charge)
AtomsBase.hasatomkey(system::RaspaSystem, x::Symbol)      = x in atomkeys(system)
Base.getindex(system::RaspaSystem, i::Integer, x::Symbol) = getfield(system, x)[i]
Base.getindex(system::RaspaSystem, ::Colon, x::Symbol)    = getfield(system, x)
atomic_charge(s::RaspaSystem) = s.atomic_charge

"""
    load_RASPA_framework(name, forcefield)

Load the framework called `name` with the given `forcefield` from the RASPA directory.
Return a [`RaspaSystem`](@ref).
"""
function load_RASPA_framework(name, forcefield)
    raspa::String = get_raspadir()
    cif = joinpath(raspa, "structures", "cif", splitext(name)[2] == ".cif" ? name : name*".cif")
    system = load_system(AtomsIO.ChemfilesParser(), cif)
    pseudoatoms = parse_pseudoatoms(joinpath(raspa, "forcefield", forcefield, "pseudo_atoms.def"))
    mass  = Vector{typeof(ATOMMASS_UNIT)}(undef, length(system))
    charges = Vector{typeof(CHARGE_UNIT)}(undef, length(system))
    for (i, atom) in enumerate(system)
        pseudo::PseudoAtomInfo = pseudoatoms[atomic_symbol(atom)]
        mass[i]    = pseudo.mass*ATOMMASS_UNIT
        charges[i] = pseudo.charge*CHARGE_UNIT
    end
    RaspaSystem(bounding_box(system), position(system), atomic_symbol(system), atomic_number(system), mass, charges)
end

"""
    load_RASPA_molecule(name::AbstractString, forcefield::AbstractString, framework_forcefield::AbstractString, framework_system::Union{AbstractSystem,Nothing}=nothing)

Load the molecule called `name` with the given `forcefield` from the RASPA directory.
The forcefield of the associated framework must be given in order to correctly set the
atomic numbers, masses and charges of the constituents of the molecule.

The framework system may optionnally be passed to set the corresponding bounding box.
Otherwise, an infinite bounding box is assumed as well as no periodic conditions.

Return a [`RaspaSystem`](@ref).
"""
function load_RASPA_molecule(name, forcefield, framework_forcefield, framework_system=nothing)
    raspa::String = get_raspadir()
    symbols, positions = parse_molecule(joinpath(raspa, "molecules", forcefield, splitext(name)[2] == ".def" ? name : name*".def"))
    pseudoatoms = parse_pseudoatoms(joinpath(raspa, "forcefield", framework_forcefield, "pseudo_atoms.def"))
    mass  = Vector{typeof(ATOMMASS_UNIT)}(undef, length(symbols))
    charges = Vector{typeof(CHARGE_UNIT)}(undef, length(symbols))
    atom_numbers = Vector{Int}(undef, length(symbols))
    for (i, atom) in enumerate(symbols)
        pseudo::PseudoAtomInfo = pseudoatoms[atom]
        mass[i]    = pseudo.mass*ATOMMASS_UNIT
        charges[i] = pseudo.charge*CHARGE_UNIT
        atom_numbers[i] = AtomsBase.element(pseudo.symbol).number
    end
    bbox = if framework_system !== nothing
        bounding_box(framework_system)
    else
        SVector{3,SVector{3,typeof(ANG_UNIT)}}([[Inf, 0.0, 0.0]*ANG_UNIT, [0.0, Inf, 0.0]*ANG_UNIT, [0.0, 0.0, Inf]*ANG_UNIT])
    end
    RaspaSystem(bbox, positions, symbols, atom_numbers, mass, charges)
end


"""
    raspa_setup(framework, forcefield_framework, molecule, forcefield_molecule, gridstep=0.15, supercell=(1,1,1))

Return a [`CrystalEnergySetup`](@ref) for studying `molecule` (with its forcefield) in
`framework` (with its forcefield), extracted from existing RASPA grids and completed with
Ewald sums precomputations.

`gridsteps` and `supercell` should match that used when creating the grids.
"""
function raspa_setup(framework, forcefield_framework, molecule, forcefield_molecule, gridstep=0.15, supercell=(1,1,1))
    raspa::String = get_raspadir()
    syst_framework = load_RASPA_framework(framework, forcefield_framework)
    gridstep_name = @sprintf "%.6f" gridstep
    supercell_name = join(supercell, 'x')
    grid_dir = joinpath(raspa, "grids", forcefield_framework, framework, gridstep_name)
    coulomb = parse_grid(joinpath(grid_dir, supercell_name, framework*"_Electrostatics_Ewald.grid"), true)

    syst_mol = load_RASPA_molecule(molecule, forcefield_molecule, forcefield_framework, syst_framework)
    trunc_or_shift = strip(first(Iterators.filter(!startswith('#'), eachline(joinpath(raspa, "forcefield", forcefield_framework, "force_field_mixing_rules.def")))))
    atomdict = IdDict{Symbol,Int}()
    atoms = syst_mol[:,:atomic_symbol]
    for atom in atoms
        get!(Returns(length(atomdict)+1), atomdict, atom)
    end
    atomsidx = [atomdict[atom] for atom in atoms]
    grids = Vector{CrystalEnergyGrid}(undef, length(atomdict))
    for (atom, i) in atomdict
        grids[i] = parse_grid(joinpath(grid_dir, string(framework, '_', atom, '_', trunc_or_shift, ".grid")), false)
    end

    ewald = initialize_ewald(syst_framework, supercell)

    return CrystalEnergySetup(syst_framework, syst_mol, coulomb, grids, atomsidx, ewald)
end
