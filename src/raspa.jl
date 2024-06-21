# Interface to RASPA2

using DoubleArrayTries: DoubleArrayTrie, CommonPrefixSearch, lookup
using AtomsIO
using Unitful, UnitfulAtomic
using Printf

export setdir_RASPA!, getdir_RASPA
export PseudoAtomListing, parse_pseudoatoms_RASPA, parse_forcefield_RASPA
export RASPASystem, setup_RASPA, setup_probe_RASPA

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

function PseudoAtomInfo(pa::PseudoAtomInfo;
                        type::Symbol=pa.type,
                        printas::Symbol=pa.printas,
                        symbol::Symbol=pa.symbol,
                        oxidation::Float64=pa.oxidation,
                        mass::Float64=pa.mass,
                        charge::Float64=pa.charge,
                        polarization::Float64=pa.polarization,
                        Bfactor::Float64=pa.Bfactor,
                        radius::Float64=pa.radius,
                        connectivity::Float64=pa.connectivity,
                        anisotropic::Float64=pa.anisotropic,
                        anisotropy_absolute::Bool=pa.anisotropy_absolute,
                        tinker_type::Int=pa.tinker_type)
    PseudoAtomInfo(type, printas, symbol, oxidation, mass, charge, polarization,Bfactor, radius, connectivity, anisotropic, anisotropy_absolute, tinker_type)
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
    parse_pseudoatoms_RASPA(forcefield, ::Nothing)

Return a [`PseudoAtomListing`](@ref) extracted from the `pseudo_atoms.def` file of the
given force field.
"""
function parse_pseudoatoms_RASPA(forcefield, ::Nothing)
    if forcefield isa AbstractString || forcefield isa ForceField || !(forcefield[2] isa PseudoAtomListing)
        ff = forcefield isa Tuple ? forcefield[1] : forcefield
        ffname = ff isa ForceField ? ff.name : ff
        parse_pseudoatoms_RASPA(joinpath(getdir_RASPA(), "forcefield", ffname, "pseudo_atoms.def"))
    else
        forcefield[2]::PseudoAtomListing
    end
end

function parse_molecule_RASPA(file)
    lines = filter!(x -> !isempty(x) && x[1] != '#', readlines(file))
    num_atoms = parse(Int, lines[4])
    num_groups = parse(Int, lines[5])
    @assert num_groups == 1
    @assert num_atoms == parse(Int, lines[7])
    @assert num_atoms == 1 || strip(lines[6]) == "rigid"
    positions = Vector{SVector{3,TÅ}}(undef, num_atoms)
    symbols = Vector{Symbol}(undef, num_atoms)
    for i in 1:num_atoms
        l = split(lines[7+i])
        symbols[i] = Symbol(l[2])
        positions[i] = if num_atoms == 1
            zero(SVector{3,Float64})*u"Å"
        else
            SVector{3,Float64}(parse(Float64, l[3]), parse(Float64, l[4]), parse(Float64, l[5]))*u"Å"
        end
    end
    (symbols, positions)
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
function load_framework_RASPA(name::AbstractString, forcefield)
    raspa::String = getdir_RASPA()
    cif = joinpath(raspa, "structures", "cif", splitext(name)[2] == ".cif" ? name : name*".cif")
    system = load_system(AtomsIO.ChemfilesParser(), cif)
    pseudoatoms = parse_pseudoatoms_RASPA(forcefield, nothing)
    mass  = Vector{typeof(1.0u"u")}(undef, length(system))
    charges = Vector{Te_au}(undef, length(system))
    for (i, atom) in enumerate(system)
        pseudo::PseudoAtomInfo = pseudoatoms[AtomsBase.atomic_symbol(atom)]
        mass[i]    = pseudo.mass*u"u"
        if forcefield isa AbstractString || forcefield[2] isa PseudoAtomListing
            charges[i] = pseudo.charge*u"e_au"
        else
            ff = forcefield[2]::ForceField
            interactions = ff[pseudo.type,pseudo.type]
            found_charges = false
            for inter in (interactions isa InteractionRule ? (interactions,) : interactions.rules)
                if inter.kind === FF.CoulombEwaldDirect
                    charges[i] = inter.params[2]*u"e_au"
                    found_charges = true
                    break
                end
            end
            if !found_charges
                charges[i] = pseudo.charge*u"e_au"
            end
        end
    end
    RASPASystem(AtomsBase.bounding_box(system), AtomsBase.position(system), AtomsBase.atomic_symbol(system), AtomsBase.atomic_number(system), mass, charges, false)
end

function load_framework_RASPA(input::AbstractMatrix{T}, ::Any) where T
    mat = T <: AbstractFloat ? input*u"Å" : uconvert.(u"Å", input)
    bbox = SVector{3}([[mat[1,1], mat[2,1], mat[3,1]],
                       [mat[1,2], mat[2,2], mat[3,2]],
                       [mat[1,3], mat[2,3], mat[3,3]]])
    RASPASystem(bbox, SVector{3,TÅ}[], Symbol[], Int[], typeof(1.0u"u")[], Te_au[], false)
end

"""
    load_molecule_RASPA(name::AbstractString, forcefield::AbstractString, framework_forcefield, framework_system::Union{AbstractSystem,Nothing}=nothing)

Load the molecule called `name` with the given `forcefield` from the RASPA directory.
The forcefield of the associated framework must be given in order to correctly set the
atomic numbers, masses and charges of the constituents of the molecule.

The framework system may optionnally be passed to set the corresponding bounding box.
Otherwise, an infinite bounding box is assumed as well as no periodic conditions.

Return a [`RASPASystem`](@ref).
"""
function load_molecule_RASPA(name::AbstractString, forcefield::AbstractString, framework_forcefield, framework_system::Union{AbstractSystem,Nothing}=nothing)
    raspa::String = getdir_RASPA()
    symbols, positions = parse_molecule_RASPA(joinpath(raspa, "molecules", forcefield, splitext(name)[2] == ".def" ? name : name*".def"))
    pseudoatoms = parse_pseudoatoms_RASPA(framework_forcefield, nothing)
    mass  = Vector{typeof(1.0u"u")}(undef, length(symbols))
    charges = Vector{Te_au}(undef, length(symbols))
    atom_numbers = Vector{Int}(undef, length(symbols))
    for (i, atom) in enumerate(symbols)
        pseudo::PseudoAtomInfo = pseudoatoms[atom]
        mass[i]    = pseudo.mass*1.0u"u"
        charges[i] = pseudo.charge*1.0u"e_au"
        el = get(AtomsBase.PeriodicTable.elements, pseudo.symbol, missing)
        atom_numbers[i] = ismissing(el) ? 0 : el.number
    end
    bbox = if framework_system !== nothing
        bounding_box(framework_system)
    else
        SVector{3,SVector{3,TÅ}}([[Inf, 0.0, 0.0]*u"Å", [0.0, Inf, 0.0]*u"Å", [0.0, 0.0, Inf]*u"Å"])
    end
    RASPASystem(bbox, positions, symbols, atom_numbers, mass, charges, true)
end

function default_system(s, ff)
    forcefield = parse_forcefield_RASPA(ff)
    s isa AbstractSystem && return s
    symbol = Symbol(get_atom_name(s))
    idx = get(forcefield.sdict, symbol, 0)
    idx == 0 && return load_molecule_RASPA(s, "Default", forcefield)
    pseudo = parse_pseudoatoms_RASPA(ff, nothing)[symbol]
    mass  = pseudo.mass*1.0u"u"
    charge = pseudo.charge*1.0u"e_au"
    el = get(AtomsBase.PeriodicTable.elements, pseudo.symbol, missing)
    atom_number = ismissing(el) ? 0 : el.number
    bbox = SA[SA[Inf, 0.0, 0.0]u"Å", SA[0.0, Inf, 0.0]u"Å", SA[0.0, 0.0, Inf]u"Å"]
    RASPASystem(bbox, [zero(SVector{3,Float64})*u"Å"], [symbol], [atom_number], [mass], [charge], true)
end

function decide_parse_block(blockfile, molecule, framework_name)::String
    framework_name isa AbstractMatrix && return ""
    newblockfile = if blockfile isa Union{Bool,AbstractString}
        blockfile
    else
        !(length(molecule) == 1 && molecule[1,:atomic_charge] > eps()*u"e_au")
    end
    if newblockfile isa Bool
        newblockfile || return ""
        return String(first(split(framework_name, '_')))
    end
    return String(newblockfile)
end

function parse_block(blockfile, framework_name, framework, molecule, spacing)
    blockpath = decide_parse_block(blockfile, molecule, framework_name)
    csetup = GridCoordinatesSetup(framework, spacing)
    isempty(blockpath) && return BlockFile(csetup)
    block = parse_blockfile(joinpath(getdir_RASPA(), "structures", "block", blockpath)*".block", csetup)
    return block
end


function grid_locations(framework, forcefield_framework, atoms, gridstep, supercell)
    framework isa AbstractMatrix && return "", String[]
    rawname, name = if forcefield_framework isa AbstractString
        forcefield_framework, forcefield_framework
    else
        forcefield_name, forcefield = forcefield_framework
        forcefield_name, string(forcefield_name, '-', hash_string(forcefield))
    end
    raspa::String = getdir_RASPA()
    rootdir = joinpath(raspa, "grids", name, framework)
    gridstep_name::String = @sprintf "%.6f" NoUnits(gridstep/u"Å")
    grid_dir = joinpath(rootdir, gridstep_name)
    supercell_name = join(supercell, 'x')
    coulomb_grid_path = joinpath(grid_dir, supercell_name, framework*"_Electrostatics_Ewald.grid")
    trunc_or_shift = strip(first(Iterators.filter(!startswith('#'), eachline(joinpath(raspa, "forcefield", rawname, "force_field_mixing_rules.def")))))
    vdws = Vector{String}(undef, length(atoms))
    for (i, atom) in enumerate(atoms)
        vdws[i] = joinpath(grid_dir, string(framework, '_', atom, '_', trunc_or_shift, ".grid"))
    end
    return coulomb_grid_path, vdws
end

function retrieve_or_create_grid(grid_path, syst_framework, forcefield::ForceField, gridstep, atom_or_eframework, mat, new, cutoff)
    (isempty(grid_path) || isinf(cutoff)) && return EnergyGrid(false)
    cutoff == 12.0u"Å" || error("Cutoff other than 12 Å or infinity is not supported.")
    iscoulomb = atom_or_eframework isa EwaldFramework
    !iscoulomb && !needsvdwgrid(forcefield, atom_or_eframework) && return EnergyGrid(true)
    text = iscoulomb ? "Coulomb grid" : "VdW grid for $atom_or_eframework"
    if new || !isfile(grid_path)
        mkpath(dirname(grid_path))
        printstyled("Creating ", text, " at ", grid_path, " ... "; color=:yellow)
        if iscoulomb
            create_grid_coulomb(grid_path, syst_framework, forcefield, gridstep, atom_or_eframework)
        else
            create_grid_vdw(grid_path, syst_framework, forcefield, gridstep, atom_or_eframework)
        end
        printstyled(text, " created.\n"; color=:cyan)
    else
        printstyled("Retrieved ", text, " at ", grid_path, ".\n"; color=:cyan)
    end
    parse_grid(grid_path, iscoulomb, mat)
end

"""
    setup_RASPA(framework, forcefield_framework, molecule, forcefield_molecule; gridstep=0.15u"Å", supercell=nothing, blockfile=nothing, new=false, cutoff=12.0u"Å")
    setup_RASPA(framework, forcefield_framework, syst_mol; gridstep=0.15u"Å", supercell=nothing, blockfile=nothing, new=false, cutoff=12.0u"Å")
    setup_RASPA(framework, forcefield, atom::Union{AbstractString,Symbol}; gridstep=0.15u"Å", supercell=nothing, blockfile=nothing, new=false, cutoff=12.0u"Å")

Return a [`CrystalEnergySetup`](@ref) for studying `molecule` (with its forcefield) in
`framework` (with its forcefield), extracted from existing RASPA grids and completed with
Ewald sums precomputations.
The molecule can alternatively be directly as an abstract system `syst_mol`, or, if it is
a single atom defined in the `forcefield`, as that very `atom`.

`forcefield_framework` can be the name of the forcefield, or a pair `(name, ff)` where `ff`
is a `ForceField`, in which case the given `ff` will be used. It can also be a pair
`(name, pal)` where `pal` is a `PseudoAtomListing`, in which case the forcefield will be
read based on its name, but the atom definitions will be taken from the given `pal`.

`gridstep` and `supercell` should match that used when creating the grids. If provided,
`gridstep` should be a floating point number and `supercell` a triplet of integers.
Otherwise, the default grid step is 0.15 Å and the default supercell is the smallest one
such that its perpendicular lengths are all above twice the `cutoff` (i.e. twice 12.0 Å by
default). The specified cutoff will be integrated in the definition of the force field.

`blockfile` can be set to `false` to allow the molecule to go everywhere in the framework.
Or it can be set to the radical (excluding the ".block" extension) of the block file in the
raspa directory to include it and prevent the molecule to go in the blocked spheres.
Setting it to `true` is equivalent to using `blockfile=first(split(framework, '_'))`.
If not provided, the default is to set it to `false` if the molecule is positively charged
and monoatomic (considered to be a small cation), or to `true` otherwise.

If `new` is set, new atomic grids will be computed. Otherwise, existing grids will be used
if they exist.
"""
function setup_RASPA(framework, forcefield_framework, syst_mol; gridstep=0.15u"Å", supercell=nothing, blockfile=nothing, new=false, cutoff=12.0u"Å")
    syst_framework = load_framework_RASPA(framework, forcefield_framework)
    if supercell isa Nothing
        supercell = isinf(cutoff) ? (1,1,1) : find_supercell(syst_framework, cutoff)
    end
    supercell::NTuple{3,Int}

    mat = stack3(bounding_box(syst_framework))

    forcefield = parse_forcefield_RASPA(forcefield_framework; cutoff)
    block = parse_block(blockfile, framework, syst_framework, syst_mol, gridstep)

    atomdict = IdDict{Symbol,Int}()
    atoms = syst_mol[:,:atomic_symbol]
    for atom in atoms
        get!(Returns(length(atomdict)+1), atomdict, atom)
    end
    atomsidx = [atomdict[atom] for atom in atoms]
    rev_atomdict = Vector{Symbol}(undef, length(atomdict))
    for (at, i) in atomdict
        rev_atomdict[i] = at
    end

    coulomb_grid_path, vdws = grid_locations(framework, forcefield_framework, rev_atomdict, gridstep, supercell)

    needcoulomb = any(!iszero(syst_mol[i,:atomic_charge])::Bool for i in 1:length(syst_mol))
    coulomb, ewald = if needcoulomb
        _eframework = initialize_ewald(syst_framework, supercell)
        if framework isa AbstractMatrix
            EnergyGrid(true)
        else
            retrieve_or_create_grid(coulomb_grid_path, syst_framework, forcefield, gridstep, _eframework, mat, new, cutoff)
        end, _eframework
    else
        EnergyGrid(true), EwaldFramework(mat)
    end

    grids = Vector{EnergyGrid}(undef, length(atomdict))
    for (atom, i) in atomdict
        grids[i] = if framework isa AbstractMatrix
            EnergyGrid(true)
        else
            retrieve_or_create_grid(vdws[i], syst_framework, forcefield, gridstep, atom, mat, new, cutoff)
        end
    end

    charges = [ustrip(u"e_au", syst_mol[i,:atomic_charge])::Float64 for i in 1:length(syst_mol)]

    CrystalEnergySetup(syst_framework, syst_mol, coulomb, charges, grids, atomsidx, ewald, forcefield, block)
end
function setup_RASPA(framework, forcefield_framework, molecule, forcefield_molecule; gridstep=0.15u"Å", supercell=nothing, blockfile=nothing, new=false, cutoff=12.0u"Å")
    syst_framework = load_framework_RASPA(framework, forcefield_framework)
    syst_mol = load_molecule_RASPA(molecule, forcefield_molecule, forcefield_framework, syst_framework)
    setup_RASPA(framework, forcefield_framework, syst_mol; gridstep, supercell, blockfile, new, cutoff)
end
function setup_RASPA(framework, forcefield, atom::Union{AbstractString,Symbol}; gridstep=0.15u"Å", supercell=nothing, blockfile=nothing, new=false, cutoff=12.0u"Å")
    syst_mol = default_system(atom, forcefield)
    setup_RASPA(framework, forcefield, syst_mol; gridstep, supercell, blockfile, new, cutoff)
end


function parse_interaction_RASPA(l, mixingrules::Bool, shift::Bool, cutoff, tailcorrection::Bool)
    splits = split(l)
    for (i, x) in enumerate(splits)
        if x[1] == '#' || (length(x)>1 && x[1] == '/' && x[2] == '/')
            resize!(splits, i-1)
            break
        end
    end
    kind = lowercase(splits[3-mixingrules])
    (Symbol(splits[1]), Symbol(splits[2-mixingrules])) => if kind == "lennard-jones"
        InteractionRule(FF.LennardJones, parse.(Float64, splits[(4-mixingrules):end]), shift, cutoff, tailcorrection)
    elseif kind == "none"
        FF.NoInteraction()
    elseif kind == "buckingham"
        InteractionRule(FF.Buckingham, parse.(Float64, splits[(4-mixingrules):end]), shift, cutoff, tailcorrection)
    elseif kind == "buckingham2"
        InteractionRuleSum([InteractionRule(FF.Buckingham, parse.(Float64, splits[(4-mixingrules):(6-mixingrules)]), shift, cutoff, tailcorrection),
                            InteractionRule(FF.HardSphere, [parse(Float64, splits[(7-mixingrules)]), 0.0], shift, cutoff, tailcorrection)])
    else
        error(lazy"$kind interaction potential not implemented")
    end
end
function parse_mixingrule_RASPA(l)
    x = lowercase(l)
    if x == "lorentz-berthelot"
        FF.LorentzBerthelot
    elseif x == "jorgensen" || x == "good-hope" || x == "geometric"
        FF.Geometric
    elseif x == "WaldmanHagler"
        FF.WaldmanHagler
    else
        error(lazy"Unknown mixing rule $l")
    end
end
function parse_shift(s)
    s2 = lowercase(s)
    s2 == "shifted" && return true
    @assert s2 == "truncated"
    return false
end
function parse_yesno(s)
    s2 = lowercase(s)
    s2 == "yes" && return true
    @assert s2 == "no"
    return false
end
function forcefield_indices(x, y, sdict)
    haskey(sdict, x) || error(lazy"Atom $x absent from force_field_mixing_rules.def")
    haskey(sdict, y) || error(lazy"Atom $y absent from force_field_mixing_rules.def")
    sdict[x], sdict[y]
end
function next_noncomment_line(lines, i)
    j = i+1
    while lines[j][1] == '#'
        j += 1
    end
    j
end

"""
    parse_forcefield_RASPA(name, [pseudoatoms::PseudoAtomListing]; cutoff=12.0u"Å", ewald_precision=!isinf(cutoff)*1e-6u"Å^-1")

Return a [`ForceField`](@ref) parsed from that called `name` in the RASPA directory. If
separately precomputed, the [`PseudoAtomListing`] can be provided.

`cutoff` specifies the interaction cutoff.

`ewald_precision` is the precision that should be used for the Ewald summation technique.
By default, it is set to 0.0 (which means not to use Ewald summation) if there is no
cutoff; when set to 0.0, the returned force field will include direct [`Coulomb`](@ref)
pair interactions between all charged species.
Otherwise, the default is to set `ewald_precision` to 1e-6 Å⁻¹. When set to a non-zero
value, the returned force field will include [`CoulombEwaldDirect`](@ref) pair interactions
between all charged species, which corresponds to the direct part of the Ewald summation.
"""
function parse_forcefield_RASPA(name, pseudoatoms_input::Union{Nothing,PseudoAtomListing}=nothing; cutoff=12.0u"Å", ewald_precision=!isinf(cutoff)*1e-6u"Å^-1")
    if !(name isa AbstractString)
        _ff = name isa ForceField ? name : name[2]
        if _ff isa ForceField
            _ff.cutoff == cutoff || error("Inconsistent cutoffs: choose between $(_ff.cutoff) from force field or given $cutoff")
            return _ff
        end
    end
    pseudoatoms = if pseudoatoms_input isa Nothing
        parse_pseudoatoms_RASPA(name, nothing)
    else
        if !(name isa AbstractString) && name[2] isa PseudoAtomListing && hash_string(name[2]) != hash_string(pseudoatoms_input)
            error("Inconsistent PseudoAtomListing given.")
        end
        pseudoatoms_input
    end::PseudoAtomListing
    rawname = name isa AbstractString ? name : name[1]
    input = Pair{Tuple{Symbol,Symbol},Union{InteractionRule,InteractionRuleSum}}[]
    general_mixingrule = FF.ErrorOnMix
    shift = true
    tailcorrection = false
    sdict = IdDict{Symbol,Int}()
    ff_mixing_rules_def = joinpath(RASPADIR[], "forcefield", rawname, "force_field_mixing_rules.def")
    if isfile(ff_mixing_rules_def)
        lines_ffmr = readlines(ff_mixing_rules_def)
        idx = next_noncomment_line(lines_ffmr, 1)
        shift = parse_shift(lines_ffmr[idx])
        idx = next_noncomment_line(lines_ffmr, idx)
        tailcorrection = parse_yesno(lines_ffmr[idx])
        idx = next_noncomment_line(lines_ffmr, idx)
        num_interactions_ffmr = parse(Int, lines_ffmr[idx])
        for i in 1:num_interactions_ffmr
            idx = next_noncomment_line(lines_ffmr, idx)
            interaction = parse_interaction_RASPA(lines_ffmr[idx], true, shift, cutoff, tailcorrection)
            sdict[interaction[1][1]] = i
            push!(input, interaction)
        end
        idx = next_noncomment_line(lines_ffmr, idx)
        general_mixingrule = parse_mixingrule_RASPA(lines_ffmr[idx])
    end
    ff = ForceField(input, general_mixingrule, cutoff, shift, tailcorrection; sdict, name=rawname)
    ff_def = joinpath(RASPADIR[], "forcefield", rawname, "force_field.def")
    if isfile(ff_def)
        lines_ff = readlines(ff_def)
        jdx = next_noncomment_line(lines_ff, 1)
        numnewrules = parse(Int, lines_ff[jdx])
        startnewrules = jdx
        jdx = next_noncomment_line(lines_ff, jdx)
        for _ in 1:numnewrules
            jdx = next_noncomment_line(lines_ff, jdx)
        end
        numnewinteractions = parse(Int, lines_ff[jdx])
        for _ in 1:numnewinteractions
            jdx = next_noncomment_line(lines_ff, jdx)
            (x1, y1), newinteraction = parse_interaction_RASPA(lines_ff[jdx], false, shift, cutoff, tailcorrection)
            i1, j1 = forcefield_indices(x1, y1, sdict)
            ff.interactions[i1,j1] = ff.interactions[j1,i1] = newinteraction
        end
        jdx = next_noncomment_line(lines_ff, jdx)
        numnewmixing = parse(Int, lines_ff[jdx])
        for _ in 1:numnewmixing
            jdx = next_noncomment_line(lines_ff, jdx)
            _x2, _y2, _newmix = split(lines_ff[jdx])
            x2 = Symbol(_x2); y2 = Symbol(_y2)
            i2, j2 = forcefield_indices(x2, y2, sdict)
            newmix = parse_mixingrule_RASPA(_newmix)
            ff.interactions[i2,j2] = ff.interactions[j2,i2] = mix_rules(ff.interactions[i2,i2], ff.interactions[j2,j2], newmix)
        end
        jdx = startnewrules
        for _ in 1:numnewrules
            jdx = next_noncomment_line(lines_ff, jdx)
            _x3, _y3, _newshift, _newtailcorrection = split(lines_ff[jdx])
            x3 = Symbol(_x3); y3 = Symbol(_y3)
            newshift = parse_shift(_newshift)
            newtailcorrection = parse_yesno(_newtailcorrection)
            i3, j3 = forcefield_indices(x3, y3, sdict)
            ff.interactions[i3,j3] = ff.interactions[j3,i3] = map_rule(ff.interactions[i3,j3]) do r
                InteractionRule(r.kind, r.params, newshift, cutoff, newtailcorrection)
            end
        end
    end
    α = if iszero(ewald_precision)
        0.0
    else
        ε = log(Float64(ewald_precision*cutoff))
        sqrt(abs(ε + log(sqrt(abs(ε)))))/Float64(cutoff/u"Å")
    end
    for (ati4, i4) in sdict
        chargei = pseudoatoms[ati4].charge
        chargei == 0.0 && continue
        for (atj4, j4) in sdict
            chargej = pseudoatoms[atj4].charge
            chargej == 0.0 && continue
            ff.interactions[i4,j4] = sum_rules(ff.interactions[i4,j4], if iszero(ewald_precision)
                InteractionRule(FF.Coulomb, [chargei, chargej], 0.0, false)
            else
                InteractionRule(FF.CoulombEwaldDirect, [α, chargei, chargej], 0.0, false)
            end)
        end
    end
    return ForceField(ff.interactions, ff.sdict, ff.symbols, ff.cutoff, rawname)
end

function setup_probe_RASPA(framework, forcefield_framework, atom)
    system = load_framework_RASPA(framework, forcefield_framework)
    forcefield = parse_forcefield_RASPA(forcefield_framework)
    ProbeSystem(system::AbstractSystem{3}, forcefield::ForceField, atom::Symbol)
end

function read_restart_RASPA(file)
    open(file) do io
        for _ in 1:29; readline(io); end
        numcomponents = parse(Int, split(readline(io))[2])
        for _ in 1:(4+11*numcomponents); readline(io); end
        positions = Vector{Vector{SVector{3,TÅ}}}[]
        l = readline(io)
        while !isempty(l)
            num_mol = parse(Int, split(l)[4])
            readline(io)
            pos = position(io)
            readline(io)
            num_atoms = 1
            while true
                parse(Int, split(readline(io))[3]) == 0 && break
                num_atoms += 1
            end
            newpos = [Vector{SVector{3,TÅ}}(undef, num_atoms) for _ in 1:num_mol]
            push!(positions, newpos)
            seek(io, pos)
            for j in 1:num_mol, k in 1:num_atoms
                splits = split(readline(io))
                @assert j == parse(Int, splits[2]) + 1
                @assert k == parse(Int, splits[3]) + 1
                newpos[j][k] = SVector{3,Float64}(parse.(Float64, @view splits[4:6]))*u"Å"
            end
            for _ in 1:(5*num_mol*num_atoms+1); readline(io); end
            l = readline(io)
        end
        positions
    end
end
