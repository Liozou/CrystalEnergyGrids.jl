export MonteCarloSetup, setup_montecarlo, baseline_energy, movement_energy, run_montecarlo!

struct MonteCarloSetup{N,T,Trng}
    step::SimulationStep{N,T}
    # step contains all the information that is not related to the framework nor to Ewald.
    # It contains all the information necessary to compute the species-species VdW energy.
    tailcorrection::Base.RefValue{TK}
    offsets::Vector{Int}
    # offsets[i] is the number of molecules belonging to a kind strictly lower than i.
    indices::Set{Tuple{Int,Int}}
    # indices is the set of all (i,j) with 1 ≤ j ≤ number of species of kind i.
    speciesblocks::Vector{BlockFile}
    # speciesblocks[i][pos] is set if pos is blocked for any atom of a species of kind i.
    atomblocks::Vector{BlockFile}
    # atomblock[ix][pos] is set if pos is blocked for all atoms of index ix in ff.
    bead::Vector{Int} # k = bead[i] is the number of the reference bead of kind i.
    mcmoves::Vector{MCMoves} # mcmoves[i] is the set of MC moves for kind i
    rng::Trng
end

function Base.show(io::IO, mc::MonteCarloSetup)
    n = length(mc.indices)
    m = length(mc.offsets)
    print(io, "Monte-Carlo setup with ", n , " atoms in ", m, " molecule kind")
    m > 1 && print(io, 's')
end


struct EwaldSystem # pseudo-AbstractSystem with only positions and charges
    position::Vector{SVector{3,TÅ}}
    atomic_charge::Vector{Te_au}
end
AtomsBase.position(x::EwaldSystem) = x.position
AtomsBase.position(x::EwaldSystem, i::Int) = x.position[i]
Base.getindex(x::EwaldSystem, ::Colon, s::Symbol) = getproperty(x, s)
Base.getindex(x::EwaldSystem, i::Int, s::Symbol) = getproperty(x, s)[i]
Base.length(x::EwaldSystem) = length(x.position)

struct IdSystem # pseudo-AbstractSystem with only atomic symbols and charges
    atomic_symbol::Vector{Symbol}
    atomic_charge::Vector{Te_au}
end
IdSystem(s::AbstractSystem{3}) = IdSystem(atomic_symbol(s), s[:,:atomic_charge])
Base.length(s::IdSystem) = length(s.atomic_charge)

function setup_montecarlo(cell::CellMatrix, csetup::GridCoordinatesSetup,
                          systems::Vector, ff::ForceField,
                          coulomb::EnergyGrid, grids::Vector{EnergyGrid},
                          blocksetup::Vector{String},
                          num_framework_atoms::Vector{Int},
                          restartpositions::Union{Nothing,Vector{Vector{Vector{SVector{3,TÅ}}}}};
                          parallel::Bool=true, rng=default_rng(), mcmoves::Vector)
    if any(≤(24.0u"Å"), perpendicular_lengths(cell.mat))
        error("The current cell has at least one perpendicular length lower than 24.0Å: please use a larger supercell")
    end
    @assert length(systems) == length(blocksetup) == length(mcmoves)

    kindsdict = Dict{Tuple{Vector{Symbol},String,MCMoves},Int}()
    systemkinds = IdSystem[]
    U = Vector{SVector{3,TÅ}} # positions of the atoms of a system
    poss = restartpositions isa Nothing ? Vector{U}[] : restartpositions
    indices = Tuple{Int,Int}[]
    rev_indices = Vector{Int}[]
    speciesblocks = BlockFile[]
    newmcmoves = MCMoves[]
    for (i, (s, block)) in enumerate(zip(systems, blocksetup))
        system, n = s isa Tuple ? s : (s, 1)
        m = length(kindsdict)+1
        mcmove = isnothing(mcmoves[i]) ? MCMoves(length(systems) == 1) : mcmoves[i]
        kind = get!(kindsdict, (atomic_symbol(system)::Vector{Symbol}, block, mcmove), m)
        if kind === m
            push!(systemkinds, IdSystem(system)::IdSystem)
            restartpositions isa Nothing && push!(poss, U[])
            push!(rev_indices, Int[])
            if isempty(block)
                push!(speciesblocks, BlockFile(csetup))
            else
                push!(speciesblocks, parse_blockfile(joinpath(getdir_RASPA(), "structures", "block", block)*".block", csetup))
            end
            push!(newmcmoves, mcmove)
        end
        push!(indices, (kind, length(poss[kind])+1))
        append!(rev_indices[kind], i for _ in 1:n)
        restartpositions isa Nothing && append!(poss[kind], copy(position(system)::Vector{SVector{3,TÅ}}) for _ in 1:n)
    end

    ffidx = [[ff.sdict[s.atomic_symbol[k]] for k in 1:length(s)] for s in systemkinds]
    charges = fill(NaN*u"e_au", length(ff.sdict))
    for (i, ids) in enumerate(ffidx)
        kindi = systemkinds[i]
        for (k, ix) in enumerate(ids)
            charges[ix] = kindi.atomic_charge[k]
        end
    end

    num_atoms = copy(num_framework_atoms)
    for (i, indices) in enumerate(ffidx)
        mult = length(poss[i])
        for id in indices
            num_atoms[id] += mult
        end
    end
    tcorrection = 0.0u"K"
    for (i, ni) in enumerate(num_atoms)
        ni == 0 && continue
        for (j, nj) in enumerate(num_atoms)
            nj == 0 && continue
            tcorrection += ni*nj*tailcorrection(ff[i,j], ff.cutoff)
        end
    end
    tcorrection *= 2π/det(ustrip.(u"Å", cell.mat))

    n = length(poss)
    offsets = Vector{Int}(undef, n)
    offsets[1] = 0
    indices_list = Tuple{Int,Int}[]
    for i in 1:(n-1)
        numi = length(poss[i])
        offsets[i+1] = offsets[i] + numi
        append!(indices_list, (i,j) for j in 1:length(poss[i]))
    end
    append!(indices_list, (n,j) for j in 1:length(poss[n]))

    buffer = MVector{3,TÅ}(undef)
    buffer2 = MVector{3,Float64}(undef)
    beads = Vector{Int}(undef, n)
    for i in 1:n
        if isempty(poss[i])
            beads[i] = 0
            continue
        end
        possi = poss[i][1]
        refpos = possi[1]
        meanpos = refpos + mean(pos - refpos for pos in possi)
        mind2 = Inf*u"Å^2"
        thisbead = 0
        for (k, pos) in enumerate(possi)
            buffer .= pos .- meanpos
            d2 = unsafe_periodic_distance2!(buffer, buffer2, cell)
            if d2 < mind2
                mind2 = d2
                thisbead = k
            end
        end
        beads[i] = thisbead
    end

    atomblocks = Vector{BlockFile}(undef, length(grids))
    @inbounds for ix in 1:length(grids)
        if isassigned(grids, ix)
            atomblocks[ix] = BlockFile(grids[ix])
        end
    end

    ewaldsystems = EwaldSystem[]
    for (i, positioni) in enumerate(poss)
        ffidxi = ffidx[i]
        charge = [charges[k] for k in ffidxi]
        rev_idx = rev_indices[i]
        bead = beads[i]
        block = speciesblocks[i]
        d = norm(sum(cell.mat; dims=2))/4
        for j in 1:length(positioni)
            s = systems[rev_idx[j]]
            n = restartpositions isa Nothing && s isa Tuple ? s[2]::Int : 1
            if n > 1
                randomize_position!(positioni[j], rng, 1:length(ffidxi), bead, block, ffidxi, atomblocks, d)
            end
            push!(ewaldsystems, EwaldSystem(positioni[j], charge))
        end
    end

    MonteCarloSetup(SimulationStep(ff, charges, poss, trues(length(ffidx)), ffidx, cell; parallel),
                    Ref(tcorrection), offsets, Set(indices_list),
                    speciesblocks, atomblocks, beads, newmcmoves, rng), indices
end

"""
    setup_montecarlo(framework, forcefield_framework::String, systems;
                     blockfiles=fill(nothing, length(systems)), gridstep=0.15u"Å",
                     supercell=nothing, new=false, restart=nothing, parallel=true,
                     mcmoves=fill(nothing, length(systems)), rng=default_rng())

Prepare a Monte Carlo simulation of the input list of `systems` in a `framework` with the
given force field.

A system can be either an `AbstractSystem`, or a pair `(s, n)` where `s` is an
`AbstractSystem` and `n` is an integer. In that case, `n` systems identical to `s` will be
added and their position and orientations randomized.

`blockfiles[i]` can be set to `false` to allow the molecule `systems[i]` to go everywhere in
the framework.
Or it can be set to the radical (excluding the ".block" extension) of the block file in the
raspa directory to include it and prevent the molecule to go in the blocked spheres.
Setting it to `true` is equivalent to using `blockfiles[i]=first(split(framework, '_'))`.
If not provided, the default is to set it to `false` if the molecule is positively charged
and monoatomic (considered to be a small cation), or to `true` otherwise.

`gridstep` is the approximate energy grid step size used.

`supercell` is a triplet of integer specifying the number of repetitions of unit cell to
use. If unspecified, it is the smallest supercell such that the perpendicular lengths are
all above 24.0 Å (i.e. twice the 12.0 Å cutoff).

If `new` is set, force recomputing all the grids. Otherwise, existing grids will be used
when available.

If `restart` is the path of a raspa restart file, the positions will be taken from the file.

`parallel` mirrors the corresponding argument to [`SimulationStep`](@ref).

`mcmoves` is the list of [`MCMoves`](@ref) for each system. The default value of `nothing`
leads to the default move probabilities: [98% translation, 2% random translation] for a
monoatomic species, or [49% translation, 49% rotation, 2% random reinsertion] else.

`rng` is the random number generator, defaulting to `Random.default_rng()`.
"""
function setup_montecarlo(framework, forcefield_framework::String, systems;
                          blockfiles=fill(nothing, length(systems)), gridstep=0.15u"Å",
                          supercell=nothing, new=false, restart=nothing, parallel=true,
                          mcmoves=fill(nothing, length(systems)), rng=default_rng())
    syst_framework = load_framework_RASPA(framework, forcefield_framework)
    ff = parse_forcefield_RASPA(forcefield_framework)
    mat = stack3(bounding_box(syst_framework))
    if supercell isa Nothing
        supercell = find_supercell(syst_framework, 12.0u"Å")
    end
    supercell::NTuple{3,Int}
    cell = CellMatrix(SMatrix{3,3,TÅ,9}(stack(bounding_box(syst_framework).*supercell)))

    csetup = GridCoordinatesSetup(syst_framework, gridstep)
    length(blockfiles) == length(systems) || error("Please provide one blockfiles element per system")
    blocksetup = [decide_parse_block(file, s isa Tuple ? s[1] : s, framework) for (file, s) in zip(blockfiles, systems)]

    needcoulomb = false
    encountered_atoms = Set{Symbol}()
    for s in systems
        system = s isa Tuple ? s[1] : s
        if !needcoulomb
            needcoulomb = any(!iszero(system[i,:atomic_charge])::Bool for i in 1:length(system))
        end
        for k in 1:length(system)
            push!(encountered_atoms, atomic_symbol(system, k))
        end
    end
    atoms = [(atom, ff.sdict[atom]) for atom in encountered_atoms]
    sort!(atoms; by=last)
    # atoms is the list of unique atom types encountered in the molecule species

    coulomb_grid_path, vdw_grid_paths = grid_locations(framework, forcefield_framework, first.(atoms), gridstep, supercell)

    coulomb = if needcoulomb
        retrieve_or_create_grid(coulomb_grid_path, syst_framework, ff, gridstep, nothing, mat, new)
    else
        EnergyGrid()
    end

    grids = if isempty(vdw_grid_paths) || isempty(framework)
        # if framework is empty, signal it by having an empty list of grids
        EnergyGrid[]
    else
        _grids = Vector{EnergyGrid}(undef, length(ff.sdict))
        for (j, (atom, i)) in enumerate(atoms)
            _grids[i] = retrieve_or_create_grid(vdw_grid_paths[j], syst_framework, ff, gridstep, atom, mat, new)
        end
        _grids
    end

    num_framework_atoms = zeros(Int, length(ff.sdict))
    Π = prod(supercell)
    for at in syst_framework
        num_framework_atoms[ff.sdict[Symbol(get_atom_name(atomic_symbol(at)))]] += Π
    end

    restartpositions = restart isa Nothing ? nothing : read_restart_RASPA(restart)

    setup_montecarlo(cell, csetup, systems, ff, coulomb, grids, blocksetup, num_framework_atoms, restartpositions; parallel, rng, mcmoves)
end

"""
    MonteCarloSetup(mc::MonteCarloSetup; parallel::Bool=mc.step.parallel)

Create a copy of `mc` that does not share its modifiable internal states (the positions and
the Ewald state). For example, the copy and the original can be used to run Monte-Carlo
simulations in parallel, and the state of one will not influence that of the other.

`parallel` specifies whether the computations on the resulting `MonteCarloSetup` should be
parallelized or not.

!!! warn
    Internal states that are semantically immutable are shared, although some of them are
    technically mutable, like the value of the atom charges for instance. Modifying such a
    state on the original or the copy can thus also modify that of the other: use
    `deppcopy` or a custom copy implementation to circumvent this issue if you plan on
    modifying states outside of the API.
"""
function MonteCarloSetup(mc::MonteCarloSetup, o::SimulationStep=mc.step; parallel::Bool=mc.step.parallel)
    rng = deepcopy(mc.rng)
    MonteCarloSetup(SimulationStep(o, :all; parallel),
                    Ref(mc.tailcorrection[]), copy(mc.offsets),
                    copy(mc.indices), deepcopy(mc.speciesblocks), deepcopy(mc.atomblocks), copy(mc.bead),
                    copy(mc.mcmoves), rng)
end

function set_position!(mc::MonteCarloSetup, (i, j), newpositions, newEiks=nothing)
    molpos = mc.step.posidx[i][j]
    for (k, newpos) in enumerate(newpositions)
        mc.step.positions[molpos[k]] = if eltype(newpositions) <: AbstractVector{<:AbstractFloat}
            newpos
        else
            NoUnits(newpos/u"Å")
        end
    end
end


choose_random_species(mc::MonteCarloSetup) = rand(mc.rng, mc.indices)


function randomize_position!(positions, rng, indices, bead, block, ffidxi, atomblocks, d)
    pos = @view positions[indices]
    posr = random_rotation(rng, random_rotation(rng, random_rotation(rng, pos, 90u"°", bead, 1), 90u"°", bead, 2), 90u"°", bead, 3)
    post = random_translation(rng, posr, d)
    positions[indices] .= post
    post
end
