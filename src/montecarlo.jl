export MonteCarloSetup, setup_montecarlo, baseline_energy, movement_energy, run_montecarlo!

struct MonteCarloSetup{N,T,Trng}
    step::SimulationStep{N,T}
    # step contains all the information that is not related to the framework nor to Ewald.
    # It contains all the information necessary to compute the species-species VdW energy.
    tailcorrection::Base.RefValue{TK}
    coulomb::EnergyGrid
    grids::Vector{EnergyGrid} # grids[ix] is the VdW grid for atom ix in ff.
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


struct IdSystem # pseudo-AbstractSystem with only atomic symbols and charges
    atomic_symbol::Vector{Symbol}
    atomic_charge::Vector{Te_au}
end
Base.length(s::IdSystem) = length(s.atomic_charge)

function setup_montecarlo(cell::CellMatrix, csetup::GridCoordinatesSetup, systems,
                          num_framework_atoms::Vector{Int};
                          parallel::Bool=true, rng=default_rng(), mcmoves::Vector)
    if any(≤(24.0u"Å"), perpendicular_lengths(cell.mat))
        error("The current cell has at least one perpendicular length lower than 24.0Å: please use a larger supercell")
    end
    @assert systems == length(mcmoves)

    coulomb = EnergyGrid()
    grids = EnergyGrid[]

    U = Vector{SVector{3,TÅ}} # positions of the atoms of a system
    poss = Vector{U}[U[[rand(SVector{3,TÅ})]] for _ in 1:systems]
    indices = Tuple{Int,Int}[(1, 1)]
    speciesblocks = BlockFile[BlockFile(csetup)]
    newmcmoves = MCMoves[MCMoves(true)]

    ffidx = [[1] for _ in 1:systems]
    charges = [NaN*u"e_au"]

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
            tcorrection += ni*nj*3.1u"K"
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

    beads = fill(1, n)

    atomblocks = BlockFile[]


    MonteCarloSetup(SimulationStep(ForceField(), charges, poss, trues(length(ffidx)), ffidx, cell; parallel),
                    Ref(tcorrection), coulomb, grids, offsets, Set(indices_list),
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
                          parallel=false,
                          mcmoves=fill(nothing, systems), rng=default_rng())
    syst_framework = load_framework_RASPA(framework, forcefield_framework)
    supercell = find_supercell(syst_framework, 12.0u"Å")
    cell = CellMatrix(SMatrix{3,3,TÅ,9}(stack(bounding_box(syst_framework).*supercell)))

    csetup = GridCoordinatesSetup(syst_framework, 0.15u"Å")

    Π = prod(supercell)
    num_framework_atoms = [Π*length(syst_framework)]

    setup_montecarlo(cell, csetup, systems, num_framework_atoms; parallel, rng, mcmoves)
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
                    Ref(mc.tailcorrection[]), deepcopy(mc.coulomb), deepcopy(mc.grids), copy(mc.offsets),
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


function randomize_position!(positions, rng, indices, bead, d)
    pos = @view positions[indices]
    posr = random_rotation(rng, random_rotation(rng, random_rotation(rng, pos, 90u"°", bead, 1), 90u"°", bead, 2), 90u"°", bead, 3)
    post = random_translation(rng, posr, d)
    positions[indices] .= post
    post
end
