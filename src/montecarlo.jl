export MonteCarloSetup, setup_montecarlo, baseline_energy, movement_energy, run_montecarlo!

struct MonteCarloSetup{T,Trng}
    step::SimulationStep{T}
    # step contains all the information that is not related to the framework nor to Ewald.
    # It contains all the information necessary to compute the species-species VdW energy.
    ewald::IncrementalEwaldContext
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
                          systems::Vector, ff::ForceField, eframework::EwaldFramework,
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

    ewald = IncrementalEwaldContext(EwaldContext(eframework, ewaldsystems))

    MonteCarloSetup(SimulationStep(ff, charges, poss, trues(length(ffidx)), ffidx, cell; parallel),
                    ewald, Ref(tcorrection), coulomb, grids, offsets, Set(indices_list),
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

    coulomb, eframework = if needcoulomb
        _eframework = initialize_ewald(syst_framework)
        retrieve_or_create_grid(coulomb_grid_path, syst_framework, ff, gridstep, _eframework, mat, new), _eframework
    else
        EnergyGrid(), EwaldFramework(mat)
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

    setup_montecarlo(cell, csetup, systems, ff, eframework, coulomb, grids, blocksetup, num_framework_atoms, restartpositions; parallel, rng, mcmoves)
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
    ewaldsystems = EwaldSystem[]
    for (i, posidxi) in enumerate(o.posidx)
        ffidxi = o.ffidx[i]
        charge = [o.charges[k] for k in ffidxi]
        append!(ewaldsystems, EwaldSystem(o.positions[molpos], charge) for molpos in posidxi)
    end
    ewald = IncrementalEwaldContext(EwaldContext(mc.ewald.ctx.eframework, ewaldsystems))
    rng = mc.rng isa TaskLocalRNG ? mc.rng : copy(mc.rng)
    MonteCarloSetup(SimulationStep(o, :all; parallel),
                    ewald, Ref(mc.tailcorrection[]), mc.coulomb, mc.grids, mc.offsets,
                    copy(mc.indices), mc.speciesblocks, mc.atomblocks, mc.bead,
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

    oldEikx, oldEiky, oldEikz = mc.ewald.Eiks
    if newEiks isa Nothing
        move_one_system!(mc.ewald, mc.offsets[i] + j, newpositions)
    else
        newEikx, newEiky, newEikz = newEiks
        kx, ky, kz = mc.ewald.eframework.kspace.ks
        kxp = kx+1; tkyp = ky+ky+1; tkzp = kz+kz+1
        jofs = 1 + mc.ewald.offsets[mc.offsets[i] + j]
        copyto!(oldEikx, 1 + jofs*kxp, newEikx, 1, length(newEikx))
        copyto!(oldEiky, 1 + jofs*tkyp, newEiky, 1, length(newEikx))
        copyto!(oldEikz, 1 + jofs*tkzp, newEikz, 1, length(newEikx))
    end
    nothing
end


struct FrameworkEnergyReport
    vdw::TK
    direct::TK
end
FrameworkEnergyReport() = FrameworkEnergyReport(0.0u"K", 0.0u"K")
Base.Float64(f::FrameworkEnergyReport) = Float64((f.vdw + f.direct)/u"K")

struct MCEnergyReport
    framework::FrameworkEnergyReport
    inter::TK
    reciprocal::TK
end
MCEnergyReport(v, d, i, r) = MCEnergyReport(FrameworkEnergyReport(v, d), i, r)
Base.Float64(e::MCEnergyReport) = Float64(e.framework) + Float64((e.inter + e.reciprocal)/u"K")
Base.show(io::IO, e::MCEnergyReport) = show(io, Float64(e)*u"K")
function Base.show(io::IO, ::MIME"text/plain", e::MCEnergyReport)
    println(io, Float64(e), " = [", e.framework.vdw, " (f VdW) + ", e.framework.direct, " (f direct)] + [", e.inter, " (internal) + ", e.reciprocal, " (reciprocal)]")
end

struct BaselineEnergyReport
    er::MCEnergyReport
    tailcorrection::TK
end
BaselineEnergyReport(f, i, r, t) = BaselineEnergyReport(MCEnergyReport(f, i, r), t)
BaselineEnergyReport(v, d, i, r, t) = BaselineEnergyReport(MCEnergyReport(v, d, i, r), t)
Base.Float64(ber::BaselineEnergyReport) = Float64(ber.er) + Float64(ber.tailcorrection/u"K")
Base.show(io::IO, ber::BaselineEnergyReport) = show(io, Float64(ber)*u"K")
function Base.show(io::IO, ::MIME"text/plain", b::BaselineEnergyReport)
    println(io, Float64(b), " = [", b.er.framework.vdw, " (f VdW) + ", b.er.framework.direct, " (f direct)] + [", b.er.inter, " (internal) + ", b.er.reciprocal, " (reciprocal)] + ", b.tailcorrection, " (tailcorrection)")
end

for op in (:+, :-)
    @eval begin
        function Base.$(op)(f1::FrameworkEnergyReport, f2::FrameworkEnergyReport)
            FrameworkEnergyReport($op(f1.vdw, f2.vdw), $op(f1.direct, f2.direct))
        end
        function Base.$(op)(e1::MCEnergyReport, e2::MCEnergyReport)
            MCEnergyReport($op(e1.framework, e2.framework), $op(e1.inter, e2.inter), $op(e1.reciprocal, e2.reciprocal))
        end
        function Base.$(op)(b1::BaselineEnergyReport, e2::MCEnergyReport)
            BaselineEnergyReport($op(b1.er, e2), b1.tailcorrection)
        end
    end
end


function framework_interactions(grids::Vector{EnergyGrid}, coulombgrid::EnergyGrid, charges::Vector{Te_au}, indices::Vector{Int}, positions)
    isempty(grids) && return FrameworkEnergyReport()
    vdw = 0.0u"K"
    direct = 0.0u"K"
    hascoulomb = coulombgrid.ewald_precision != -Inf
    for (k, pos) in enumerate(positions)
        ix = indices[k]
        vdw += interpolate_grid(grids[ix], pos)
        if hascoulomb
            coulomb = interpolate_grid(coulombgrid, pos)
            direct += ifelse(coulomb == 1e100u"K", coulomb, Float64(charges[ix]/u"e_au")*coulomb)
        end
    end
    return FrameworkEnergyReport(vdw, direct)
end
function framework_interactions(mc::MonteCarloSetup, indices::Vector{Int}, positions)
    framework_interactions(mc.grids, mc.coulomb, mc.step.charges, indices, positions)
end

"""
    framework_interactions(mc::MonteCarloSetup, i, positions)

Energy contribution of the interaction between the framework and a molecule of system kind
`i` at the given `positions`.

Only return the Van der Waals and the direct part of the Ewald summation. The reciprocal
part can be obtained with [`compute_ewald`](@ref).
"""
function framework_interactions(mc::MonteCarloSetup, i::Int, positions)
    framework_interactions(mc, mc.step.ffidx[i], positions)
end

"""
    baseline_energy(mc::MonteCarloSetup)

Compute the energy of the current configuration.
"""
function baseline_energy(mc::MonteCarloSetup)
    reciprocal = compute_ewald(mc.ewald)
    vdw = compute_vdw(mc.step)
    fer = FrameworkEnergyReport()
    isempty(mc.grids) && return BaselineEnergyReport(fer, vdw, reciprocal, mc.tailcorrection[])
    for (i, indices) in enumerate(mc.step.ffidx)
        molposi = mc.step.posidx[i]
        for molpos in molposi
            fer += framework_interactions(mc, indices, @view mc.step.positions[molpos])
        end
    end
    return BaselineEnergyReport(fer, vdw, reciprocal, mc.tailcorrection[])
end

"""
    movement_energy(mc::MonteCarloSetup, (i, j), positions=nothing)

Compute the energy contribution of the `j`-th molecule of kind `i` when placed at
`positions`. If not provided, `positions` is the current position for that molecule.

The energy difference between the new position for the molecule and the current one is
`movement_energy(mc, (i,j), positions) - movement_energy(mc, (i,j))`.

!!! warning
    `baseline_energy(mc)` must have been called at least once before, otherwise the
    computation of the Ewald part will error.
    See also [`single_contribution_ewald`](@ref).
"""
function movement_energy(mc::MonteCarloSetup, idx, positions=nothing)
    i, j = idx
    ij = mc.offsets[i]+j
    poss = if positions isa Nothing
        molpos = mc.step.posidx[i][j]
        @view mc.step.positions[molpos]
    else
        positions
    end
    singlereciprocal = @spawn single_contribution_ewald($mc.ewald, $ij, $positions)
    fer = @spawn framework_interactions($mc, $i, $poss)
    singlevdw = single_contribution_vdw(mc.step, (i,j), poss)
    MCEnergyReport(fetch(fer), singlevdw, fetch(singlereciprocal))
end

"""
    update_mc!(mc::MonteCarloSetup, idx::Tuple{Int,Int}, positions::Vector{SVector{3,TÅ}})

Following a call to [`movement_energy(mc, idx, positions)`](@ref), update the internal
state of `mc` so that the species of index `idx` is now at `positions`.
"""
function update_mc!(mc::MonteCarloSetup, (i,j)::Tuple{Int,Int}, positions::Vector{SVector{3,TÅ}})
    update_ewald_context!(mc.ewald)
    L = mc.step.posidx[i][j]
    mc.step.positions[L] .= positions
    nothing
end


function inblockpocket(block::BlockFile, atomblocks::Vector{BlockFile}, ffidxi::Vector{Int}, newpos::Vector{SVector{3,TÅ}})
    for (j, pos) in enumerate(newpos)
        block[pos] && return true
        ablock = atomblocks[ffidxi[j]]
        ablock[pos + ablock.csetup.Δ./2] && return true
    end
    return false
end

choose_random_species(mc::MonteCarloSetup) = rand(mc.rng, mc.indices)
function compute_accept_move(before::MCEnergyReport, after::MCEnergyReport, T, rng)
    b = Float64(before)
    a = Float64(after)
    a < b && return true
    e = exp((b-a)*u"K"/T)
    return rand(rng) < e
end


function randomize_position!(positions, rng, indices, bead, block, ffidxi, atomblocks, d)
    pos = @view positions[indices]
    for _ in 1:30
        posr = random_rotation(rng, random_rotation(rng, random_rotation(rng, pos, 90u"°", bead, 1), 90u"°", bead, 2), 90u"°", bead, 3)
        for _ in 1:1000
            post = random_translation(rng, posr, d)
            if !inblockpocket(block, atomblocks, ffidxi, post)
                positions[indices] .= post
                return post
            end
        end
    end
    error(lazy"Could not find a suitable position for molecule $((i,j))!")
end

#=
"""
    randomize_position!(mc::MonteCarloSetup, idx, update_ewald=true)

Put the species at the given index to a random position and orientation.

!!! warn
    Using `update_ewald=false` will result in an inconsistent state for `mc` which will
    error on [`run_montecarlo`](@ref) unless [`baseline_energy(mc)`](@ref) is called.

!!! warn
    `update_ewald=true` requires a prior call to [`baseline_energy(mc)`](@ref).
"""
function randomize_position!(mc::MonteCarloSetup, (i,j), update_ewald=true)
    d = norm(sum(mc.step.mat; dims=2))/4
    newpos = randomize_position!(mc.step.positions, mc.step.posidx[i][j], mc.bead[i], mc.speciesblocks[i], mc.step.ffidx[i], mc.atomblocks, d)
    if update_ewald
        single_contribution_ewald(mc.ewald, mc.offsets[i] + j, newpos)
        update_mc!(mc, (i,j), newpos)
    else
        mc.ewald.last[] = -1 # signal the inconsistent state
    end
    newpos
end
=#
