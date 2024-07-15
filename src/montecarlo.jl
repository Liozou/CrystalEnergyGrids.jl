export MonteCarloSetup, setup_montecarlo, baseline_energy, movement_energy

struct MonteCarloSetup{T,Trng}
    step::SimulationStep{T}
    # step contains all the information that is not related to the framework nor to Ewald.
    # It contains all the information necessary to compute the species-species VdW energy.
    ewald::IncrementalEwaldContext
    flatidx::Vector{Vector{Int}}
    # flatidx[i][j] is the index ij used to refer to the j-th species of kind i in ewald.
    revflatidx::Vector{Tuple{Int,Int}} # revflatidx[ij] == (i,j) iff flatidx[i][j] == ij
    tailcorrection::TailCorrection
    coulomb::EnergyGrid
    grids::Vector{EnergyGrid} # grids[ix] is the VdW grid for atom ix in ff.
    speciesblocks::Vector{BlockFile}
    # speciesblocks[i][pos] is set if pos is blocked for any atom of a species of kind i.
    atomblocks::Vector{BlockFile}
    # atomblock[ix][pos] is set if pos is blocked for all atoms of index ix in ff.
    bead::Vector{Int} # k = bead[i] is the number of the reference bead of kind i.
    models::Vector{Vector{SVector{3,TÅ}}}
    # models[i] is a list of atom positions compatible with a species of kind i
    mcmoves::Vector{MCMoves} # mcmoves[i] is the set of MC moves for kind i
    gcmcdata::GCMCData
    rng::Trng
end

function Base.show(io::IO, mc::MonteCarloSetup)
    n = length(mc.revflatidx)
    m = length(mc.mcmoves)
    print(io, "Monte-Carlo setup with ", n , " molecules among ", m, " molecule kind")
    m > 1 && print(io, 's')
end

function Base.:(==)(mc1::MonteCarloSetup, mc2::MonteCarloSetup)
    mc1.step == mc2.step &&
    mc1.ewald == mc2.ewald &&
    mc1.flatidx == mc2.flatidx &&
    mc1.revflatidx == mc2.revflatidx &&
    mc1.tailcorrection == mc2.tailcorrection &&
    mc1.coulomb == mc2.coulomb &&
    length(mc1.grids) == length(mc2.grids) &&
    all(i -> begin a = isassigned(mc1.grids, i); a == isassigned(mc2.grids, i) && (!a || mc1.grids[i] == mc2.grids[i]) end, 1:length(mc1.grids)) &&
    mc1.speciesblocks == mc2.speciesblocks &&
    length(mc1.atomblocks) == length(mc2.atomblocks) &&
    all(i -> begin a = isassigned(mc1.atomblocks, i); a == isassigned(mc2.atomblocks, i) && (!a || mc1.atomblocks[i] == mc2.atomblocks[i]) end, 1:length(mc1.atomblocks)) &&
    mc1.bead == mc2.bead &&
    mc1.models == mc2.models &&
    mc1.mcmoves == mc2.mcmoves &&
    mc1.gcmcdata == mc2.gcmcdata
end

isdefined_ewald(mc::MonteCarloSetup) = isdefined_ewald(mc.ewald)

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
                          parallel::Bool=true, rng=default_rng(), mcmoves::Vector, pff=ff)
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
    models = Vector{SVector{3,TÅ}}[]
    idxsystem = 0
    # charge_flat is a signal that the system should be expanded based on charges
    charge_flag::Union{Nothing,NTuple{4,Int}} = nothing
    while idxsystem < length(systems)
        idxsystem += 1
        s = systems[idxsystem]
        block = blocksetup[idxsystem]
        _system, n = s isa Tuple ? s : (s, 1)
        system = default_system(_system, pff)
        m = length(kindsdict)+1
        mcmove = isnothing(mcmoves[idxsystem]) ? MCMoves(length(s) == 1) : mcmoves[idxsystem]
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
            push!(models, copy(position(system)::Vector{SVector{3,TÅ}}))
        end
        push!(indices, (kind, length(poss[kind])+1))
        append!(rev_indices[kind], idxsystem for _ in 1:n)
        restartpositions isa Nothing && append!(poss[kind], copy(position(system)::Vector{SVector{3,TÅ}}) for _ in 1:n)
        if n == -1
            charge_flag isa Nothing || error("At most one species can have their number specified as -1")
            charge_flag = (idxsystem, kind, length(rev_indices[kind]), length(poss[kind]))
        end
    end

    ffidx = [[ff.sdict[s.atomic_symbol[k]] for k in 1:length(s)] for s in systemkinds]
    charges = fill(NaN*u"e_au", length(ff.sdict))
    for (i, ids) in enumerate(ffidx)
        kindi = systemkinds[i]
        for (k, ix) in enumerate(ids)
            charges[ix] = kindi.atomic_charge[k]
        end
    end

    if charge_flag isa NTuple{4,Int}
        idxsystem′, kind′, i_revidxk′, i_possk′ = charge_flag
        _system′ = first(systems[idxsystem])
        system′ = default_system(_system′, pff)
        tot_charge = eframework.net_charges_framework + sum(sum(charges[ix] for ix in ffidx[i]; init=0.0u"e_au")*length(poss[i]) for i in 1:length(poss); init=0.0u"e_au")
        this_charge = sum(uconvert(u"e_au", system′[i,:atomic_charge])::Te_au for i in 1:length(system′); init=0.0u"e_au")
        iszero(this_charge) && error("Cannot use a number equal to -1 on a neutral species")
        n′ = round(Int, -tot_charge/this_charge)
        n′ ≥ 0 || error(lazy"Cannot compensate the total charge of $tot_charge with a species of the charge $this_charge")
        systems[idxsystem] = (_system′, n′)
        splice!(rev_indices[kind′], (i_revidxk′+1):i_revidxk′, idxsystem′ for _ in 1:n′)
        restartpositions isa Nothing && splice!(poss[kind′], (i_possk′+1):i_possk′, copy(position(system′)::Vector{SVector{3,TÅ}}) for _ in 1:n′)
    end

    cellvolume = det(cell.mat)
    λ = ustrip(u"Å^-3", 2π/cellvolume)
    tailcorrection = TailCorrection(ff, ffidx, num_framework_atoms, λ, length.(poss))

    n = length(poss)
    flatidx = Vector{Vector{Int}}(undef, length(poss))
    revflatidx = Tuple{Int,Int}[]
    for i in 1:n
        len = length(poss[i])
        flatidx[i] = collect((length(revflatidx)+1):(length(revflatidx)+len))
        append!(revflatidx, (i,j) for j in 1:len)
    end

    buffer = MVector{3,TÅ}(undef)
    buffer2 = MVector{3,Float64}(undef)
    beads = Vector{Int}(undef, n)
    for i in 1:n
        possi = models[i]
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
        models[i] .-= (models[i][thisbead],)
    end

    atomblocks = Vector{BlockFile}(undef, length(grids))
    @inbounds for ix in 1:length(grids)
        if isassigned(grids, ix)
            atomblocks[ix] = BlockFile(grids[ix])
        end
    end

    ewaldsystems = [EwaldSystem[] for _ in poss]
    for (i, positioni) in enumerate(poss)
        ffidxi = ffidx[i]
        charge = [charges[k] for k in ffidxi]
        rev_idx = rev_indices[i]
        bead = beads[i]
        block = speciesblocks[i]
        for j in 1:length(positioni)
            s = systems[rev_idx[j]]
            n = restartpositions isa Nothing && s isa Tuple ? s[2]::Int : 1
            if n > 1
                randomize_position!(positioni[j], rng, 1:length(ffidxi), bead, block, ffidxi, atomblocks, cell.mat)
            end
            push!(ewaldsystems[i], EwaldSystem(positioni[j], charge))
        end
    end

    # Collect empty entries separately to create a proper EwaldSystem
    emptysystems = findall(isempty, ewaldsystems)
    for i_empty in emptysystems
        push!(ewaldsystems[i_empty], EwaldSystem(models[i_empty], charges[ffidx[i_empty]]))
    end

    ewald = IncrementalEwaldContext(EwaldContext(eframework, ewaldsystems, emptysystems))

    gcmc = GCMCData(ff, ffidx, speciesblocks, eframework.mat)

    MonteCarloSetup(SimulationStep(ff, charges, poss, trues(length(ffidx)), ffidx, cell; parallel),
                    ewald, flatidx, revflatidx, tailcorrection, coulomb, grids,
                    speciesblocks, atomblocks, beads, models, newmcmoves, gcmc, rng), indices
end

"""
    setup_montecarlo(framework, pff, systems;
                     blockfiles=fill(nothing, length(systems)), gridstep=0.15u"Å",
                     supercell=nothing, new=false, restart=nothing, parallel=true,
                     mcmoves=fill(nothing, length(systems)), rng=default_rng(),
                     cutoff=12.0u"Å")

Prepare a Monte Carlo simulation of the input list of `systems` in a `framework` with the
given force field `pff`.

`pff` can be either a [`ForceField`](@ref), its name only, or a tuple with a name and
either a `ForceField` or a [`PseudoAtomListing`](@ref). If a `PseudoAtomListing` is
provided, it takes priority over the "pseudo_atoms.def" of the force field.

A system can be either an `AbstractSystem`, or a pair `(s, n)` where `s` is an
`AbstractSystem` and `n` is an integer. In that case, `n` systems identical to `s` will be
added and their position and orientations randomized.
Using `n = -1` is allowed for up to one species: in that case, the value of `n` will be
chosen so that the resulting system is electrically neutral. Using `n = -1` is thus
forbidden on neutral species.

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
all above twice the `cutoff`. The default cutoff being 12 Å, the default minimal
perpendicular length is 24 Å.

If `new` is set, force recomputing all the grids. Otherwise, existing grids will be used
when available.

If `restart` is the path of a raspa restart file, the positions will be taken from the file.

`parallel` mirrors the corresponding argument to [`SimulationStep`](@ref).

`mcmoves` is the list of [`MCMoves`](@ref) for each system. The default value of `nothing`
leads to the default move probabilities: [98% translation, 2% random translation] for a
monoatomic species, or [49% translation, 49% rotation, 2% random reinsertion] else.

`rng` is the random number generator, defaulting to `Random.default_rng()`.
"""
function setup_montecarlo(framework, pff, systems;
                          blockfiles=fill(nothing, length(systems)), gridstep=0.15u"Å",
                          supercell=nothing, new=false, restart=nothing, parallel=true,
                          mcmoves=fill(nothing, length(systems)), rng=default_rng(),
                          cutoff=12.0u"Å")
    if framework isa AbstractString && endswith(framework, ".cif")
        framework = framework[1:prevind(framework, end-3)]
    end
    syst_framework = framework isa AbstractSystem ? framework : load_framework_RASPA(framework, pff)
    ff = _ff(pff; cutoff)
    mat = stack3(bounding_box(syst_framework))
    if supercell isa Nothing
        supercell = find_supercell(syst_framework, cutoff)
    end
    supercell::NTuple{3,Int}
    cell = CellMatrix(SMatrix{3,3,TÅ,9}(stack(bounding_box(syst_framework).*supercell)))

    csetup = GridCoordinatesSetup(syst_framework, gridstep)
    length(blockfiles) == length(systems) || error("Please provide one blockfiles element per system")
    blocksetup = [decide_parse_block(file, s isa Tuple ? s[1] : s, framework) for (file, s) in zip(blockfiles, systems)]

    needcoulomb = false
    encountered_atoms = Set{Symbol}()
    for s in systems
        system = default_system((s isa Tuple ? s[1] : s), pff)
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

    coulomb_grid_path, vdw_grid_paths = grid_locations(framework, pff, ff, first.(atoms), gridstep, supercell)

    coulomb, eframework = if needcoulomb
        _eframework = initialize_ewald(syst_framework)
        retrieve_or_create_grid(coulomb_grid_path, syst_framework, ff, gridstep, _eframework, mat, new, cutoff), _eframework
    else
        EnergyGrid(true), EwaldFramework(mat)
    end

    grids = if isempty(vdw_grid_paths) || isempty(framework)
        # if framework is empty, signal it by having an empty list of grids
        EnergyGrid[]
    else
        _grids = Vector{EnergyGrid}(undef, length(ff.sdict))
        for (j, (atom, i)) in enumerate(atoms)
            _grids[i] = retrieve_or_create_grid(vdw_grid_paths[j], syst_framework, ff, gridstep, atom, mat, new, cutoff)
        end
        _grids
    end

    num_framework_atoms = zeros(Int, length(ff.sdict))
    Π = prod(supercell)
    for at in syst_framework
        num_framework_atoms[ff.sdict[Symbol(get_atom_name(atomic_symbol(at)))]] += Π
    end

    restartpositions = restart isa Nothing ? nothing : read_restart_RASPA(restart)

    setup_montecarlo(cell, csetup, systems, ff, eframework, coulomb, grids, blocksetup, num_framework_atoms, restartpositions; parallel, rng, mcmoves, pff)
end

"""
    MonteCarloSetup(mc::MonteCarloSetup, o::SimulationStep=mc.step; parallel::Bool=mc.step.parallel, mcmoves::AbstractVector{MCMoves}=copy(mc.mcmoves))

Create a copy of `mc` that does not share its modifiable internal states (the positions and
the Ewald state). For example, the copy and the original can be used to run Monte-Carlo
simulations in parallel, and the state of one will not influence that of the other.

The step elements of `mc` are actually copied from `o`, which can be set to a different
`SimulationStep` if need be.

`parallel` specifies whether the computations on the resulting `MonteCarloSetup` should be
parallelized or not.

`mcmoves` specifies the new [`MCMoves`](@ref) of the setup.

!!! warning
    Internal states that are semantically immutable are shared, although some of them are
    technically mutable, like the value of the atom charges for instance. Modifying such a
    state on the original or the copy can thus also modify that of the other: use
    `deppcopy` or a custom copy implementation to circumvent this issue if you plan on
    modifying states outside of the API.
"""
function MonteCarloSetup(mc::MonteCarloSetup, o::SimulationStep=mc.step;
                         parallel::Bool=mc.step.parallel,
                         mcmoves::AbstractVector{MCMoves}=copy(mc.mcmoves),
                         rng = mc.rng isa TaskLocalRNG ? mc.rng : copy(mc.rng))
    ewaldsystems = Vector{EwaldSystem}[]
    for (i, posidxi) in enumerate(o.posidx)
        ffidxi = o.ffidx[i]
        charge = [o.charges[k] for k in ffidxi]
        push!(ewaldsystems, [EwaldSystem(o.positions[molpos], charge) for molpos in posidxi])
    end
    emptysystems = findall(isempty, ewaldsystems)
    for i_empty in emptysystems
        push!(ewaldsystems[i_empty], EwaldSystem(mc.models[i_empty], o.charges[o.ffidx[i_empty]]))
    end
    ewald = IncrementalEwaldContext(EwaldContext(mc.ewald.ctx.eframework, ewaldsystems, emptysystems))
    MonteCarloSetup(SimulationStep(o, :all; parallel),
                    ewald, map(copy, mc.flatidx), copy(mc.revflatidx), copy(mc.tailcorrection),
                    mc.coulomb, mc.grids, mc.speciesblocks, mc.atomblocks, mc.bead, mc.models,
                    mcmoves, mc.gcmcdata, rng)
end


struct NanoSystem # pseudo-AbstractSystem with only positions and atom symbols
    position::Vector{SVector{3,TÅ}}
    atomic_symbols::Vector{Symbol}
end
AtomsBase.position(x::NanoSystem) = x.position
AtomsBase.atomic_symbol(x::NanoSystem, i) = x.atomic_symbols[i]
Base.length(x::NanoSystem) = length(x.position)
function NanoSystem(mc::MonteCarloSetup, i::Int)
    position = mc.models[i]
    atomic_symbols = [mc.step.ff.symbols[x] for x in mc.step.ffidx[i]]
    NanoSystem(position, atomic_symbols)
end

"""
    move_one_system!(mc::MonteCarloSetup, idx, positions)

Move the species given by its index `idx` in `mc` to `positions`.
"""
function move_one_system!(mc::MonteCarloSetup, idx, positions)
    if isdefined_ewald(mc) && mc.ewald.last[] != -1
        i, j = idx
        flatidxi = mc.flatidx[i]
        ij = j == length(flatidxi) + 1 ? -i : mc.flatidx[i][j]
        single_contribution_ewald(mc.ewald, ij, positions)
    end
    update_mc!(mc, idx, positions)
end

struct FrameworkEnergyReport
    vdw::TK
    direct::TK
end
FrameworkEnergyReport() = FrameworkEnergyReport(0.0u"K", 0.0u"K")
Base.Number(f::FrameworkEnergyReport) = f.vdw + f.direct
Base.Float64(f::FrameworkEnergyReport) = ustrip(u"K", Number(f))::Float64
function Base.isapprox(f1::FrameworkEnergyReport, f2::FrameworkEnergyReport; atol::Real=0, rtol::Real=(atol>0 ? 0 : sqrt(eps())))
    isapprox(ustrip(u"K", f1.vdw), ustrip(u"K", f2.vdw); atol, rtol) &&
    isapprox(ustrip(u"K", f1.direct), ustrip(u"K", f2.direct); atol, rtol)
end
function Base.isapprox(f::FrameworkEnergyReport, x::Number; atol::Real=0, rtol::Real=(atol>0 ? 0 : sqrt(eps())))
    isapprox(Float64(f), x; atol, rtol)
end
Base.:(-)(f::FrameworkEnergyReport) = FrameworkEnergyReport(-f.vdw, -f.direct)
Base.:(/)(f::FrameworkEnergyReport, n::Int) = FrameworkEnergyReport(f.vdw/n, f.direct/n)

struct MCEnergyReport
    framework::FrameworkEnergyReport
    inter::TK
    reciprocal::TK
end
MCEnergyReport(v::TK, d::TK, i::TK, r::TK) = MCEnergyReport(FrameworkEnergyReport(v, d), i, r)
MCEnergyReport(v::T, d::T, i::T, r::T) where {T<:Real} = MCEnergyReport(v*u"K", d*u"K", i*u"K", r*u"K")
Base.Number(e::MCEnergyReport) = Number(e.framework) + e.inter + e.reciprocal
Base.Float64(e::MCEnergyReport) = ustrip(u"K", Number(e))::Float64
Base.show(io::IO, e::MCEnergyReport) = show(io, Float64(e)*u"K")
function Base.show(io::IO, ::MIME"text/plain", e::MCEnergyReport)
    println(io, Float64(e), " = [", e.framework.vdw, " (f VdW) + ", e.framework.direct, " (f direct)] + [", e.inter, " (internal) + ", e.reciprocal, " (reciprocal)]")
end
function Base.isapprox(e1::MCEnergyReport, e2::MCEnergyReport; atol::Real=0, rtol::Real=(atol>0 ? 0 : sqrt(eps())))
    isapprox(e1.framework, e2.framework; atol, rtol) &&
    isapprox(ustrip(u"K", e1.inter), ustrip(u"K", e2.inter); atol, rtol) &&
    isapprox(ustrip(u"K", e1.reciprocal), ustrip(u"K", e2.reciprocal); atol, rtol)
end
function Base.isapprox(e::MCEnergyReport, x::Number; atol::Real=0, rtol::Real=(atol>0 ? 0 : sqrt(eps())))
    isapprox(Float64(e), x; atol, rtol)
end
Base.:(-)(e::MCEnergyReport) = MCEnergyReport(-e.framework, -e.inter, -e.reciprocal)
Base.:(/)(e::MCEnergyReport, n::Int) = MCEnergyReport(e.framework/n, e.inter/n, e.reciprocal/n)

struct BaselineEnergyReport
    er::MCEnergyReport
    tailcorrection::TK
end
BaselineEnergyReport(f::FrameworkEnergyReport, i::TK, r::TK, t::TK) = BaselineEnergyReport(MCEnergyReport(f, i, r), t)
BaselineEnergyReport(v::TK, d::TK, i::TK, r::TK, t::TK) = BaselineEnergyReport(MCEnergyReport(v, d, i, r), t)
BaselineEnergyReport(v::T, d::T, i::T, r::T, t::T) where {T<:Real} = BaselineEnergyReport(MCEnergyReport(v, d, i, r), t*u"K")
Base.Number(ber::BaselineEnergyReport) = Number(ber.er) + ber.tailcorrection
Base.Float64(ber::BaselineEnergyReport) = ustrip(u"K", Number(ber))::Float64
Base.show(io::IO, ber::BaselineEnergyReport) = show(io, Float64(ber)*u"K")
function Base.show(io::IO, ::MIME"text/plain", b::BaselineEnergyReport)
    println(io, Float64(b), " = [", b.er.framework.vdw, " (f VdW) + ", b.er.framework.direct, " (f direct)] + [", b.er.inter, " (internal) + ", b.er.reciprocal, " (reciprocal)] + ", b.tailcorrection, " (tailcorrection)")
end
function Base.isapprox(b1::BaselineEnergyReport, b2::BaselineEnergyReport; atol::Real=0, rtol::Real=(atol>0 ? 0 : sqrt(eps())))
    isapprox(b1.er, b2.er; atol, rtol) &&
    isapprox(ustrip(u"K", b1.tailcorrection), ustrip(u"K", b2.tailcorrection); atol, rtol)
end
function Base.isapprox(b::BaselineEnergyReport, x::Number; atol::Real=0, rtol::Real=(atol>0 ? 0 : sqrt(eps())))
    isapprox(Float64(b), x; atol, rtol)
end
Base.:(-)(b::BaselineEnergyReport) = BaselineEnergyReport(-b.er, -b.tailcorrection)
Base.:(/)(b::BaselineEnergyReport, n::Int) = BaselineEnergyReport(b.er/n, b.tailcorrection/n)

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
        function Base.$(op)(b1::BaselineEnergyReport, b2::BaselineEnergyReport)
            BaselineEnergyReport($op(b1.er, b2.er), $op(b1.tailcorrection, b2.tailcorrection))
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
    baseline_energy(hs::SiteHopping)

Compute the energy of the current configuration.
"""
baseline_energy

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
`movement_energy(mc, (i,j), positions) - movement_energy(mc, (i,j))`. See also
[`combined_movement_energy`](@ref) for a slightly more optimized version of that
computation.

If `j` is a non-attributed molecule number for kind `i`, `positions` must be provided. Note
that, in that case, the tail correction contribution will not be included.

!!! warning
    `baseline_energy(mc)` must have been called at least once before, otherwise the
    computation of the Ewald part will error.
    See also [`single_contribution_ewald`](@ref).
"""
function movement_energy(mc::MonteCarloSetup, idx, positions=nothing)
    i, j = idx
    flatidxi = mc.flatidx[i]
    ij = j == length(flatidxi) + 1 ? -i : mc.flatidx[i][j]
    poss = if positions isa Nothing
        molpos = mc.step.posidx[i][j]
        @view mc.step.positions[molpos]
    else
        positions
    end
    @spawnif mc.step.parallel begin
        singlereciprocal = @spawn single_contribution_ewald($mc.ewald, $ij, $positions)
        fer = @spawn framework_interactions($mc, $i, $poss)
        singlevdw = single_contribution_vdw(mc.step, idx, poss)
        MCEnergyReport(fetch(fer)::FrameworkEnergyReport, singlevdw, fetch(singlereciprocal)::TK)
    end
end

"""
    combined_movement_energy(mc, idx, newpos)

Equivalent to `(movement_energy(mc, idx, nothing), movement_energy(mc, idx, newpos))` but
slightly more parallel-efficient.

See also [`movement_energy`](@ref) for more information.
"""
function combined_movement_energy(mc, idx, newpos)
    i, j = idx
    ij = mc.flatidx[i][j]
    molpos = mc.step.posidx[i][j]
    positions = @view mc.step.positions[molpos]
    @spawnif mc.step.parallel begin
        singlevdw_before_task = @spawn single_contribution_vdw($mc.step, $(i,j), $positions)
        singlereciprocal_after = @spawn single_contribution_ewald($mc.ewald, $ij, $newpos)
        singlereciprocal_before = @spawn single_contribution_ewald($mc.ewald, $ij, nothing)
        fer_before = @spawn framework_interactions($mc, $i, $positions)
        fer_after = @spawn framework_interactions($mc, $i, $newpos)
        singlevdw_before = fetch(singlevdw_before_task)::TK
        singlevdw_after = single_contribution_vdw(mc.step, (i,j), newpos)
        before = MCEnergyReport(fetch(fer_before)::FrameworkEnergyReport, singlevdw_before, fetch(singlereciprocal_before)::TK)
        after = MCEnergyReport(fetch(fer_after)::FrameworkEnergyReport, fetch(singlevdw_after)::TK, fetch(singlereciprocal_after)::TK)
    end
    before, after
end


"""
    update_mc!(mc::MonteCarloSetup, idx::Tuple{Int,Int}, positions::Vector{SVector{3,TÅ}}, [swapinfo::SwapInformation])

Following a call to [`movement_energy(mc, idx, positions)`](@ref), update the internal
state of `mc` so that the species of index `idx` is now at `positions`.
"""
function update_mc!(mc::MonteCarloSetup, (i,j)::Tuple{Int,Int}, positions::Vector{SVector{3,TÅ}}, swapinfo::SwapInformation=SwapInformation(0.0u"K", 0, false, false))
    if swapinfo.isswap
        if swapinfo.isinsertion
            add_one_system!(mc, i, positions) # TODO: optimize by making mc aware that the energy computation has already been done
        else
            remove_one_system!(mc, i, j) # TODO: same optimization as for insertion
        end
    else
        update_ewald_context!(mc.ewald)
        L = mc.step.posidx[i][j]
        mc.step.positions[L] .= positions
    end
    nothing
end


function inblockpocket(block::BlockFile, atomblocks::Vector{BlockFile}, ffidxi::Vector{Int}, newpos::AbstractVector{SVector{3,TÅ}})
    for (j, pos) in enumerate(newpos)
        block[pos] && return true
        if !isempty(atomblocks) # atomblocks
            ablock = atomblocks[ffidxi[j]]
            ablock[pos + ablock.csetup.Δ./2] && return true
        end
    end
    return false
end


"""
    energy_point_mc!(mc::MonteCarloSetup, positions, current, before)

Compute the energy associated with the last molecule of kind 1 in `mc` at the given
`positions`.

`current` is the energy of `mc` in its current state. `before` is the contribution to the
energy of the moving molecule.

Return the energy and new values for `current` and `before`, up to date with the new
version of `mc`. In many cases, neither `mc`, `current` nor `before` are modified.

If `current` is above `1e50u"K"`, the energy will be computed with [`baseline_energy`](@ref).

!!! warning
    This function modifies `mc`. In particular, it is assumed that `mc` is task-local, i.e.
    `mc` must not be modified concurrently.
"""
function energy_point_mc!(mc::MonteCarloSetup, positions, current, before)
    if inblockpocket(mc.speciesblocks[1], mc.atomblocks, first(mc.step.ffidx), positions)
        return 1e100u"K", current, before
    end
    # current = Inf*u"K" # FIXME: remove this and uncomment the previous lines after understanding and fixing why the results differ
    j = length(mc.flatidx[1])
    ij = mc.flatidx[1][j]
    if mc.ewald.last[] == -1 || current ≥ 1e50*u"K"
        @label compute_baseline
        molpos = mc.step.posidx[1][j]
        mc.step.positions[molpos] = positions
        isdefined_ewald(mc) && move_one_system!(mc.ewald.ctx, ij, positions)
        ret = Number(baseline_energy(mc))
        if ret < Inf*u"K"
            oldpos = @view mc.step.positions[molpos]
            fer_before = framework_interactions(mc, 1, oldpos)
            singlevdw_before = single_contribution_vdw(mc.step, (1,j), oldpos)
            singlereciprocal_before = single_contribution_ewald(mc.ewald, ij, nothing)
            newbefore = Number(MCEnergyReport(fer_before, singlevdw_before, singlereciprocal_before))
            return ret, ret, newbefore
        end
        return ret, Inf*u"K", NaN*u"K"
    end
    @assert !isnan(before)
    fer = framework_interactions(mc, 1, positions)
    singlereciprocal = single_contribution_ewald(mc.ewald, ij, positions)
    singlevdw = single_contribution_vdw(mc.step, (1,j), positions)
    after = Number(MCEnergyReport(fer, singlevdw, singlereciprocal))
    result = current + after - before
    if result < current && (abs(result) < max(abs(current), abs(before))/10)
        @goto compute_baseline
    end
    if result < (current < 0.0u"K" ? 2*current : current/2)
        update_mc!(mc, (1,j), positions)
        return result, result, after
    end
    return result, current, before
end


choose_random_species(mc::MonteCarloSetup) = rand(mc.rng, mc.revflatidx)
function compute_accept_move(before, after, T, mc, swapinfo::Union{SwapInformation,Nothing}=nothing)
    b = Number(before)
    a = Number(after)
    if swapinfo isa SwapInformation && swapinfo.isswap
        tcchange = modify_species_dryrun(mc.tailcorrection, swapinfo.i, ifelse(swapinfo.isinsertion, 1, -1))
        return compute_accept_move_swap(a+tcchange-b, T, mc, swapinfo)
    end
    a < b && return true
    e = exp((b-a)/T)
    return rand(mc.rng) < e
end


function randomize_position!(positions, rng, indices, bead, block, ffidxi, atomblocks, mat)
    pos = @view positions[indices]
    for _ in 1:30
        posr = random_rotation(rng, random_rotation(rng, random_rotation(rng, pos, 90u"°", bead, 1), 90u"°", bead, 2), 90u"°", bead, 3)
        for _ in 1:1000
            post = random_translation(rng, posr, mat)
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

!!! warning
    Using `update_ewald=false` will result in an inconsistent state for `mc` which will
    error on [`run_montecarlo`](@ref) unless [`baseline_energy(mc)`](@ref) is called.

!!! warning
    `update_ewald=true` requires a prior call to [`baseline_energy(mc)`](@ref).
"""
function randomize_position!(mc::MonteCarloSetup, (i,j), update_ewald=true)
    newpos = randomize_position!(mc.step.positions, mc.step.posidx[i][j], mc.bead[i], mc.speciesblocks[i], mc.step.ffidx[i], mc.atomblocks, mc.step.mat)
    if update_ewald
        single_contribution_ewald(mc.ewald, mc.offsets[i] + j, newpos)
        update_mc!(mc, (i,j), newpos)
    else
        mc.ewald.last[] = -1 # signal the inconsistent state
    end
    newpos
end
=#

"""
    remove_one_system!(mc::MonteCarloSetup, i::Int, j::Int)

Remove the `j`-th molecule of kind i.

Return an index `j′` such that the molecule which used to be the `j′`-th molecule of kind
`i` is now the `j`-th molecule of kind `i`. `j′ == j` if there is no modification.
"""
function remove_one_system!(mc::MonteCarloSetup, i::Int, j::Int)
    # 0. Sanity check
    hasewald = !isempty(mc.ewald.ctx.charges)
    hasewald && abs(sum(mc.ewald.ctx.charges[i])) ≥ 1e-8u"e_au" && error("Cannot remove a charged species.")

    # 1. Update step
    # see comment on the method in remove_one_system! in ewald.jl
    oldn = length(mc.step.atoms)
    prev_atoms = sort!(mc.step.posidx[i][j])
    n = length(prev_atoms) # number of atoms of the molecule to delete
    lastidx = length(prev_atoms)
    firstidx = 1
    for k in 1:n
        l = oldn + 1 - k
        if l == prev_atoms[lastidx]
            lastidx -= 1
        else
            first_l = prev_atoms[firstidx]
            mc.step.positions[first_l] = mc.step.positions[l]
            new_i, new_j, new_k = mc.step.atoms[first_l] = mc.step.atoms[l]
            mc.step.posidx[new_i][new_j][new_k] = first_l
            firstidx += 1
        end
    end
    resize!(mc.step.positions, length(mc.step.positions) - n)
    resize!(mc.step.atoms, length(mc.step.atoms) - n)

    # 2. Ewald
    ewij = mc.flatidx[i][j]
    oldewij = hasewald ? remove_one_system!(mc.ewald, ewij) : length(mc.revflatidx)
    oldewi, oldewj = mc.revflatidx[oldewij]
    mc.flatidx[oldewi][oldewj] = ewij
    mc.revflatidx[ewij] = (oldewi, oldewj)

    # 3. Swap the j-th and lastj-th molecules and shrink the vectors
    lastj = length(mc.step.posidx[i])
    @assert lastj == length(mc.flatidx[i])
    if j != lastj
        lastij = mc.flatidx[i][lastj]
        mc.flatidx[i][j] = lastij
        mc.revflatidx[lastij] = (i,j)
        mc.step.posidx[i][j] = mc.step.posidx[i][lastj]
        for (k, l) in enumerate(mc.step.posidx[i][j])
            mc.step.atoms[l] = (i, j, k)
        end
    end
    pop!(mc.flatidx[i])
    pop!(mc.revflatidx)
    pop!(mc.step.posidx[i])

    # 4. Tail correction
    modify_species!(mc.tailcorrection, i, -1)

    lastj
end

"""
    add_one_system!(mc::MonteCarloSetup, i::Int, positions=nothing; uninitialize=isnothing(positions))

Add a new molecule of kind `i`. Return its index `j`.

If `positions` is unset, the atom positions will be that of the model molecule used to
create `mc`. Otherwise, these specify the new atomic positions.

If `uninitialize` is set, then the Ewald state is not updated and a call to
`baseline_energy` will be required to obtain correct energies. Otherwise, the Ewald state
is maintained up to date if it had already been initialized.
"""
function add_one_system!(mc::MonteCarloSetup, i::Int, positions=nothing; uninitialize=false)
    # 0. Sanity check
    hasewald = !isempty(mc.ewald.ctx.charges)
    hasewald && abs(sum(mc.ewald.ctx.charges[i])) ≥ 1e-8u"e_au" && error("Cannot add a charged species.")

    # 1. Update step
    oldn = length(mc.step.atoms)
    n = length(mc.step.ffidx[i]) # number of atoms of the molecule to add
    @assert n == length(mc.models[i])
    push!(mc.step.posidx[i], collect((oldn+1):(oldn+n)))
    j = length(mc.step.posidx[i])
    append!(mc.step.atoms, (i, j, k) for k in 1:n)
    poss = positions isa Nothing ? mc.models[i] : positions
    append!(mc.step.positions, poss)

    # 2. Ewald
    ij = length(mc.revflatidx)+1
    if hasewald
        if uninitialize
            mc.ewald.last[] = -1 # reset the Ewald state
        end
        ewij = add_one_system!(mc.ewald, i, poss)
        @assert ewij == ij
    end
    push!(mc.flatidx[i], ij)
    push!(mc.revflatidx, (i,j))

    # 3. Tail correction
    modify_species!(mc.tailcorrection, i, 1)

    j
end
