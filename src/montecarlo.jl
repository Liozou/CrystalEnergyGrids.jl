export setup_montecarlo, baseline_energy, movement_energy

struct MonteCarloSimulation
    ff::ForceField
    ewald::IncrementalEwaldContext
    cell::CellMatrix
    tailcorrection::Base.RefValue{typeof(1.0u"K")}
    coulomb::EnergyGrid
    grids::Vector{EnergyGrid} # grids[i] is the VdW grid for atom i in ff.
    offsets::Vector{Int}
    # offsets[i] is the number of molecules belonging to a kind strictly lower than i.
    indexof::Vector{Tuple{Int,Int}}
    # indexof[m] is the index of the m-th species of the system. In particular, if there
    # are at least `j` species of kind `i`, then indexof[offsets[i]+j] == (i,j)
    idx::Vector{Vector{Int}}
    # idx[i][k] is ff.sdict[atomic_symbol(systemkinds[i], k)], i.e. the numeric index
    # in the force field for the k-th atom in a system of kind i.
    charges::Vector{typeof(1.0u"e_au")}
    # charges[ix] is the charge of the atom whose of index ix in ff, i.e. the charge of the
    # k-th atom in a system of kind i is charges[idx[i][k]].
    positions::Vector{Vector{Vector{SVector{3,typeof(1.0u"Å")}}}}
    # positions[i][j][k] is the position of the k-th atom in the j-th system of kind i.
    isrigid::BitVector # always set to true
    speciesblocks::Vector{BlockFile}
    # speciesblocks[i][pos] is set if pos is blocked for any atom of a species of kind i.
    atomblocks::Vector{BlockFile}
    # atomblock[ix][pos] is set if pos is blocked for all atoms of index ix in ff.
    bead::Vector{Int} # k = bead[i] is the number of the reference bead of kind i.
end

function SimulationStep(mc::MonteCarloSimulation)
    SimulationStep(mc.ff, mc.charges, mc.positions, mc.isrigid, mc.idx, mc.cell)
end


struct EwaldSystem # pseudo-AbstractSystem with only positions and charges
    position::Vector{SVector{3,typeof(1.0u"Å")}}
    atomic_charge::Vector{typeof(1.0u"e_au")}
end
AtomsBase.position(x::EwaldSystem) = x.position
AtomsBase.position(x::EwaldSystem, i::Int) = x.position[i]
Base.getindex(x::EwaldSystem, ::Colon, s::Symbol) = getproperty(x, s)
Base.getindex(x::EwaldSystem, i::Int, s::Symbol) = getproperty(x, s)[i]
Base.length(x::EwaldSystem) = length(x.position)

struct IdSystem # pseudo-AbstractSystem with only atomic symbols and charges
    atomic_symbol::Vector{Symbol}
    atomic_charge::Vector{typeof(1.0u"e_au")}
end
IdSystem(s::AbstractSystem{3}) = IdSystem(atomic_symbol(s), s[:,:atomic_charge])
Base.length(s::IdSystem) = length(s.atomic_charge)

function setup_montecarlo(cell::CellMatrix, csetup::GridCoordinatesSetup,
                          systems, ff::ForceField, eframework::EwaldFramework,
                          coulomb::EnergyGrid, grids::Vector{EnergyGrid},
                          blocksetup::Vector{String},
                          num_framework_atoms::Vector{Int})
    if any(≤(24.0u"Å"), perpendicular_lengths(cell.mat))
        error("The current cell has at least one perpendicular length lower than 24.0Å: please use a larger supercell")
    end

    kindsdict = Dict{Tuple{Vector{Symbol},String},Int}()
    systemkinds = IdSystem[]
    U = Vector{SVector{3,typeof(1.0u"Å")}} # positions of the atoms of a system
    poss = Vector{U}[]
    indices = Tuple{Int,Int}[]
    rev_indices = Vector{Int}[]
    speciesblocks = BlockFile[]
    for (i, (s, block)) in enumerate(zip(systems, blocksetup))
        system, n = s isa Tuple ? s : (s, 1)
        m = length(kindsdict)+1
        kind = get!(kindsdict, (atomic_symbol(system)::Vector{Symbol}, block), m)
        if kind === m
            push!(systemkinds, IdSystem(system)::IdSystem)
            push!(poss, U[])
            push!(rev_indices, Int[])
            if isempty(block)
                push!(speciesblocks, BlockFile(csetup))
            else
                push!(speciesblocks, parse_blockfile(joinpath(getdir_RASPA(), "structures", "block", block)*".block", csetup))
            end
        end
        push!(indices, (kind, length(poss[kind])+1))
        append!(rev_indices[kind], i for _ in 1:n)
        append!(poss[kind], copy(position(system)::Vector{SVector{3,typeof(1.0u"Å")}}) for _ in 1:n)
    end

    idx = [[ff.sdict[s.atomic_symbol[k]] for k in 1:length(s)] for s in systemkinds]
    charges = fill(NaN*u"e_au", length(ff.sdict))
    for (i, ids) in enumerate(idx)
        kindi = systemkinds[i]
        for (k, ix) in enumerate(ids)
            charges[ix] = kindi.atomic_charge[k]
        end
    end

    num_atoms = copy(num_framework_atoms)
    for (i, indices) in enumerate(idx)
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
    tcorrection *= 2π/det(ustrip.(cell.mat))

    n = length(poss)
    offsets = Vector{Int}(undef, n)
    offsets[1] = 0
    indexof = Tuple{Int,Int}[]
    for i in 1:(n-1)
        numi = length(poss[i])
        offsets[i+1] = offsets[i] + numi
        append!(indexof, (i,j) for j in 1:length(poss[i]))
    end
    append!(indexof, (n,j) for j in 1:length(poss[n]))

    buffer = MVector{3,typeof(1.0u"Å")}(undef)
    buffer2 = MVector{3,Float64}(undef)
    beads = Vector{Int}(undef, n)
    for i in 1:n
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
        idxi = idx[i]
        charge = [charges[k] for k in idxi]
        rev_idx = rev_indices[i]
        bead = beads[i]
        block = speciesblocks[i]
        d = norm(sum(cell.mat; dims=2))/4
        for j in 1:length(positioni)
            s = systems[rev_idx[j]]
            n = s isa Tuple ? s[2]::Int : 1
            if n > 1
                randomize_position!(positioni, j, bead, block, idxi, atomblocks, d)
            end
            push!(ewaldsystems, EwaldSystem(positioni[j], charge))
        end
    end

    ewald = IncrementalEwaldContext(EwaldContext(eframework, ewaldsystems))

    MonteCarloSimulation(ff, ewald, cell, Ref(tcorrection), coulomb, grids, offsets, indexof,
                         idx, charges, poss, trues(length(idx)), speciesblocks, atomblocks, beads), indices
end

"""
    setup_montecarlo(framework, forcefield_framework::String, systems;
                     blockfiles=fill(nothing, length(systems)), gridstep=0.15u"Å", supercell=nothing, new=false)

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
"""
function setup_montecarlo(framework, forcefield_framework::String, systems;
                          blockfiles=fill(nothing, length(systems)), gridstep=0.15u"Å", supercell=nothing, new=false)
    syst_framework = load_framework_RASPA(framework, forcefield_framework)
    ff = parse_forcefield_RASPA(forcefield_framework)
    mat = stack3(bounding_box(syst_framework))
    if supercell isa Nothing
        supercell = find_supercell(syst_framework, 12.0u"Å")
    end
    supercell::NTuple{3,Int}
    cell = CellMatrix(SMatrix{3,3,typeof(1.0u"Å"),9}(stack(bounding_box(syst_framework).*supercell)))

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

    setup_montecarlo(cell, csetup, systems, ff, eframework, coulomb, grids, blocksetup, num_framework_atoms)
end

"""
    MonteCarloSimulation(mc::MonteCarloSimulation)

Create a copy of `mc` that does not share its modifiable internal states (the positions and
the Ewald state). For example, the copy and the original can be used to run Monte-Carlo
simulations in parallel, and the state of one will not influence that of the other.

!!! warn
    Internal states that are semantically immutable are shared, although some of them are
    technically mutable, like the value of the atom charges for instance. Modifying such a
    state on the original or the copy can thus also modify that of the other: use
    `deppcopy` or a custom copy implementation to circumvent this issue if you plan on
    modifying states outside of the API.
"""
function MonteCarloSimulation(mc::MonteCarloSimulation)
    positions = [[[pos for pos in poss] for poss in positioni] for positioni in mc.positions]
    ewaldsystems = EwaldSystem[]
    for (i, positioni) in enumerate(positions)
        idxi = mc.idx[i]
        charge = [mc.charges[k] for k in idxi]
        append!(ewaldsystems, EwaldSystem(poss, charge) for poss in positioni)
    end
    ewald = IncrementalEwaldContext(EwaldContext(mc.ewald.ctx.eframework, ewaldsystems))
    MonteCarloSimulation(mc.ff, ewald, mc.cell, Ref(mc.tailcorrection[]), mc.coulomb,
                         mc.grids, mc.offsets, mc.indexof, mc.idx, mc.charges, positions,
                         mc.isrigid, mc.speciesblocks, mc.atomblocks, mc.bead)
end

function set_position!(mc::MonteCarloSimulation, (i, j), newpositions, newEiks=nothing)
    mc.positions[i][j] = if eltype(positions) <: AbstractVector{<:AbstractFloat}
        newpositions
    else
        (NoUnits.(newpositions[j]/u"Å") for j in 1:length(newpositions))
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
    vdw::typeof(1.0u"K")
    direct::typeof(1.0u"K")
end
FrameworkEnergyReport() = FrameworkEnergyReport(0.0u"K", 0.0u"K")
Base.Float64(f::FrameworkEnergyReport) = Float64((f.vdw + f.direct)/u"K")

struct MCEnergyReport
    framework::FrameworkEnergyReport
    inter::typeof(1.0u"K")
    reciprocal::typeof(1.0u"K")
end
Base.Float64(e::MCEnergyReport) = Float64(e.framework) + Float64((e.inter + e.reciprocal)/u"K")
Base.show(io::IO, e::MCEnergyReport) = show(io, Float64(e)*u"K")
function Base.show(io::IO, ::MIME"text/plain", e::MCEnergyReport)
    println(io, Float64(e), " = [", e.framework.vdw, " (f VdW) + ", e.framework.direct, " (f direct)] + [", e.inter, " (internal) + ", e.reciprocal, " (reciprocal)]")
end

struct BaselineEnergyReport
    er::MCEnergyReport
    tailcorrection::typeof(1.0u"K")
end
BaselineEnergyReport(f, i, r, t) = BaselineEnergyReport(MCEnergyReport(f, i, r), t)
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


function framework_interactions(grids::Vector{EnergyGrid}, coulombgrid::EnergyGrid, charges::Vector{typeof(1.0u"e_au")}, indices::Vector{Int}, positions)
    isempty(grids) && return FrameworkEnergyReport()
    n = length(indices)
    vdw = 0.0u"K"
    direct = 0.0u"K"
    hascoulomb = coulombgrid.ewald_precision != -Inf
    for k in 1:n
        ix = indices[k]
        pos = positions[k]
        vdw += interpolate_grid(grids[ix], pos)
        if hascoulomb
            coulomb = interpolate_grid(coulombgrid, pos)
            direct += ifelse(coulomb == 1e100u"K", coulomb, Float64(charges[ix]/u"e_au")*coulomb)
        end
    end
    return FrameworkEnergyReport(vdw, direct)
end
function framework_interactions(mc::MonteCarloSimulation, indices::Vector{Int}, positions)
    framework_interactions(mc.grids, mc.coulomb, mc.charges, indices, positions)
end

"""
    framework_interactions(mc::MonteCarloSimulation, i, positions)

Energy contribution of the interaction between the framework and a molecule of system kind
`i` at the given `positions`.

Only return the Van der Waals and the direct part of the Ewald summation. The reciprocal
part can be obtained with [`compute_ewald`](@ref).
"""
function framework_interactions(mc::MonteCarloSimulation, i::Int, positions)
    framework_interactions(mc, mc.idx[i], positions)
end

"""
    baseline_energy(mc::MonteCarloSimulation)

Compute the energy of the current configuration.
"""
function baseline_energy(mc::MonteCarloSimulation)
    reciprocal = compute_ewald(mc.ewald)
    vdw = compute_vdw(SimulationStep(mc))
    fer = FrameworkEnergyReport()
    isempty(mc.grids) && return BaselineEnergyReport(fer, vdw, reciprocal, mc.tailcorrection[])
    for (i, indices) in enumerate(mc.idx)
        poss_i = mc.positions[i]
        for poss in poss_i
            fer += framework_interactions(mc, indices, poss)
        end
    end
    return BaselineEnergyReport(fer, vdw, reciprocal, mc.tailcorrection[])
end

"""
    movement_energy(mc::MonteCarloSimulation, (i, j), positions=nothing)

Compute the energy contribution of the `j`-th molecule of kind `i` when placed at
`positions`. If not provided, `positions` is the current position for that molecule.

The energy difference between the new position for the molecule and the current one is
`movement_energy(mc, (i,j), positions) - movement_energy(mc, (i,j))`.

!!! warning
    `baseline_energy(mc)` must have been called at least once before, otherwise the
    computation of the Ewald part will error.
    See also [`single_contribution_ewald`](@ref).
"""
function movement_energy(mc::MonteCarloSimulation, idx, positions=nothing)
    i, j = idx
    k = mc.offsets[i]+j
    poss = positions isa Nothing ? mc.positions[idx[1]][idx[2]] : positions
    singlevdw = @spawn single_contribution_vdw(SimulationStep(mc), (i,j), poss)
    fer = @spawn framework_interactions(mc, i, poss)
    singlereciprocal = single_contribution_ewald(mc.ewald, k, positions)
    MCEnergyReport(fetch(fer), fetch(singlevdw), singlereciprocal)
end

"""
    update_mc!(mc::MonteCarloSimulation, idx, positions)

Following a call to [`movement_energy(mc, idx, positions)`](@ref), update the internal
state of `mc` so that the species of index `idx` is now at `positions`.
"""
function update_mc!(mc::MonteCarloSimulation, idx, positions)
    update_ewald_context!(mc.ewald)
    mc.positions[idx[1]][idx[2]] = positions
    nothing
end


function inblockpocket(block::BlockFile, atomblocks::Vector{BlockFile}, idx::Vector{Int}, newpos::Vector{SVector{3,typeof(1.0u"Å")}})
    for (j, pos) in enumerate(newpos)
        block[pos] && return true
        ablock = atomblocks[idx[j]]
        ablock[pos + ablock.csetup.Δ./2] && return true
    end
    return false
end

function random_translation(positions::Vector{SVector{3,typeof(1.0u"Å")}}, dmax::typeof(1.0u"Å"))
    r = SVector{3}(((2*rand()-1)*dmax) for _ in 1:3)
    [poss + r for poss in positions]
end
function random_rotation(positions::Vector{SVector{3,typeof(1.0u"Å")}}, θmax, bead, _r=nothing)
    θ = θmax*(2*rand()-1)
    s, c = sincos(θ)
    r = _r isa Nothing ? rand(1:3) : _r
    mat = if r == 1
        SMatrix{3,3}(1, 0, 0, 0, c, s, 0, -s, c)
    elseif r == 2
        SMatrix{3,3}(c, 0, -s, 0, 1, 0, s, 0, c)
    else
        SMatrix{3,3}(c, s, 0, -s, c, 0, 0, 0, 1)
    end
    refpos = positions[bead]
    [refpos + mat*(poss - refpos) for poss in positions]
end

choose_random_species(mc::MonteCarloSimulation) = rand(mc.indexof)
function compute_accept_move(before::MCEnergyReport, after::MCEnergyReport, T)
    b = Float64(before)
    a = Float64(after)
    a < b && return true
    e = exp((b-a)*u"K"/T)
    return rand() < e
end


function randomize_position!(positioni, j, bead, block, idx, atomblocks, d)
    pos = positioni[j]
    for _ in 1:30
        posr = random_rotation(random_rotation(random_rotation(pos, 90u"°", bead, 1), 90u"°", bead, 2), 90u"°", bead, 3)
        for _ in 1:1000
            post = random_translation(posr, d)
            if !inblockpocket(block, atomblocks, idx, post)
                positioni[j] = post
                return post
            end
        end
    end
    error(lazy"Could not find a suitable position for molecule $((i,j))!")
end

"""
    randomize_position!(mc::MonteCarloSimulation, idx, update_ewald=true)

Put the species at the given index to a random position and orientation.

!!! warn
    Using `update_ewald=false` will result in an inconsistent state for `mc` which will
    error on [`run_montecarlo`](@ref) unless [`baseline_energy(mc)`](@ref) is called.

!!! warn
    `update_ewald=true` requires a prior call to [`baseline_energy(mc)`](@ref).
"""
function randomize_position!(mc::MonteCarloSimulation, (i,j), update_ewald=true)
    d = norm(sum(mc.cell.mat; dims=2))/4
    post = randomize_position!(mc.positions[i], j, mc.bead[i], mc.blocks[i], d)
    if update_ewald
        single_contribution_ewald(mc.ewald, mc.offsets[i] + j, mc.speciesblocks[i], mc.idx[i], mc.atomblocks, d)
        update_mc!(mc, (i,j), post)
    else
        mc.ewald.last[] = -1 # signal the inconsistent state
    end
    post
end



"""
    run_montecarlo!(mc::MonteCarloSimulation, T, nsteps::Int)

Run a Monte-Carlo simulation at temperature `T` (given in K) during `nsteps`.
"""
function run_montecarlo!(mc::MonteCarloSimulation, T, nsteps::Int)
    energy = baseline_energy(mc)
    reports = [energy]
    running_update = @spawn nothing
    oldpos = SVector{3,typeof(1.0)}[]
    old_idx = (0,0)
    before = 0.0u"K"
    after = 0.0u"K"
    accepted = false
    for k in 1:nsteps
        idx = choose_random_species(mc)
        currentposition = (accepted&(old_idx==idx)) ? oldpos : mc.positions[idx[1]][idx[2]]
        idxi = mc.idx[idx[1]]
        newpos = if length(idxi) == 1
            random_translation(currentposition, 1.3u"Å")
        else
            if rand() < 0.5
                random_translation(currentposition, 1.3u"Å")
            else
                random_rotation(currentposition, 30.0u"°", mc.bead[idx[1]])
            end
        end
        inblockpocket(mc.speciesblocks[idx[1]], mc.atomblocks, idxi, newpos) && continue
        wait(running_update)
        if old_idx == idx
            if accepted
                before = after
            # else
            #     before = fetch(before_task) # simply keep its previous value
            end
            after = movement_energy(mc, idx, newpos)
        else
            before_task = @spawn movement_energy(mc, idx)
            after = movement_energy(mc, idx, newpos)
            before = fetch(before_task)
        end
        old_idx = idx
        accepted = compute_accept_move(before, after, T)
        oldpos = newpos
        if accepted
            running_update = @spawn update_mc!(mc, idx, oldpos) # do not use newpos since it can be changed in the next iteration before the Task is run
            diff = after - before
            if abs(Float64(diff.framework)) > 1e50 # an atom left a blocked pocket
                fetch(running_update)
                energy = baseline_energy(mc) # to avoid underflows
            else
                energy += diff
            end

            # wait(running_update)
            # shadow = MonteCarloSimulation(mc)
            # if !isapprox(Float64(baseline_energy(shadow)), Float64(energy), rtol=1e-4)
            #     println("discrepancy ! ", Float64(baseline_energy(shadow)), " vs ", Float64(energy))
            #     @show idx, k
            #     push!(reports, energy)
            #     return reports
            # end
            # println(Float64(energy), " vs ", Float64(baseline_energy(shadow)))
            # @show shadow.ewald.ctx.Eiks[1][30]
            # display(energy)
            # display(baseline_energy(shadow))
            # println(mc.positions)
        end
        push!(reports, energy)
    end
    wait(running_update)
    reports
end
