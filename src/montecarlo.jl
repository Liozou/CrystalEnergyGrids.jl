export setup_montecarlo, move_cost

struct MonteCarloSimulation{TFramework}
    ff::ForceField
    ewald::IncrementalEwaldContext
    framework::TFramework
    coulomb::EnergyGrid
    grids::Vector{EnergyGrid} # grids[i] is the VdW grid for atom i in ff
    offsets::Vector{Int}
    # offsets[i] is the number of molecules belonging to a kind strictly lower than i
    idx::Vector{Vector{Int}}
    # idx[i][k] is ff.sdict[atomic_symbol(systemkinds[i], k)], i.e. the numeric index
    # in the force field for the k-th atom in a system of kind i
    charges::Vector{typeof(1.0u"e_au")}
    # charges[ix] is the charge of the atom whose of index ix in ff, i.e. the charge of the
    # k-th atom in a system of kind i is charges[idx[i][k]].
    positions::Vector{Vector{Vector{SVector{3,typeof(1.0u"Å")}}}}
    # positions[i][j][k] is the position of the k-th atom in the j-th system of kind i
end

function SimulationStep(mc::MonteCarloSimulation)
    SimulationStep(mc.ff, mc.charges, mc.positions, trues(length(mc.idx)), mc.idx)
end

function setup_montecarlo(framework::AbstractSystem{3}, systems::Vector{T}, ff::ForceField,
                          eframework::EwaldFramework, coulomb::EnergyGrid, grids::Vector{EnergyGrid}) where {T<:AbstractSystem{3}}
    kindsdict = Dict{Vector{Symbol},Int}()
    systemkinds = T[]
    U = Vector{SVector{3,typeof(1.0u"Å")}} # positions of the atoms of a system
    poss = Vector{U}[]
    indices = Tuple{Int,Int}[]
    for system in systems
        n = length(kindsdict)+1
        kind = get!(kindsdict, atomic_symbol(system), n)
        if kind === n
            push!(systemkinds, system)
            push!(poss, U[])
        end
        push!(poss[kind], position(system))
        push!(indices, (kind, length(poss[kind])))
    end

    idx = [[ff.sdict[atomic_symbol(s, k)] for k in 1:length(s)] for s in systemkinds]
    charges = fill(NaN*u"e_au", length(ff.sdict))
    for (i, ids) in enumerate(idx)
        kindi = systemkinds[i]
        for (k, ix) in enumerate(ids)
            charges[ix] = kindi[k,:atomic_charge]
        end
    end


    n = length(poss)
    offsets = Vector{Int}(undef, n)
    offsets[1] = 0
    for i in 1:(n-1)
        offsets[i+1] = offsets[i] + length(poss[i])
    end

    ewald = IncrementalEwaldContext(EwaldContext(eframework, systems))

    MonteCarloSimulation(ff, ewald, framework, coulomb, grids, offsets, idx, charges, poss), indices
end

function setup_montecarlo(framework, forcefield_framework::String, systems::Vector{<:AbstractSystem{3}};
                          gridstep=0.15, supercell=nothing, new=false)
    syst_framework = load_framework_RASPA(framework, forcefield_framework)
    ff = parse_forcefield_RASPA(forcefield_framework)
    mat = stack3(bounding_box(syst_framework))
    if supercell isa Nothing
        supercell = find_supercell(syst_framework, 12.0u"Å")
    end
    supercell::NTuple{3,Int}

    needcoulomb = false
    encountered_atoms = Set{Symbol}()
    for system in systems
        if !needcoulomb
            needcoulomb = any(!iszero(system[i,:atomic_charge])::Bool for i in 1:length(system))
        end
        for k in 1:length(system)
            push!(encountered_atoms, atomic_symbol(system, k))
        end
    end
    atoms = [(atom, ff.sdict[atom]) for atom in encountered_atoms]
    sort!(atoms; by=last)

    coulomb_grid_path, vdw_grid_paths = grid_locations(framework, forcefield_framework, first.(atoms), gridstep, supercell)

    coulomb, eframework = if needcoulomb
        _eframework = initialize_ewald(syst_framework)
        retrieve_or_create_grid(coulomb_grid_path, syst_framework, ff, gridstep, _eframework, mat, new), _eframework
    else
        EnergyGrid(), EwaldFramework(mat)
    end

    grids = Vector{EnergyGrid}(undef, length(ff.sdict))
    for (j, (atom, i)) in enumerate(atoms)
        grids[i] = retrieve_or_create_grid(vdw_grid_paths[j], syst_framework, ff, gridstep, atom, mat, new)
    end

    setup_montecarlo(syst_framework, systems, ff, eframework, coulomb, grids)
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
function Base.show(io::IO, ::MIME"text/plain", e::MCEnergyReport)
    println(io, Float64(e), " = ", e.framework.vdw, " + ", e.framework.direct, " + ", e.inter, " + ", e.reciprocal)
end

for op in (:+, :-)
    @eval begin
        function Base.$(op)(f1::FrameworkEnergyReport, f2::FrameworkEnergyReport)
            FrameworkEnergyReport($op(f1.vdw, f2.vdw), $op(f1.direct, f2.direct))
        end
        function Base.$(op)(e1::MCEnergyReport, e2::MCEnergyReport)
            MCEnergyReport($op(e1.framework, e2.framework), $op(e1.inter, e2.inter), $op(e1.reciprocal, e2.reciprocal))
        end
    end
end

"""
    framework_interactions(mc::MonteCarloSimulation, i, positions)

Energy contribution of the interaction between the framework and a molecule of system kind
`i` at the given `positions`.

Only return the Van der Waals and the direct part of the Ewald summation. The reciprocal
part can be obtained with [`compute_ewald`](@ref).
"""
function framework_interactions(mc::MonteCarloSimulation, indices::Vector{Int}, positions)
    isempty(mc.framework) && return FrameworkEnergyReport()
    n = length(indices)
    vdw = 0.0u"K"
    direct = 0.0u"K"
    hascoulomb = mc.coulomb.ε_Ewald != -Inf
    for k in 1:n
        ix = indices[k]
        pos = positions[k]
        vdw += interpolate_grid(mc.grids[ix], pos)
        hascoulomb && (direct += Float64(mc.charges[ix]/u"e_au")*interpolate_grid(mc.coulomb, pos))
    end
    return FrameworkEnergyReport(vdw, direct)
end
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
    isempty(mc.framework) && return MCEnergyReport(fer, vdw, reciprocal)
    for (i, indices) in enumerate(mc.idx)
        poss_i = mc.positions[i]
        for poss in poss_i
            fer += framework_interactions(mc, indices, poss)
        end
    end
    return MCEnergyReport(fer, vdw, reciprocal)
end

"""
    movement_energy(mc::MonteCarloSimulation, (i, j), positions=mc.positions[i][j])

Compute the energy contribution of the `j`-th molecule of kind `i` when placed at
`positions`. If not provided, `positions` is the current position for that molecule.

The energy difference between the new position for the molecule and the current one is
`movement_energy(mc, (i,j), positions) - movement_energy(mc, (i,j))`.
"""
function movement_energy(mc::MonteCarloSimulation, idx, positions=mc.positions[idx[1]][idx[2]])
    i, j = idx
    k = mc.offsets[i]+j
    singlereciprocal = single_contribution_ewald(mc.ewald, k, positions)
    singlevdw = single_contribution_vdw(SimulationStep(mc), (i,j), positions)
    fer = framework_interactions(mc, i, positions)
    MCEnergyReport(fer, singlevdw, singlereciprocal)
end
