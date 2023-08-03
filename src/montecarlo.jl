export setup_montecarlo, move_cost

struct MonteCarloSimulation{TFramework,TMolecules}
    ff::ForceField
    ewald::EwaldContext
    framework::TFramework
    coulomb::EnergyGrid
    grids::Vector{EnergyGrid} # grids[i] is the VdW grid for atom i in ff
    systemkinds::Vector{TMolecules}
    idx::Vector{Vector{Int}}
    # idx[i][k] is ff.sdict[atomic_symbol(systemkinds[i], k)], i.e. the numeric index
    # in the force field for the k-th atom in a system of kind i
    positions::Vector{Vector{Vector{SVector{3,typeof(1.0u"Å")}}}}
    # positions[i][j][k] is the position of the k-th atom in the j-th system of kind i
end

function setup_montecarlo(framework::AbstractSystem{3}, systems::Vector{T}, ff::ForceField,
                          eframework::EwaldFramework, coulomb::EnergyGrid, grids::Vector{EnergyGrid}) where {T<:AbstractSystem{3}}
    kindsdict = Dict{Vector{Symbol},Int}()
    systemkinds = T[]
    U = Vector{SVector{3,typeof(1.0u"Å")}} # positions of the atoms of a system
    poss = Vector{U}[]
    indices = Tuple{Int,Int}[]
    needcoulomb = false
    for system in systems
        n = length(kindsdict)+1
        kind = get!(kindsdict, atomic_symbol(system), n)
        if kind === n
            push!(systemkinds, system)
            push!(poss, U[])
        end
        push!(poss[kind], position(system))
        push!(indices, (kind, length(poss[kind])))
        if !needcoulomb
            needcoulomb = any(!iszero(system[i,:atomic_charge])::Bool for i in 1:length(system))
        end
    end

    ewald = EwaldContext(eframework, systems)
    idx = [[ff.sdict[atomic_symbol(s, k)] for k in 1:length(s)] for s in systemkinds]

    MonteCarloSimulation(ff, ewald, framework, coulomb, grids, systemkinds, idx, poss)
end

function setup_montecarlo(framework::String, forcefield_framework::String, systems::Vector{<:AbstractSystem{3}};
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

function move_cost(mc::MonteCarloSimulation, idx, newpositions)
    
end
