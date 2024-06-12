using Random: randperm

export SiteHopping

"""
    SiteHopping

Structure representing a finite set of sites between which a mobile species can jump.
A `SiteHopping` instance stores both the sites and the energies relative to its possible
population, so it depends on both the sites and the framework.
"""
struct SiteHopping{Trng}
    population::Vector{UInt16}
    sites::Vector{Vector{SVector{3,TÅ}}}
    frameworkenergies::Vector{TK}
    interactions::Matrix{TK}
    mat::SMatrix{3,3,TÅ,9}
    atomnames::Vector{Symbol}
    tailcorrection::TK
    rng::Trng
end

function SiteHopping(sites, framework, forcefield_framework, syst_mol, population=-1;
                  supercell=nothing, rng=default_rng(), cutoff=12.0u"Å")
    length(sites) > typemax(UInt16) && error("Cannot handle more than $(typemax(UInt16)) different sites")
    print_charge_warning = PRINT_CHARGE_WARNING[]
    PRINT_CHARGE_WARNING[] = false
    n = length(sites)
    frameworkenergies = Vector{TK}(undef, n)
    ewaldenergies = Vector{TK}(undef, n)
    mc1 = first(setup_montecarlo(framework, forcefield_framework, [syst_mol]; supercell, blockfiles=[false], cutoff, rng))
    mat = mc1.step.mat
    sites = eltype(eltype(eltype(sites))) <: AbstractFloat ? [[mat * at for at in s] for s in sites] : sites
    base1 = baseline_energy(mc1).er - movement_energy(mc1, (1,1))
    for (i, s) in enumerate(sites)
        Δ1 = base1 + movement_energy(mc1, (1,1), s)
        ewaldenergies[i] = Δ1.inter + Δ1.reciprocal
        frameworkenergies[i] = Number(Δ1)
    end
    mc2 = first(setup_montecarlo(framework, forcefield_framework, [(syst_mol, 2)]; supercell, blockfiles=[false], cutoff, rng))
    root = baseline_energy(mc2).er
    interactions = Matrix{TK}(undef, n, n)
    for (i, si) in enumerate(sites)
        interactions[i,i] = 1e100u"K"
        before, after = combined_movement_energy(mc2, (1,1), si)
        update_mc!(mc2, (1,1), si)
        root = root - before + after
        base2 = root - movement_energy(mc2, (1,2))
        for j in (i+1:n)
            sj = sites[j]
            Δ2 = base2 + movement_energy(mc2, (1,2), sj)
            interactions[i,j] = interactions[j,i] = Δ2.inter + Δ2.reciprocal - ewaldenergies[i] - ewaldenergies[j]
        end
    end
    m = population isa Integer ? population : length(population)
    PRINT_CHARGE_WARNING[] = print_charge_warning
    mcpop = first(setup_montecarlo(framework, forcefield_framework, [(syst_mol, m)]; supercell, blockfiles=[false], cutoff, rng))
    tailcorrection = mcpop.tailcorrection[]
    actualpopulation = if population isa Integer
        m = length(mcpop.revflatidx)
        actualpopulation = randperm(UInt16(n))
        resize!(actualpopulation, m)
    else
        copy(population)
    end
    M = maximum(actualpopulation)
    M > length(sites) && error(lazy"Cannot populate site $M among $(length(sites)) sites")
    atomnames = mcpop.step.ff.symbols[mcpop.step.ffidx[1]]
    SiteHopping(actualpopulation, sites, frameworkenergies, interactions, mat, atomnames,
                tailcorrection, rng)
end

function Base.show(io::IO, sh::SiteHopping)
    print(io, "SiteHopping setup with ", length(sh.population) , " molecules among ", length(sh.sites), " sites")
end

function baseline_energy(sh::SiteHopping)
    inter = 0.0u"K"
    for (i, si) in enumerate(sh.population), j in (i+1:length(sh.population))
        sj = sh.population[j]
        inter += sh.interactions[si,sj]
    end
    sh.tailcorrection + sum(@view sh.frameworkenergies[sh.population]) + inter
end

function movement_energy(sh::SiteHopping, i, si=sh.population[i])
    inter = 0.0u"K"
    for (j, sj) in enumerate(sh.population)
        i == j && continue
        inter += sh.interactions[si,sj]
    end
    sh.frameworkenergies[si] + inter
end
combined_movement_energy(sh::SiteHopping, i, si) = (movement_energy(sh, i), movement_energy(sh, i, si))

