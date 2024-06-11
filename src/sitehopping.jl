export SiteHopping

"""
    SiteHopping

Structure representing a finite set of sites between which a mobile species can jump.
A `SiteHopping` instance stores both the sites and the energies relative to its possible
population, so it depends on both the sites and the framework.
"""
struct SiteHopping{Trng}
    population::Vector{Int}
    sites::Vector{Vector{SVector{3,TÅ}}}
    frameworkenergies::Vector{TK}
    interactions::Matrix{TK}
    rng::Trng
    tailcorrection::TK
end

function SiteHopping(sites, framework, forcefield_framework, syst_mol, population=-1;
                  supercell=nothing, rng=default_rng(), cutoff=12.0u"Å")
    print_charge_warning = PRINT_CHARGE_WARNING[]
    PRINT_CHARGE_WARNING[] = false
    mc1 = first(setup_montecarlo(framework, forcefield_framework, [syst_mol]; supercell, blockfiles=[false], cutoff, rng))
    baseline_energy(mc1)
    n = length(sites)
    frameworkenergies = Vector{TK}(undef, n)
    ewaldenergies = Vector{TK}(undef, n)
    sites = eltype(eltype(eltype(sites))) <: AbstractFloat ? [[mc1.step.mat * at for at in s] for s in sites] : sites
    for (i, s) in enumerate(sites)
        site = eltype(s) === AbstractFloat ? mc1.step.mat * s : s
        movement_energy(mc1, (1,1), site); update_mc!(mc1, (1,1), site)
        baseer = baseline_energy(mc1).er
        ewaldenergies[i] = baseer.inter + baseer.reciprocal
        frameworkenergies[i] = Number(baseer)
    end
    mc2 = first(setup_montecarlo(framework, forcefield_framework, [(syst_mol, 2)]; supercell, blockfiles=[false], cutoff, rng))
    baseline_energy(mc2)
    interactions = Matrix{TK}(undef, n, n)
    for (i, si) in enumerate(sites)
        interactions[i,i] = 1e100u"K"
        movement_energy(mc2, (1,1), si); update_mc!(mc2, (1,1), si)
        for j in (i+1:n)
            sj = sites[j]
            movement_energy(mc2, (1,2), sj); update_mc!(mc2, (1,2), sj)
            baseer2 = baseline_energy(mc2).er
            interactions[i,j] = interactions[j,i] = baseer2.inter + baseer2.reciprocal - ewaldenergies[i] - ewaldenergies[j]
        end
    end
    m = population isa Integer ? population : length(population)
    PRINT_CHARGE_WARNING[] = print_charge_warning
    mcpop = first(setup_montecarlo(framework, forcefield_framework, [(syst_mol, m)]; supercell, blockfiles=[false], cutoff, rng))
    tailcorrection = mcpop.tailcorrection[]
    actualpopulation = if population isa Integer
        m = length(mcpop.revflatidx)
        actualpopulation = collect(1:m)
    else
        copy(population)
    end
    SiteHopping(actualpopulation, sites, frameworkenergies, interactions, rng, tailcorrection)
end

function baseline_energy(sh::SiteHopping)
    inter = 0.0u"K"
    for (i, si) in enumerate(sh.population), j in (i+1:length(sh.population))
        sj = sh.population[j]
        inter += sh.interactions[si,sj]
    end
    sh.tailcorrection + sum(@view sh.frameworkenergies[sh.population]) + inter
end
