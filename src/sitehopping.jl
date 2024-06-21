using Random: randperm, shuffle!

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
    m = population isa Integer ? population : length(population)
    mcpop = first(setup_montecarlo(framework, forcefield_framework, [(syst_mol, m)]; supercell, blockfiles=[false], cutoff, rng))
    mat = mcpop.step.mat
    sites = eltype(eltype(eltype(sites))) <: AbstractFloat ? [[mat * at for at in s] for s in sites] : sites
    n = length(sites)
    atomnames = mcpop.step.ff.symbols[mcpop.step.ffidx[1]]

    storeinteractions = joinpath(scratchspace, string(join(atomnames), '-', hash_string(string(hash_string(mcpop.step.ff), '-', hash_string(sites)))))
    storeframework = string(storeinteractions, '-', hash_string(string(hash_string(mcpop.coulomb.grid), '-', isempty(mcpop.grids) ? "0" : hash_string([mcpop.grids[k].grid for k in mcpop.step.ffidx[1]]))))
    frameworkenergies, interactions = if isfile(storeframework) && isfile(storeinteractions)
        printstyled("Retrieved framework site energies at ", storeframework, " and site interactions at ", storeinteractions, ".\n"; color=:cyan)
        _frameworkenergies = deserialize(storeframework)::Vector{TK}
        _interactions = deserialize(storeinteractions)::Matrix{TK}
        (_frameworkenergies, _interactions)
    else
        frameworkenergiesC = Vector{TK}(undef, n)
        ewaldenergies = Vector{TK}(undef, n)
        print_charge_warning = PRINT_CHARGE_WARNING[]
        PRINT_CHARGE_WARNING[] = false
        mc1 = first(setup_montecarlo(framework, forcefield_framework, [syst_mol]; supercell, blockfiles=[false], cutoff, rng))
        PRINT_CHARGE_WARNING[] = print_charge_warning
        printstyled("Creating framework site energies at ", storeframework, " ... "; color=:yellow)
        base1 = baseline_energy(mc1).er - movement_energy(mc1, (1,1))
        for (i, s) in enumerate(sites)
            Δ1 = base1 + movement_energy(mc1, (1,1), s)
            ewaldenergies[i] = Δ1.inter + Δ1.reciprocal
            frameworkenergiesC[i] = Number(Δ1)
        end
        serialize(storeframework, frameworkenergiesC)
        printstyled("Done.\n"; color=:cyan)
        if isfile(storeinteractions)
            _interactions2 = deserialize(storeinteractions)::Matrix{TK}
            printstyled("Retrieved site interactions at ", storeinteractions, ".\n"; color=:cyan)
            frameworkenergiesC, _interactions2
        else
            PRINT_CHARGE_WARNING[] = false
            mc2 = first(setup_montecarlo(framework, forcefield_framework, [(syst_mol, 2)]; supercell, blockfiles=[false], cutoff, rng))
            PRINT_CHARGE_WARNING[] = print_charge_warning
            printstyled("Creating site interactions at ", storeinteractions, " ... "; color=:yellow)
            root = baseline_energy(mc2).er
            interactionsC = Matrix{TK}(undef, n, n)
            for (i, si) in enumerate(sites)
                interactionsC[i,i] = 1e100u"K"
                before, after = combined_movement_energy(mc2, (1,1), si)
                update_mc!(mc2, (1,1), si)
                root = root - before + after
                base2 = root - movement_energy(mc2, (1,2))
                for j in (i+1:n)
                    sj = sites[j]
                    Δ2 = base2 + movement_energy(mc2, (1,2), sj)
                    interactionsC[i,j] = interactionsC[j,i] = Δ2.inter + Δ2.reciprocal - ewaldenergies[i] - ewaldenergies[j]
                end
            end
            serialize(storeinteractions, interactionsC)
            printstyled("Done.\n"; color=:cyan)
            frameworkenergiesC, interactionsC
        end
    end
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
    SiteHopping(actualpopulation, sites, frameworkenergies, interactions, mat, atomnames,
                tailcorrection, rng)
end

function Base.show(io::IO, sh::SiteHopping)
    print(io, "SiteHopping setup with ", length(sh.population) , " molecules among ", length(sh.sites), " sites")
end
function Base.copy(sh::SiteHopping)
    SiteHopping(copy(sh.population), sh.sites, sh.frameworkenergies, sh.interactions, sh.mat, sh.atomnames, copy(sh.tailcorrection), copy(sh.rng))
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

"""
    extremal_placements(sh::SiteHopping, sitemap::AbstractVector{T}) where T

Return a list of pairs `priorities => energy` that give the energy corresponding to placing
the mobile species of `sh` according to `priorities`. The list is sorted by `energy`.

`sitemap` is the list such that `sitemap[i]` is the name of site `i`. Given this, a
`priorities` equal to `[p1, p2, ..., pN]` means: fill all the sites whose name is `p1`,
then fill all the sites whose name is `p2`, etc., until trying to fill the sites whose name
is `pN` and placing the last species. If the sites `pN` are not filled, the returned energy
is the minimum over 100 random assignments among those sites.
"""
function extremal_placements(sh::SiteHopping, sitemap::AbstractVector{T}) where {T}
    allsites = unique!(sort(sitemap))
    N = length(allsites)
    allsitesdict = Dict(n => i for (i,n) in enumerate(allsites))
    newsitemap = [allsitesdict[name] for name in sitemap]
    sites_per_name = [UInt16[] for _ in 1:N]
    for (i, site) in enumerate(newsitemap)
        push!(sites_per_name[site], i)
    end
    ret = Pair{Vector{Int},Tuple{TK,Vector{UInt16}}}[]
    newsh = copy(sh)
    m = length(sh.population)
    empty!(newsh.population)
    _extremal_placements!(ret, newsh, Int[], BitSet(), m, N, newsitemap, sites_per_name)
    sort!([[allsites[i] for i in t] => (e, sort!(l)) for (t,(e,l)) in ret]; by=first∘last)
end

function _extremal_placements!(ret, sh::SiteHopping, priorities, used, m, N, newsitemap, sites_per_name)
    for site in 1:N
        site in used && continue
        idxs = sites_per_name[site]
        n = length(idxs)
        currpopsize = length(sh.population)
        currpopsize + n == m || isempty(priorities) || site > priorities[end] || continue
        push!(priorities, site)
        r = m - currpopsize
        if r >= n
            append!(sh.population, idxs)
            if currpopsize + n == m
                push!(ret, copy(priorities) => (baseline_energy(sh), copy(sh.population)))
            else
                _extremal_placements!(ret, copy(sh), copy(priorities), union(used, site), m, N, newsitemap, sites_per_name)
            end
        else
            perm = randperm(n)
            append!(sh.population, @views idxs[perm[1:r]])
            min_e = baseline_energy(sh)
            min_pop = sh.population[end-r+1:end]
            RETRIES = 100
            for _ in 2:RETRIES
                resize!(sh.population, currpopsize)
                shuffle!(perm)
                append!(sh.population, @views idxs[perm[1:r]])
                @assert length(sh.population) == m
                new_e = baseline_energy(sh)
                if new_e < min_e
                    min_pop = sh.population[end-r+1:end]
                    min_e = new_e
                end
            end
            push!(ret, copy(priorities) => (min_e, [(@view sh.population[1:end-r]); min_pop]))
        end

        pop!(priorities)
        resize!(sh.population, currpopsize)
    end
end
