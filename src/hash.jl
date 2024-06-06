import StableHashTraits

StableHashTraits.hash_method(::CrystalEnergySetup) = StableHashTraits.StructHash()
StableHashTraits.hash_method(::EnergyGrid) = StableHashTraits.StructHash(filter(!=(:grid))∘fieldnames∘typeof => getfield)
StableHashTraits.hash_method(::GridCoordinatesSetup) = StableHashTraits.StructHash()
StableHashTraits.hash_method(::BlockFile) = StableHashTraits.StructHash(filter(!=(:block))∘fieldnames∘typeof => getfield)
StableHashTraits.hash_method(::EwaldFramework) = StableHashTraits.FnHash(x -> (x.kspace.ks, x.kspace.num_kvecs, x.precision))

StableHashTraits.hash_method(::EwaldContext) = StableHashTraits.FnHash(x -> (x.eframework, x.charges, x.atoms, x.speciesof, ustrip(u"K", x.static_contribution[])))
StableHashTraits.hash_method(::SimulationStep) = StableHashTraits.FnHash(x -> (ustrip.(u"e_au", x.charges), x.atoms, [[[ustrip(u"Å", k) for k in j] for j in i] for i in x.positions], ustrip.(u"Å", x.mat)))
StableHashTraits.hash_method(::MonteCarloSetup) = StableHashTraits.FnHash(x -> (x.step, x.ewald.ctx, x.coulomb, x.speciesblocks, x.bead))

StableHashTraits.hash_method(::InteractionRule) = StableHashTraits.StructHash()
StableHashTraits.hash_method(::InteractionRuleSum) = StableHashTraits.StructHash()
StableHashTraits.hash_method(::ForceField) = StableHashTraits.StructHash()

StableHashTraits.hash_method(::PseudoAtomInfo) = StableHashTraits.StructHash()
StableHashTraits.hash_method(::PseudoAtomListing) = StableHashTraits.FnHash(pal -> begin
    kvs = collect(pal.exact)
    I = sortperm(kvs; by=last)
    (first.(kvs)[I], pal.info[I])
end)

function hash_string(x)
    h = StableHashTraits.stable_hash(x, StableHashTraits.HashVersion{2}())
    bytes2hex(h)
end
