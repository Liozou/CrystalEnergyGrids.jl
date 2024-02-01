import StableHashTraits

StableHashTraits.hash_method(::CrystalEnergySetup) = StableHashTraits.StructHash()
StableHashTraits.hash_method(::EnergyGrid) = StableHashTraits.StructHash(filter(!=(:grid))∘fieldnames∘typeof => getfield)
StableHashTraits.hash_method(::GridCoordinatesSetup) = StableHashTraits.StructHash()
StableHashTraits.hash_method(::BlockFile) = StableHashTraits.StructHash(filter(!=(:block))∘fieldnames∘typeof => getfield)
StableHashTraits.hash_method(::EwaldFramework) = StableHashTraits.FnHash(x -> (x.kspace.ks, x.kspace.num_kvecs, x.precision))

StableHashTraits.hash_method(::EwaldContext) = StableHashTraits.FnHash(x -> (x.eframework, x.charges, x.atoms, ustrip(u"K", x.static_contribution)))
StableHashTraits.hash_method(::SimulationStep) = StableHashTraits.FnHash(x -> (ustrip.(u"e_au", x.charges), x.atoms, [[[ustrip(u"Å", k) for k in j] for j in i] for i in x.positions], ustrip.(u"Å", x.mat)))
StableHashTraits.hash_method(::MonteCarloSetup) = StableHashTraits.FnHash(x -> (x.step, x.ewald.ctx, x.coulomb, x.speciesblocks, x.bead))

function hash_string(x)
    h = StableHashTraits.stable_hash(x, StableHashTraits.HashVersion{2}())
    bytes2hex(h)
end
