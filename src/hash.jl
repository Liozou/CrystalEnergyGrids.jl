import StableHashTraits

StableHashTraits.hash_method(::CrystalEnergySetup) = StableHashTraits.StructHash()
StableHashTraits.hash_method(::EnergyGrid) = StableHashTraits.StructHash(filter(!=(:grid))∘fieldnames∘typeof => getfield)
StableHashTraits.hash_method(::GridCoordinatesSetup) = StableHashTraits.StructHash()
StableHashTraits.hash_method(::BlockFile) = StableHashTraits.StructHash(filter(!=(:block))∘fieldnames∘typeof => getfield)
StableHashTraits.hash_method(::EwaldFramework) = StableHashTraits.FnHash(x -> (x.kspace.ks, x.kspace.num_kvecs, x.precision))

function hash_string(x)
    h = StableHashTraits.stable_hash(x, StableHashTraits.HashVersion{2}())
    bytes2hex(h)
end
