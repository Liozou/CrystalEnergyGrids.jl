import StableHashTraits: transformer
using StableHashTraits: Transformer, omit_fields, pick_fields, stable_hash, HashVersion

# transformer(::Type{<:CrystalEnergySetup}) = # DEFAULT
transformer(::Type{EnergyGrid}) = Transformer(omit_fields(:grid))
# transformer(::Type{GridCoordinatesSetup}) = # DEFAULT
transformer(::Type{BlockFile}) = Transformer(omit_fields(:block))
transformer(::Type{EwaldFramework}) = Transformer(x -> (x.kspace.ks, x.kspace.num_kvecs, x.precsision); hoist_type=true)

function transformer(::Type{EwaldContext})
    Transformer(x -> (x.eframework, x.charges, x.atoms, x.speciesof, ustrip(u"K", x.static_contribution[])); hoist_type=true)
end
function transformer(::Type{SimulationStep{T}}) where T
    Transformer(x -> (ustrip.(u"e_au", x.charges), x.atoms, [[[ustrip(u"Å", k) for k in j] for j in i] for i in x.positions], ustrip.(u"Å", x.mat)); hoist_type=true)
end
function transformer(::Type{MonteCarloSetup{T,Trng}}) where {T,Trng}
    Transformer(x -> (x.step, x.ewald.ctx, x.coulomb, x.speciesblocks, x.bead); hoist_type=true)
end

# transformer(::Type{InteractionRule}) = # DEFAULT
# transformer(::Type{InteractionRuleSum}) = # DEFAULT
# transformer(::Type{ForceField}) = # DEFAULT

# transformer(::Type{PseudoAtomInfo}) = # DEFAULT
transformer(::Type{PseudoAtomListing}) = Transformer(pal -> begin
    kvs = collect(pal.exact)
    I = sortperm(kvs; by=last)
    (first.(kvs)[I], pal.info[I])
end; hoist_type=true)

function hash_string(x)
    h = stable_hash(x, HashVersion{4}())
    bytes2hex(h)
end
