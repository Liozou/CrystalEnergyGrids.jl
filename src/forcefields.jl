struct ForceField
    interactions::Matrix{Union{InteractionRule,InteractionRuleSum}}
end

"""
    ForceField(input, mixing::FF.MixingRule=FF.ErrorOnMix)

Build a force field from a description of the pairwise interactions.

`input` must be a list of triplets of the form `(:A, :B, rule)` where
- `:A` and `:B` are the `Symbol` of two species as returned by `AtomsBase.atomic_symbol`
- `rule` is their interaction potential, which can be either an [`InteractionRule`](@ref)
  or an [`InteractionRuleSum`](@ref)

Missing pairs of species will be attributed an interaction according to the given
`mixing` rule.
"""
function ForceField(input, mixing::FF.MixingRule=FF.ErrorOnMix)
end
