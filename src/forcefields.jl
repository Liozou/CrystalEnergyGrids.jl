# Definition of a force field

export ForceField

struct ForceField
    interactions::Matrix{Union{InteractionRule,InteractionRuleSum}}
    sdict::IdDict{Symbol,Int}
    symbols::Vector{Symbol}
    cutoff::TÅ
    name::String
end

struct DoubleDefinedInteractionRule <: Exception
    a::Symbol
    b::Symbol
    A::InteractionRule
    B::InteractionRule
end
function Base.showerror(io::IO, e::DoubleDefinedInteractionRule)
    print(io, "Interaction between species ", e.a, " and ", e.b, " defined at multiple occurences: ", e.A, " vs ", e.B)
end

struct AsymetricSelfInteractionRule <: Exception
    s::Symbol
    rule::InteractionRule
end
function Base.showerror(io::IO, e::AsymetricSelfInteractionRule)
    ofs = e.rule.kind === FF.CoulombEwaldDirect
    print(io, e.rule.kind, " interaction between two species ", e.s, " should be defined by a single parameter specific to ",
              e.s, " but two were given: ", e.rule.params[1+ofs], " ≠ ", e.rule.params[2+ofs])
end

struct IncompatibleMixing <: Exception
    A::InteractionRule
    B::InteractionRule
    mixing::FF.MixingRule
end
function Base.showerror(io::IO, e::IncompatibleMixing)
    print(io, "Cannot mix interactions of kind ", e.A.kind, " and ", e.B.kind, " using ", e.mixing, " mixing rules")
end


function _mix_rules(A::InteractionRule, B::InteractionRule, mixing::FF.MixingRule)
    A, B = minmax(A, B) # commented conditions are implied by this ordering
    if B.kind === FF.NoInteraction || B.kind === FF.UndefinedInteraction # || A.kind === B.kind
        InteractionRule(B.kind, Float64[])
    elseif B.kind === FF.HardSphere # && A.kind === FF.HardSphere
        FF.HardSphere(A.params[1], B.params[1])
    elseif A.kind === FF.HardSphere # && B.kind !== FF.HardSphere
        FF.HardSphere(A.params[1], 0.0)
    elseif B.kind === FF.CoulombEwaldDirect || B.kind === FF.Coulomb
        A.kind === B.kind || error("Cannot mix cutoff and no-cutoff coulomb interactions")
        if B.kind === FF.CoulombEwaldDirect
            A.params[1] == B.params[1]  || error("Cannot use two Ewald summations with different cutoffs")
            FF.CoulombEwaldDirect(A.params[1], A.params[2], B.params[2])
        else
            FF.Coulomb(A.params[1], B.params[1])
        end
    elseif A.kind === FF.CoulombEwaldDirect || A.kind === FF.Coulomb
        # && B.kind !== FF.CoulombEwaldDirect && B.kind != FF.Coulomb
        FF.NoInteraction()
    elseif B.kind === FF.LennardJones # && A.kind === FF.LennardJones
        εA, σA = A.params
        εB, σB = B.params
        if mixing === FF.LorentzBerthelot
            FF.LennardJones(sqrt(εA*εB), (σA+σB)/2)
        elseif mixing === FF.WaldmanHagler
            σA3 = σA^3
            σB3 = σB^3
            σAB6 = (σA3*σA3 + σB3*σB3)/2
            FF.LennardJones(sqrt(εA*εB)*σA3*σB3/σAB6, σAB6^(1/6))
        elseif mixing === FF.Geometric
            FF.LennardJones(sqrt(εA*εB), sqrt(σA*σB))
        else
            @assert false
        end
    elseif A.kind !== B.kind
        FF.UndefinedInteraction()
    elseif B.kind === FF.Buckingham
        FF.Buckingham(sqrt(A.params[1]*B.params[1]), sqrt(A.params[2]*B.params[2]), sqrt(A.params[3]*B.params[3]))
    elseif B.kind === FF.Monomial
        if A.params[2] != B.params[2]
            FF.UndefinedInteraction()
        else
            FF.Geometric(sqrt(A.params[1]*B.params[1]), A.params[2])
        end
    elseif B.kind === FF.Exponential
        FF.Expontential(sqrt(A.params[1]*B.params[1]), sqrt(A.params[2]*B.params[2]))
    end
end

function reduce_rule_sum!(rules::Vector{InteractionRule})
    sort!(rules)
    i = length(rules)
    while i > 0 && rules[i].kind === FF.NoInteraction
        i -= 1
    end
    i == 0 && return FF.NoInteraction()
    i == 1 && return rules[1]
    rules[i].kind === FF.UndefinedInteraction && return FF.UndefinedInteraction()
    resize!(rules, i)
    k = 1
    while k <= i && rules[k].kind === FF.HardSphere
        k += 1
    end
    if k > 1
        rule0 = popfirst!(rules)
        r1, r2 = rule0.params
        for _ in 2:(k-1)
            r = popfirst!(rules)
            r1 = max(r1, r.params[1])
            r2 = max(r2, r.params[2])
        end
        @assert r1 == rule0.params[1]
        newhardsphere = FF.HardSphere(r1, r2)
        k > i && return newhardsphere
        pushfirst!(rules, newhardsphere)
    end
    return InteractionRuleSum(rules)
end

function asymetric_mix_rules(rule::InteractionRule, is::InteractionRuleSum, mixing::FF.MixingRule)
    (rule.kind === FF.NoInteraction) || (rule.kind === FF.UndefinedInteraction) && return rule
    interactions = [_mix_rules(rule, r, mixing) for r in is.rules]
    return reduce_rule_sum!(interactions)
end

function mix_rules(rulei, rulej, mixing::FF.MixingRule)
    if rulei isa InteractionRule
        rulej isa InteractionRule && return _mix_rules(rulei, rulej, mixing)
        return asymetric_mix_rules(rulei, rulej, mixing)
    end
    rulej isa InteractionRule && return asymetric_mix_rules(rulej, rulei, mixing)
    i = 1
    j = 1
    n = length(rulei.rules)
    m = length(rulej.rules)
    rules = InteractionRule[]
    while i ≤ n && m ≤ j
        A = rulei.rules[i]
        B = rulei.rules[j]
        if A.kind === B.kind
            push!(rules, _mix_rules(A, B, mixing))
            i += 1
            j += 1
        elseif A.kind === FF.HardSphere
            push!(rules, FF.HardSphere(A.params[1], 0.0))
            i += 1
        elseif B.kind === FF.HardSphere
            push!(rules, FF.HardSphere(B.params[1], 0.0))
            j += 1
        elseif A.kind === FF.Coulomb || A.kind === FF.CoulombEwaldDirect
            i += 1
        elseif B.kind === FF.Coulomb || B.kind === FF.CoulombEwaldDirect
            j += 1
        else
            return FF.UndefinedInteraction()
        end
    end
    reduce_rule_sum!(rules)
end

function map_rule(f, rule)
    if rule isa InteractionRule
        f(rule)
    else
        InteractionRuleSum([f(r) for r in rule.rules])
    end
end


"""
    ForceField(input, mixing::FF.MixingRule=FF.ErrorOnMix, cutoff=12.0u"Å", shift::Bool=true, tailcorrection::Bool=!shift, sdict::Union{Nothing,IdDict{Symbol,Int}}=nothing)

Build a force field from a description of the pairwise interactions.

`input` must be an iterable over pairs of the form `(:A, :B) => rule` where
- `:A` and `:B` are the `Symbol` of two species as returned by `AtomsBase.atomic_symbol`
- `rule` is their interaction potential, which can be either an [`InteractionRule`](@ref)
or an [`InteractionRuleSum`](@ref)

Pair interactions are only considered between particles whose distance is below the `cutoff`.
The `cutoff` can be set to `Inf` to remove it.

Missing pairs of species are attributed an interaction according to the given `mixing` rule,
see [`FF.MixingRule`](@ref).
If no appropriate mixing rule exists for the specified interaction, an [`FF.UndefinedInteraction`](@ref))
is set and will error if ever called.

If `tailcorrection` is set, a tail-correction contribution to the energy is computed.
If `shift` is set, the pair potentials are shifted to equal zero at the cutoff.

If provided, `sdict` should be an `IdDict{Symbol,Int}` linking each species symbol to its
unique identifier. Identifiers populate the range `1:n` (where `n == length(sdict)`).
"""
function ForceField(input, mixing::FF.MixingRule=FF.ErrorOnMix, cutoff=12.0u"Å", shift::Bool=true, tailcorrection::Bool=!shift;
                    sdict::Union{Nothing,IdDict{Symbol,Int}}=nothing, name::String="(unnamed)")
    cut = @convertifnotfloat u"Å" cutoff
    smap::IdDict{Symbol,Int} = if sdict isa IdDict{Symbol,Int}
        sdict
    else
        allatoms = Vector{Symbol}(undef, 2*length(input))
        for (i, ((a, b), _)) in enumerate(input)
            allatoms[2*i-1] = a
            allatoms[2*i] = b
        end
        sort!(allatoms)
        unique!(allatoms)
        IdDict([s => i for (i,s) in enumerate(allatoms)])
    end
    n = length(smap)
    symbols = Vector{Symbol}(undef, n)
    for (s, i) in smap
        symbols[i] = Symbol(identify_atom(s))
    end
    interactions = Matrix{Union{InteractionRule,InteractionRuleSum}}(undef, n, n)
    done = falses(n, n)
    for ((a, b), rule) in input
        i = smap[a]
        j = smap[b]
        if i == j
            map_rule(rule) do r
                ofs = r.kind === FF.CoulombEwaldDirect
                if (ofs || r.kind === FF.HardSphere || r.kind === FF.Coulomb) && r.params[1+ofs] != r.params[2+ofs]
                    throw(AsymetricSelfInteractionRule(a, r))
                end
                r
            end
        end
        if done[i,j] && interactions[i,j] != rule
            throw(DoubleDefinedInteractionRule(a, b, interactions[i,j], rule))
        end
        interactions[i,j] = interactions[j,i] = rule
        done[i,j] = done[j,i] = true
    end
    for i in 1:n, j in (i+1):n
        done[i,j] && continue
        if mixing === FF.IgnoreInteraction
            interactions[i,j] = interactions[j,i] = FF.NoInteraction()
            continue
        elseif mixing === FF.ErrorOnMix || !done[i,i] || !done[j,j]
            interactions[i,j] = interactions[j,i] = FF.UndefinedInteraction()
            continue
        end
        rulei = interactions[i,i]
        rulej = interactions[j,j]
        interactions[i,j] = interactions[j,i] = mix_rules(rulei, rulej, mixing)
    end

    for i in 1:n, j in (i+1):n
        interactions[i,j] = interactions[j,i] = map_rule(interactions[i,j]) do r
            InteractionRule(r.kind, r.params, shift, cut, tailcorrection)
        end
    end

    ForceField(interactions, smap, symbols, cut*u"Å", name)
end

function Base.show(io::IO, ::MIME"text/plain", ff::ForceField)
    print(io, "Force field ", ff.name)
    if isinf(ff.cutoff)
        print(io, " with no cutoff")
    else
        print(io, " with cutoff of ", ff.cutoff)
    end
    println(io, " (ᵀ = tail correction, ₛ = shifted)")
    revdict = Dict([i => at for (at, i) in ff.sdict])
    for i in 1:length(revdict)
        ati = revdict[i]
        for j in i:length(revdict)
            atj = revdict[j]
            print(io, " * ", ati, " - ", atj, " => ")
            print(IOContext(io, :inforcefield=>true), ff.interactions[i,j])
            println(io)
        end
    end
end

Base.getindex(ff::ForceField, i::Integer, j::Integer) = ff.interactions[i,j]
Base.getindex(ff::ForceField, a::Symbol, b::Symbol) = ff[ff.sdict[a], ff.sdict[b]]

# function (ff::ForceField)(i::Integer, j::Integer, distance)
#     if distance isa (Quantity{T,Unitful.𝐋^2,U} where {T,U})
#         distance >= ff.cutoff*ff.cutoff && return 0.0u"K"
#     elseif distance isa (Quantity{T,Unitful.𝐋,U} where {T,U})
#         distance >= ff.cutoff && return 0.0u"K"
#     else
#         NoUnits(distance/u"Å") > ff.cutoff && return 0.0u"K"
#     end
#     return ff[i,j](distance)
# end
# (ff::ForceField)(a::Symbol, b::Symbol, distance) = ff(ff.sdict[a], ff.sdict[b], distance)

function derivatives_nocutoff(ff::ForceField, i::Integer, j::Integer, distance)
    derivativesGrid(ff.interactions[i,j], distance)
end

function needsvdwgrid(ff::ForceField, atom)
    i = ff.sdict[Symbol(get_atom_name(atom))]
    for inters in @view ff.interactions[:,i]
        for inter in (inters isa InteractionRule ? (inters,) : inters.rules)
            if inter.kind != FF.NoInteraction && inter.kind != FF.CoulombEwaldDirect
                return true
            end
        end
    end
    return false
end
