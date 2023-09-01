## Definition of pair interactions and mixing rules
export FF, InteractionRule, InteractionRuleSum

using SpecialFunctions: erfc

module FF
    """
        InteractionKind

    Nature of interaction that can be computed between two species.

    See also [`InteractionRule`](@ref).

    Current list:
    [`LennardJones`](@ref),
    [`Coulomb`](@ref),
    [`HardSphere`](@ref),
    [`Buckingham`](@ref),
    [`NoInteraction`](@ref),
    [`Monomial`](@ref),
    [`Exponential`](@ref)
    """
    @enum InteractionKind begin
        HardSphere
        CoulombEwaldDirect
        Coulomb
        LennardJones
        Buckingham
        Monomial
        Exponential
        UndefinedInteraction
        NoInteraction
    end

    """
        LennardJones

    Lennard-Jones `InteractionKind` between two bodies. Equals `4ε((r/σ)^12 - (r/σ)^6)`

    ## Parameters
    - `ε` in K
    - `σ` in Å
    """
    LennardJones

    """
        Coulomb

    Coulombic `InteractionKind` between two point charges. Equals `α × q₁ q₂ / r` with
    `α = e²/(4πε₀) ≈ 167101.08 atomic units`.

    ## Parameters
    - `q₁` in elementary charges.
    - `q₂` in elementary charges (if unspecified, defaults to `q₁`)
    """
    Coulomb

    """
        CoulombEwaldDirect

    Direct (i.e. real-space) contribution to the Ewald summation to take into account
    coulombic interaction with no cutoff. Equals `q₁ q₂ × erfc(α × r) / r`.

    ## Parameters
    - `α` in Å⁻¹ only depends on the required precision `p` for the Ewald summation and the
      direct cutoff `c` used. Its value is `sqrt(abs(log(c*p*t)))/c` with
      `t = sqrt(abs(log(c*p)))`
    - `q₁` in elementary charges.
    - `q₂` in elementary charges (if unspecified, defaults to `q₁`)
    """
    CoulombEwaldDirect

    """
        HardSphere

    Hard sphere interaction between two bodies. Equals `+∞` if `r < R₁ + R₂`, `0` else

    ## Parameter
    - `R₁` in Å
    - `R₂` in Å (if unspecified, defaults to `R₂`)
    """
    HardSphere

    """
        Buckingham

    Buckingham `InteractionKind` between two bodies. Equals `A exp(-B r) - C/r^6`

    ## Parameters
    - `A` in K
    - `B` in Å⁻¹
    - `C` in K⋅Å⁶
    """
    Buckingham

    """
        Monomial

    Monomial function useful as building block for constructing more complex interactions
    with an [`InteractionRuleSum`](@ref). Equals `A/r^n`

    ## Parameters
    - `A` in K⋅Åⁿ
    - `n` (no unit)

    !!! warning
        Check that you are not constructing an interaction rule which is already provided
        like [`LennardJones`](@ref) or [`Buckingham`](@ref) to ensure appropriate mixing
        rules behaviour.
    """
    Monomial

    """
        Exponential

    Exponential function useful as building block for constructing more complex interactions
    with an [`InteractionRuleSum`](@ref). Equals `A exp(-B r)`

    ## Parameters
    - `A` in K
    - `B` in Å⁻¹

    !!! warning
        Check that you are not constructing an interaction rule which is already provided
        like [`Buckingham`](@ref) to ensure appropriate mixing rules behaviour.
    """
    Exponential

    """
        NoInteraction

    Represents an explicit absence of interaction between two species.

    Attempting to combine it with an other [`InteractionKind`](@ref) through an
    [`InteractionRuleSum`](@ref) will result in an error
    """
    NoInteraction


    """
        MixingRule

    Mixing rule used when building a [`ForceField`](@ref).

    Current list:
    [`LorentzBerthelot`](@ref),
    [`WaldmanHagler`](@ref),
    [`Geometric`](@ref),
    [`IgnoreInteraction`](@ref),
    [`ErrorOnMix`](@ref)
    """
    @enum MixingRule begin
        LorentzBerthelot
        WaldmanHagler
        Geometric
        IgnoreInteraction
        ErrorOnMix
    end

    """
        LorentzBerthelot

    Lorentz-Berthelot mixing rule between two [`LennardJones`](@ref) (or in general two [`Mie`](@ref))
    interaction potentials, defined by `εᵢⱼ = √(εᵢ εⱼ)` and `σᵢⱼ = (σᵢ + σⱼ)/2`
    """
    LorentzBerthelot

    """
        WaldmanHagler

    Waldman-Hagler mixing rule between two [`LennardJones`](@ref) interaction potentials,
    defined by `εᵢⱼ = √(εᵢ εⱼ) × 2σᵢ³σⱼ³/(σᵢ⁶ + σⱼ⁶)` and `σᵢⱼ = ((σᵢ⁶ + σⱼ⁶)/2)^(1/6)`
    """
    WaldmanHagler

    """
        Geometric

    Good-Hope / Jorgensen mixing rule between two [`LennardJones`](@ref) (or in general two
    [`Mie`](@ref)) potentials, defined by `εᵢⱼ = √(εᵢ εⱼ)` and `σᵢⱼ = √(σᵢσⱼ)`
    """
    Geometric

    """
        IgnoreInteraction

    Using this mixing rule specifies that any pair of species lacking an explicit
    interaction will be attributed a [`NoInteraction`](@ref) rule.
    """
    IgnoreInteraction

    """
        ErrorOnMix

    Using this mixing rule specifies that all species pair should be explicity given an
    interaction rule. Failing this condition will result in an expected error.
    """
    ErrorOnMix
end

"""
    InteractionRule

Specific pairwise interaction rule defined by an [`InteractionKind`](@ref FF.InteractionKind)
and a list of parameters.

Construct one by calling the `InteractionKind` with the appropriate arguments. The nature
and number of required and optional arguments are documented for each `InteractionKind`.

Calling the rule with a pair distance computes the interaction energy between two particles
for this rule at this distance (in Å).

Multiple interaction rules can be combined to form more complex rules through an [`InteractionRuleSum`](@ref).

## Example
```jldoctest
julia> σ = 124.07u"K"; ε = 0.338u"nm" # Argon LJ parameters from doi:10.1021/jp803753h

julia> lj_argon = FF.LennardJones(σ, ε)
LennardJones(124.07, 3.38)

julia> lj_argon == FF.LennardJones(124.07, 3.38) # default units are documented: here K and Å
true

julia> typeof(lj_argon)
InteractionRule

julia> lj_argon(4.5u"Å) # energy in K
-73.11312068282488
```
"""
struct InteractionRule
    kind::FF.InteractionKind
    params::Vector{Float64} # units are documented for each interaction kind
    shift::Float64 # in K
    tailcorrection::Bool
end
InteractionRule(kind::FF.InteractionKind, params::Vector{Float64}) = InteractionRule(kind, params, 0.0, true)
function Base.show(io::IO, ik::InteractionRule)
    print(io, ik.kind, '(')
    join(io, ik.params, ", ")
    print(io, ')')
    if ik.kind !== FF.NoInteraction && ik.kind !== FF.UndefinedInteraction &&
       ik.kind !== FF.Coulomb && ik.kind !== FF.CoulombEwaldDirect && get(io, :inforcefield, false)
        if ik.tailcorrection
            print(io, 'ᵀ')
        end
        if ik.shift != 0.0
            print(io, 'ₛ')
        end
    end
end
Base.:(==)(a::InteractionRule, b::InteractionRule) = a.kind == b.kind && a.params == b.params
Base.isless(a::InteractionRule, b::InteractionRule) = isless(a.kind, b.kind) || (a.kind === b.kind && isless(a.params, b.params))

struct InvalidParameterNumber <: Exception
    ik::FF.InteractionKind
    got::Int
    expected::Vector{Int}
end
function Base.showerror(io::IO, x::InvalidParameterNumber)
    print(io, "Cannot construct a ", x.ik, " interaction rule from ", x.got, " parameters: expected ")
    join(io, x.expected, ", ", " or ")
    print(io, " parameters.")
end
struct UndefinedInteractionError <: Exception end
Base.showerror(io::IO, ::UndefinedInteractionError) = print(io, "Undefined interaction")

macro convertifnotfloat(unit, x)
    val = esc(x)
    quote
        Float64($val isa Unitful.AbstractQuantity ? NoUnits($val/$unit) : $val)
    end
end

function (ik::FF.InteractionKind)(args...)
    n = length(args)
    if ik === FF.LennardJones
        n == 2 || throw(InvalidParameterNumber(ik, n, [2]))
        ε = @convertifnotfloat u"K" args[1]
        σ = @convertifnotfloat u"Å" args[2]
        InteractionRule(ik, [ε, σ])
    elseif ik === FF.Monomial
        n == 2 || throw(InvalidParameterNumber(ik, n, [2]))
        InteractionRule(ik, [Float64(x) for x in args])
    elseif ik === FF.Exponential
        n == 2 || throw(InvalidParameterNumber(ik, n, [2]))
        Aexp = @convertifnotfloat u"K" args[1]
        Bexp = @convertifnotfloat u"Å^-1" args[2]
        InteractionRule(ik, [Aexp, Bexp])
    elseif ik === FF.CoulombEwaldDirect
        n > 1 || throw(InvalidParameterNumber(ik, n, [1, 2]))
        α = args[1]
        rcoulombdirect = if n == 3
            c1 = @convertifnotfloat u"e_au" args[2]
            c2 = @convertifnotfloat u"e_au" args[3]
            InteractionRule(ik, [α, c1, c2])
        elseif n == 2
            conly = @convertifnotfloat u"e_au" args[2]
            InteractionRule(ik, [α, conly, conly])
        else
            throw(InvalidParameterNumber(ik, n, [1, 2]))
        end
        if rcoulombdirect.params[2] == 0.0 || rcoulombdirect.params[3] == 0.0
            InteractionRule(FF.NoInteraction, Float64[])
        else
            rcoulombdirect
        end
    elseif ik === FF.Coulomb
        rcoulomb = if n == 2
            q1 = @convertifnotfloat u"e_au" args[1]
            q2 = @convertifnotfloat u"e_au" args[2]
            InteractionRule(ik, [q1, q2])
        elseif n == 1
            q1only = @convertifnotfloat u"e_au" args[1]
            InteractionRule(ik, [q1only, q1only])
        else
            throw(InvalidParameterNumber(ik, n, [1, 2]))
        end
        if rcoulomb.params[1] == 0.0 || rcoulomb.params[2] == 0.0
            InteractionRule(FF.NoInteraction, Float64[])
        else
            rcoulomb
        end
    elseif ik === FF.HardSphere
        if n == 2
            r1 = @convertifnotfloat u"Å" args[1]
            r2 = @convertifnotfloat u"Å" args[2]
            InteractionRule(ik, [r1, r2])
        elseif n == 1
            r1only = @convertifnotfloat u"Å" args[1]
            InteractionRule(ik, [r1only, r1only])
        else
            throw(InvalidParameterNumber(ik, n, [1, 2]))
        end
    elseif ik === FF.Buckingham
        n == 3 || throw(InvalidParameterNumber(ik, n, [3]))
        Abuck = @convertifnotfloat u"K" args[1]
        Bbuck = @convertifnotfloat u"Å^-1" args[2]
        Cbuck = @convertifnotfloat u"Å^6*K" args[3]
        InteractionRule(ik, [Abuck, Bbuck, Cbuck])
    elseif ik === FF.NoInteraction || ik === FF.UndefinedInteraction
        n == 0 || throw(InvalidParameterNumber(ik, n, [0]))
        InteractionRule(ik, Float64[])
    else
        @assert false # logically impossible
    end
end

function _shifted_InteractionRule(kind::FF.InteractionKind, params::Vector{Float64}, cutoff, tailcorrection::Bool=false)
    cut = @convertifnotfloat u"Å" cutoff
    rule = InteractionRule(kind, params, 0.0, false)
    InteractionRule(kind, params, rule(cut), tailcorrection)
end
function _shifted_InteractionRule(rule::InteractionRule, cutoff)
    _shifted_InteractionRule(rule.kind, rule.params, cutoff, rule.tailcorrection)
end

function InteractionRule(kind::FF.InteractionKind, params::Vector{Float64}, shift::Bool, cutoff, tailcorrection::Bool=!shift)
    if shift
        _shifted_InteractionRule(kind, params, cutoff, tailcorrection)
    else
        InteractionRule(kind, params, 0.0, tailcorrection)
    end
end

function (rule::InteractionRule)(r)
    if rule.kind === FF.LennardJones
        x6lj = (rule.params[2]/r)^6
        4*rule.params[1]*x6lj*(x6lj - 1)
    elseif rule.kind === FF.CoulombEwaldDirect
        NoUnits(COULOMBIC_CONVERSION_FACTOR*ENERGY_TO_KELVIN/u"K")*rule.params[2]*rule.params[3]*erfc(rule.params[1]*r)/r
    elseif rule.kind === FF.Coulomb
        NoUnits(COULOMBIC_CONVERSION_FACTOR*ENERGY_TO_KELVIN/u"K")*rule.params[1]*rule.params[2]/r
    elseif rule.kind === FF.HardSphere
        ifelse(r < rule.params[1] + rule.params[2], Inf, 0.0)
    elseif rule.kind === FF.Buckingham
        rule.params[1]*exp(-rule.params[2]*r) - rule.params[3]/(r^6)
    elseif rule.kind === FF.NoInteraction
        0.0
    elseif rule.kind === FF.Monomial
        rule.params[1]/r^rule.params[2]
    elseif rule.kind === FF.Exponential
        rule.params[1]*exp(-rule.params[2]*r)
    elseif rule.kind === FF.UndefinedInteraction
        throw(UndefinedInteractionError())
    else
        @assert false # logically impossible
    end - rule.shift
end

(rule::InteractionRule)(input::Quantity{T,Unitful.𝐋,U} where {T,U}) = rule(NoUnits(input/u"Å"))*u"K"

function (rule::InteractionRule)(input::Quantity{T,Unitful.𝐋^2,U} where {T,U})
    r2 = NoUnits(input/u"Å^2")
    (if rule.kind === FF.LennardJones
        σ2 = rule.params[2]^2
        x6lj = (σ2/r2)^3
        4*rule.params[1]*x6lj*(x6lj - 1)
    elseif rule.kind === FF.HardSphere
        ifelse(r2 < (rule.params[1] + rule.params[2])^2, Inf, 0.0)
    elseif rule.kind === FF.NoInteraction
        0.0
    elseif rule.kind === FF.Monomial
        rule.params[1]/r2^(rule.params[2]/2)
    else
        return rule(sqrt(r2))*u"K"
    end - rule.shift)*u"K"
end

function tailcorrection(rule::InteractionRule, cutoff)
    cut = @convertifnotfloat u"Å" cutoff
    (if !rule.tailcorrection || rule.kind === FF.NoInteraction || rule.kind === FF.HardSphere || rule.kind === FF.CoulombEwaldDirect || isinf(cutoff)
        0.0
    elseif rule.kind === FF.LennardJones
        ε, σ = rule.params
        xlj3 = (σ/cut)^3
        xlj9 = xlj3^3
        (4/3)*ε*σ^3*(xlj9/3 - xlj3)
    elseif rule.kind === FF.Buckingham
        A, B, C = rule.params
        A*exp(-B*cut)*(2.0 + B*cut*(2.0 + B*cut))/B^3 - C/(3*cut^3)
    elseif rule.kind === FF.Coulomb
        error("Coulomb direct pair interaction cannot have a tail correction: use Ewald summation.")
    elseif rule.kind === FF.UndefinedInteraction
        throw(UndefinedInteractionError())
    else
        @assert false
    end)*u"K"
end

function derivativesGrid(rule::InteractionRule, d2)
    r2 = @convertifnotfloat u"Å^2" d2
    if rule.kind === FF.LennardJones
        ε, σ = rule.params
        x6 = (σ^2/r2)^3
        r4 = r2*r2
        value = 4*ε*x6*(x6 - 1)
        ∂1 = 24*ε*(x6*(1 - 2*x6))/r2
        ∂2 = 96*ε*(x6*(7*x6 - 2))/r4
        ∂3 = 384*ε*(x6*(5 - 28*x6))/(r4*r4)
    elseif rule.kind === FF.Coulomb
        error("Coulomb interactions should not be taken into account in VdW grids.")
    elseif rule.kind === FF.HardSphere
        value = ifelse(r2 < (rule.params[1] + rule.params[2])^2, Inf, 0.0)
        ∂1 = ∂2 = ∂3 = 0.0
    elseif rule.kind === FF.Buckingham
        A, B, C = rule.params
        r4 = r2*r2
        r = sqrt(r2)
        r6 = r4*r2
        x6 = C/r6
        xe = A*exp(-B*r)
        value = xe - x6
        ∂1 = -B*xe/r + 6*x6/r2
        ∂2 = -48*x6/r4 + B*xe*(1+B*r)/(r2*r)
        ∂3 = -(3*B*r + B*B*r2 + 3)*B*xe*r/r6 + 480*C/(r6*r6)
    elseif rule.kind === FF.NoInteraction || rule.kind == FF.CoulombEwaldDirect
        # CoulombEwaldDirect is not included in the VdW grid but directly taken into
        # account in the Ewald grid
        return 0.0, 0.0, 0.0, 0.0
    elseif rule.kind === FF.Monomial
        error("VdW grid not implemented for Monomial")
    elseif rule.kind === FF.Exponential
        error("VdW grid not implemented for Exponential")
    elseif rule.kind === FF.UndefinedInteraction
        throw(UndefinedInteractionError())
    else
        @assert false # logically impossible
    end
    return value-rule.shift, ∂1, ∂2, ∂3
end


struct RepeatedRuleKind <: Exception
    x::InteractionRule
    y::InteractionRule
end
function Base.showerror(io::IO, e::RepeatedRuleKind)
    print(io, "Defining an `InteractionRuleSum` of two rules with the same kind is forbidden: attempted to sum ",
              e.x, " with ", e.y)
end

"""
    InteractionRuleSum

Represents the sum of multiple interactions between the same pair of species.

!!! note
    Attempting to sum any interaction with a [`NoInteraction`](@ref) will yield in an error.
    It is also forbidden to put a single `InteractionRule` in the argument list of an
    [`InteractionRuleSum`](@ref).

## Example
```jldoctest
julia> buck = FF.Buckingham(5.581e7, 3.9850, 9.167e5)
Buckingham(5.581e7, 3.985, 916700.0)

julia> buck(2.6) # Buckingham potential
-1201.4909774834446

julia> buck(0.85) # Unphysical attraction of the Buckingham potential at short distance
-544139.1256276625

julia> inter = InteractionRuleSum([FF.HardSphere(1.5), buck])
InteractionRuleSum([HardSphere(0.6), Buckingham(5.581e7, 3.985, 916700.0))

julia> inter(2.6) # Buckingham interaction in the physical part
-1201.4909774834446

julia> inter(0.85) # a hard sphere model corrects the Buckingham potential at short distance
Inf
```
"""
struct InteractionRuleSum
    rules::Vector{InteractionRule}

    function InteractionRuleSum(rules::Vector{InteractionRule})
        isempty(rules) && error("`InteractionRuleSum(InteractionRule[])` is ill-defined. Use `FF.NoInteraction` or `FF.UndefinedInteraction` if need be.")
        length(rules) == 1 && error("Defining `InteractionRuleSum([rule]) is forbidden. Directly use `rule` instead.")
        srules = sort(rules)
        srules[end].kind == FF.NoInteraction && error("Summing any rule with a `FF.NoInteraction` is forbidden. Sum the other rules instead.")
        srules[end].kind == FF.UndefinedInteraction && error("Attempting to sum a rule with a `FF.UndefinedInteraction`.")
        kinds = [rule.kind for rule in srules]
        unique!(kinds)
        if length(kinds) != length(srules)
            @assert length(kinds) < length(srules)
            idx = length(srules)
            for (i, (kind, rule)) in enumerate(zip(kinds, srules))
                if kind != rule.kind
                    idx = i
                    break
                end
            end
            throw(RepeatedRuleKind(srules[idx-1], srules[idx]))
        end
        new(srules)
    end
end
function Base.show(io::IO, f::InteractionRuleSum)
    print(io, typeof(f), '(', '[')
    join(io, f.rules, ", ")
    print(io, ')')
end

function (f::InteractionRuleSum)(r)
    ret = r isa Quantity ? 0.0u"K" : 0.0
    for x in f.rules
        ret += x(r)
    end
    ret
end

function derivativesGrid(f::InteractionRuleSum, d2)
    r2 = @convertifnotfloat u"Å^2" d2
    value = ∂1 = ∂2 = ∂3 = 0.0
    for x in f.rules
        v, p1, p2, p3 = derivativesGrid(x, r2)
        value += v
        ∂1 += p1
        ∂2 += p2
        ∂3 += p3
    end
    value, ∂1, ∂2, ∂3
end

function sum_one_rule(l::Vector{InteractionRule}, r::InteractionRule)
    if r.kind === FF.NoInteraction
        l
    elseif r.kind === FF.UndefinedInteraction
        [r]
    else
        vcat(l, r)
    end
end
function sum_rules(r1, r2)
    if r1 isa InteractionRuleSum
        if r2 isa InteractionRuleSum
            InteractionRuleSum(vcat(r1.rules, r2.rules))
        else
            InteractionRuleSum(sum_one_rule(r1.rules, r2))
        end
    elseif r2 isa InteractionRuleSum
        InteractionRuleSum(sum_one_rule(r2.rules, r1))
    else
        if r1.kind === FF.NoInteraction
            r2
        elseif r1.kind === FF.UndefinedInteraction
            r1
        else
            l = sum_one_rule([r1], r2)
            if length(l) == 1
                l[1]
            else
                InteractionRuleSum(l)
            end
        end
    end
end

function tailcorrection(f::InteractionRuleSum, cutoff)
    ret = 0.0u"K"
    for x in f.rules
        ret += tailcorrection(x, cutoff)
    end
    ret
end
