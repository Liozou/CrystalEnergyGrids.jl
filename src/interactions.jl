export FF, InteractionRule, InteractionRuleSum

module FF
    @enum InteractionKind begin
        Coulomb
        HardSphere
        LennardJones
        Buckingham
        Monomial
        Exponential
    end

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
        HardSphere

    Hard sphere interaction between two bodies. Equals `+∞` if `r < R`, `0` else

    ## Parameter
    - `R` in Å
    """
    HardSphere

    """
        LennardJones

    Lennard-Jones `InteractionKind` between two bodies. Equals `4ε((r/σ)^12 - (r/σ)^6)`

    ## Parameters
    - `ε` in K
    - `σ` in Å
    """
    LennardJones

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

    !!! warn
        Check that you are not constructing an interaction rule which is already provided
        like [`LennardJones`](@ref) or [`Buckingham`](@ref) as those will have better
        performance.
    """
    Monomial

    """
        Exponential

    Exponential function useful as building block for constructing more complex interactions
    with an [`InteractionRuleSum`](@ref). Equals `A exp(-B r)`

    ## Parameters
    - `A` in K
    - `B` in Å⁻¹

    !!! warn
        Check that you are not constructing an interaction rule which is already provided
        like [`Buckingham`](@ref) as those will have better performance.
    """
    Exponential
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
julia> σ = 124.07; ε = 3.38 # Argon LJ parameters in K and Å from doi:10.1021/jp803753h

julia> lj_argon = FF.LennardJones(σ, ε)
LennardJones(124.07, 3.38)

julia> typeof(lj_argon)
InteractionRule

julia> lj_argon(4.5) # energy in K
-73.11312068282488
```
"""
struct InteractionRule
    kind::FF.InteractionKind
    params::Vector{Float64}
end
function Base.show(io::IO, ik::InteractionRule)
    print(io, ik.kind, '(')
    join(io, ik.params, ", ")
    print(io, ')')
end

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

function (ik::FF.InteractionKind)(args...)
    n = length(args)
    if ik === FF.Coulomb
        if n == 2
            InteractionRule(FF.Coulomb, [Float64(x) for x in args])
        elseif n == 1
            InteractionRule(FF.Coulomb, [args[1], args[1]])
        else
            throw(InvalidParameterNumber(ik, n, [1, 2]))
        end
    elseif ik === FF.HardSphere
        n == 1 || throw(InvalidParameterNumber(ik, n, [1]))
        InteractionRule(ik, [Float64(args[1])])
    elseif ik === FF.LennardJones || ik === FF.Monomial || ik === FF.Exponential
        n == 2 || throw(InvalidParameterNumber(ik, n, [2]))
        InteractionRule(ik, [Float64(x) for x in args])
    elseif ik == FF.Buckingham
        n == 3 || throw(InvalidParameterNumber(ik, n, [3]))
        InteractionRule(ik, [Float64(x) for x in args])
    else
        @assert false # logically impossible
    end
end

# (ik::FF.InteractionKind)(args...) = complete_rule(ik, [Float64(x) for x in args])
# """
#     complete_rule(::Val{ik}, args::Vector{Float64}) where {ik}

# !!! warn
#     This function should only be extended with a new method if you intend to add a new
#     [`InteractionKind`](@ref FF.InteractionKind)

# Called when constructing an [`InteractionRule`](@ref): given `ik`, an `InteractionKind` and a
# list of parameters converted to `Float64`, should return an `InteractionRule` with the same
# `InteractionKind` and the appropriate parameters. In practise, this serves as an extendable
# way of handling optional arguments.

# ## Example

# ```julia
# julia> cation = FF.Coulomb(1) # coulomb interaction between charges of +1 e
# FF.Coulomb(1.0, 1.0)

# julia> cation === complete_rule(Val(FF.Coulomb), [1.0])
# true
# ```
# """
# complete_rule(::Val{ik}, args::Vector{Float64}) where {ik} = InteractionRule(ik, args)
# function complete_rule(::Val{FF.Coulomb}, args::Vector{Float64})
#     if length(args) == 2
#         InteractionRule(FF.Coulomb, args)
#     elseif length(args) == 1
#         InteractionRule(FF.Coulomb, [args[1], args[1]])
#     else
#         error(lazy"Inappropriate number of arguments when constructing a Coulomb interaction ($(length(args)), expected 1 or 2)")
#     end
# end

function (rule::InteractionRule)(r)
    if rule.kind === FF.Coulomb
        (COULOMBIC_CONVERSION_FACTOR*ENERGY_TO_KELVIN)*rule.params[1]*rule.params[2]/r
    elseif rule.kind === FF.HardSphere
        ifelse(r < rule.params[1], Inf, 0.0)
    elseif rule.kind === FF.LennardJones
        x6lj = (rule.params[2]/r)^6
        4*rule.params[1]*x6lj*(x6lj - 1)
    elseif rule.kind === FF.Buckingham
        rule.params[1]*exp(-rule.params[2]*r) - rule.params[3]/(r^6)
    elseif rule.kind === FF.Monomial
        rule.params[1]/r^rule.params[2]
    elseif rule.kind === FF.Exponential
        rule.params[1]*exp(-rule.params[2]*r)
    else
        @assert false # logically impossible
    end
end


"""
    InteractionRuleSum

Represents the sum of multiple interactions between the same pair of species.

## Example
```jldoctest
julia> buck = FF.Buckingham(5.581e7, 3.9850, 9.167e5)
Buckingham(5.581e7, 3.985, 916700.0)

julia> buck(2.6) # Buckingham potential
-1201.4909774834446

julia> buck(0.85) # Unphysical attraction of the Buckingham potential at short distance
-544139.1256276625

julia> inter = InteractionRuleSum([FF.HardSphere(1.5), buck])
InteractionRuleSum([HardSphere(1.3), Buckingham(5.581e7, 3.985, 916700.0))

julia> inter(2.6) # Buckingham interaction in the physical part
-1201.4909774834446

julia> inter(0.85) # a hard sphere model corrects the Buckingham potential at short distance
Inf
```
"""
struct InteractionRuleSum
    rules::Vector{InteractionRule}
end
function Base.show(io::IO, f::InteractionRuleSum)
    print(io, typeof(f), '(', '[')
    join(io, f.rules, ", ")
    print(io, ')')
end
function (f::InteractionRuleSum)(r)
    ret = 0.0
    for x in f.rules
        ret += x(r)
    end
    ret
end
