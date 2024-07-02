const mcmovenames = (
    :translation,
    :rotation,
    :random_translation,
    :random_rotation,
    :random_reinsertion,
    :swap,
)

"""
    MCMoves

List of MC moves and their respective probability.
The possible moves are:
- `:translation`: displacement along a random direction with a uniformly chosen distance. The
  maximum distance is adjusted at the end of each cycle to have a average acceptance rate
  around 0.5.
- `:rotation`: rotation around the reference atom of a random angle chosen uniformly. The
  maximum angle is adjusted at the end of each cycle to have a average acceptance rate
  around 0.5.
- `:random_translation`: displacement along a random direction. The distance is chosen at
  random between 0 and the maximal distance between two atoms in the unit cell.
- `:random_rotation`: rotation around the reference atom of a random angle. The angle is
  chosen at random between 0° and 180° along each axis.
- `:random_reinsertion`: combination of `:random_translation` and `:random_rotation`,
  effectively equivalent to removing the molecule and trying to insert it again somewhere.
- `:swap`: either deletion of the molecule or insertion of a similar molecule at a random
  placement in the system. Swap moves will also occur even if there are no molecule of that
  kind in the system, to allow insertion of molecules from an empty starting point.

An `MCMoves` object can be defined by giving the name of the MC moves and their
corresponding probabilities as keyword arguments (the probabilities) will be normalized.
For instance:
```julia
julia> x = MCMoves(; translation=2, random_rotation=0.5, random_reinsertion=2.5)
MCMoves(; translation=0.4, random_rotation=0.1, random_reinsertion=0.5)
```

Such an object `x::MCMoves` can be called with a random number sampled in `[0,1)` to yield
the name of a move according to its probability. For example:
```julia
julia> x(rand())
:translation

julia> x(rand())
:random_reinsertion

julia> x(rand())
:random_reinsertion

julia> x(rand())
:translation

julia> x(rand())
:random_reinsertion
```
"""
struct MCMoves
    # the stored values are the cumulative probability
    cumulatives::NTuple{length(mcmovenames)-1,Float64}
end
function MCMoves(monoatomic::Bool)
    MCMoves(if monoatomic
        (0.5, 0.5, 1.0, 1.0, 1.0) # 50% translation, 50% random translation
    else
        (0.33, 0.66, 0.66, 0.66, 1.0) # 33% translation, 33% rotation, 34% random reinsertion
    end)
end

mutable struct AccumulateKwargs{T}
    const tot::Float64
    const kwargs::T
    last_i::Int
    setkwargs::Int
    c::Float64
end
AccumulateKwargs(tot, kwargs) = AccumulateKwargs(tot, kwargs, 0, 0, 0.0)
function (acc::AccumulateKwargs)(i::Int)
    @assert (acc.last_i += 1) == i
    v = Float64(get(acc.kwargs, @inbounds(mcmovenames[i]), 0.0)/acc.tot)::Float64
    acc.setkwargs += v != 0.0
    acc.c += v
end
function MCMoves(; kwargs...)
    tot = Float64(sum(last, kwargs; init=0.0))::Float64
    acc = AccumulateKwargs(tot, kwargs)
    ret = MCMoves(ntuple(acc, Val(length(mcmovenames)-1)))
    @assert acc.c ≤ 1.0
    if acc.setkwargs + haskey(kwargs, last(mcmovenames)) != length(kwargs)
        wrongkwargs = filter(x -> first(x) ∉ mcmovenames, kwargs)
        error(lazy"""The following kwarg$(length(wrongkwargs)==1 ? " is" : "s are") invalid for MCMoves: \"$(join([first(x) for x in wrongkwargs], "\\\", \\\"", "\\\" and \\\""))\".
        Choose among \"$(join(mcmovenames, "\\\", \\\"", "\\\" and \\\""))\"""")
    end
    ret
end
function Base.show(io::IO, m::MCMoves)
    lastcumul = 0.0
    print(io, MCMoves, "(; ")
    for (i, c) in enumerate(m.cumulatives)
        c == lastcumul && continue
        lastcumul == 0.0 || print(io, ", ")
        @printf io "%s=%.9g" mcmovenames[i] (c - lastcumul)
        lastcumul = c
    end
    if !isapprox(lastcumul, 1)
        lastcumul == 0.0 || print(io, ", ")
        @printf io "%s=%.9g" last(mcmovenames) (1 - lastcumul)
    end
    print(io, ')')
end

# Sample a move at random
function (m::MCMoves)(r::Float64)
    for (i, c) in enumerate(m.cumulatives)
        r < c && return mcmovenames[i]
    end
    return last(mcmovenames)
end

@inline function Base.getindex(m::MCMoves, s::Symbol)
    @inbounds if s === first(mcmovenames)
        first(m.cumulatives)
    elseif s === last(mcmovenames)
        1.0 - last(m.cumulatives)
    else
        i = findfirst(==(s), mcmovenames)
        if i isa Nothing
            error(lazy"No MC move corresponding to name $s")
        else
            m.cumulatives[i] - m.cumulatives[i-1]
        end
    end
end

Base.:(==)(m1::MCMoves, m2::MCMoves) = m1.cumulatives == m2.cumulatives
Base.hash(m::MCMoves) = hash(m.cumulatives, hash(MCMoves))


function random_translation(rng, positions::AbstractVector{SVector{3,TÅ}}, dmax_or_mat::Union{TÅ,AbstractMatrix})
    r = if dmax_or_mat isa TÅ
        SVector{3}(((2*rand(rng)-1)*dmax_or_mat) for _ in 1:3)
    else
        dmax_or_mat*(rand(rng, SVector{3,Float64}) .- 0.5)
    end
    [poss + r for poss in positions]
end
function random_rotation(rng, positions::AbstractVector{SVector{3,TÅ}}, θmax, bead, _r=nothing)
    bead == 0 && return positions
    θ = θmax*(2*rand(rng)-1)
    s, c = sincos(θ)
    r = _r isa Nothing ? rand(rng, 1:3) : _r
    mat = if r == 1
        SMatrix{3,3}(1, 0, 0, 0, c, s, 0, -s, c)
    elseif r == 2
        SMatrix{3,3}(c, 0, -s, 0, 1, 0, s, 0, c)
    else
        SMatrix{3,3}(c, s, 0, -s, c, 0, 0, 0, 1)
    end
    refpos = positions[bead]
    [refpos + mat*(poss - refpos) for poss in positions]
end
