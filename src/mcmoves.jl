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
    cumulatives::NTuple{4,Float64}
end
const mcmovenames = (
    :translation,
    :rotation,
    :random_translation,
    :random_rotation,
    :random_reinsertion,
)
function MCMoves(monoatomic::Bool)
    MCMoves(if monoatomic
        (0.98, 0.98, 1.0, 1.0) # 98% translation, 2% random translation
    else
        (0.49, 0.98, 0.98, 0.98) # 49% translation, 49% rotation, 2% random reinsertion
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

function (m::MCMoves)(r::Float64)
    for (i, c) in enumerate(m.cumulatives)
        r < c && return mcmovenames[i]
    end
    return last(mcmovenames)
end

Base.:(==)(m1::MCMoves, m2::MCMoves) = m1.cumulatives == m2.cumulatives
Base.hash(m::MCMoves) = hash(m.cumulatives, hash(MCMoves))


function random_translation(rng, positions::AbstractVector{SVector{3,TÅ}}, dmax::TÅ)
    r = SVector{3}(((2*rand(rng)-1)*dmax) for _ in 1:3)
    [poss + r for poss in positions]
end