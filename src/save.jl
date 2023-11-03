using CodecZstd, TranscodingStreams

function _save_init(io::IO, step::SimulationStep{N,T}) where {N,T}
    s = Serializer(io)
    Serialization.writeheader(s)
    serialize(s, SimulationStep{N,T})
    serialize(s, step.ff)
    serialize(s, step.charges)
    serialize(s, step.isrigid)
    serialize(s, step.ffidx)
    serialize(s, step.parallel)
    nothing
end

function save_init(path, step::SimulationStep)
    if splitext(path)[2] == ".zst"
        open(ZstdCompressorStream, path, "w") do io
            _save_init(io, step)
        end
    else
        open(path, "w") do io
            _save_init(io, step)
        end
    end
end

function write_vector(io::IO, x::Vector{T}) where T
    bytes = 0
    bytes += write(io, length(x))
    if T <: Vector
        for a in x
            bytes += write_vector(io, a)
        end
    else
        bytes += write(io, x)
    end
    bytes
end

function read_vector(io::IO, ::Type{Vector{T}}) where T
    l = read(io, Int)
    x = Vector{T}(undef, l)
    if T <: Vector
        for i in 1:l
            x[i] = read_vector(io, T)
        end
    else
        read!(io, x)
    end
    x
end


function save(path::AbstractString, step::SimulationStep)
    if splitext(path)[2] == ".zst"
        open(ZstdCompressorStream, path, "a") do io
            write(io, step.mat) + write_vector(io, step.atoms) + write_vector(io, step.positions)
        end
    else
        open(path, "a") do io
            write(io, step.mat) + write_vector(io, step.atoms) + write_vector(io, step.positions)
        end
    end
end


struct StreamSimulationStep{N,T}
    io::Union{IOStream,TranscodingStream{ZstdDecompressor,IOStream}}
    ff::ForceField
    charges::Vector{Te_au}
    isrigid::BitVector
    ffidx::Vector{Vector{Int}}
    parallel::Bool
    runnable::Bool
end

function _load_init(io::IO, ::Type{SimulationStep{N,T}}; runnable=true) where {N,T}
    s = Serializer(io)
    deserialize(s) == SimulationStep{N,T} || error(lazy"$path does not correspond to a SimulationStep")
    ff = deserialize(s)::ForceField
    charges = deserialize(s)::Vector{Te_au}
    isrigid = deserialize(s)::BitVector
    ffidx = deserialize(s)::Vector{Vector{Int}}
    parallel = deserialize(s)::Bool
    StreamSimulationStep{N,T}(io, ff, charges, isrigid, ffidx, parallel, runnable)
end

function load_init(path, ::Type{SimulationStep{N,T}}; runnable=true) where {N,T}
    if splitext(path)[2] == ".zst"
        _load_init(ZstdDecompressorStream(open(path)), SimulationStep{N,T}; runnable)
    else
        _load_init(open(path), SimulationStep{N,T}; runnable)
    end
end

function StreamSimulationStep(path; runnable=true)
    T = if splitext(path)[2] == ".zst"
        open(ZstdDecompressorStream, path) do io
            deserialize(io)
        end
    else
        deserialize(path)
    end
    load_init(path, T; runnable)
end


function _load(io::IO, stream::StreamSimulationStep{N,T}) where {N,T}
    mat = read(io, SMatrix{3,3,TÅ,9})
    atoms = read_vector(io, Vector{Tuple{Int,Int,Int}})
    positions = read_vector(io, Vector{SVector{3,TÅ}})
    if stream.runnable
        # the position of the k-th atom in the j-th system of kind i is positions[l] where
        # atoms[l] == (i,j,k) and posidx[i][j][k] == l
        posidx = [Vector{Int}[] for _ in stream.ffidx]
        not_in_freespecies = [BitSet() for _ in stream.ffidx]
        for (l, (i,j,k)) in enumerate(atoms)
            I = posidx[i]
            if length(I) < j
                lenJ = length(stream.ffidx[i])
                append!(I, Vector{Int}(undef, lenJ) for _ in j-length(I))
            end
            J = I[j]
            J[k] = l
            push!(not_in_freespecies[i], j)
        end
        freespecies = [setdiff(1:length(posidxi), not_in_i) for (posidxi, not_in_i) in zip(posidx, not_in_freespecies)]
    else
        posidx = Vector{Vector{Vector{Int}}}(undef, length(stream.ffidx))
        freespecies = Vector{Vector{Int}}(undef, length(stream.ffidx))
    end
    psystem = PeriodicSystem(; xpositions=positions,
                               ypositions=SVector{3,TÅ}[],
                               unitcell=mat,
                               parallel=stream.parallel,
                               cutoff=stream.ff.cutoff, output=0.0u"K")
    eof(io) && close(io)
    SimulationStep{N,T}(stream.ff, stream.charges, psystem, atoms, posidx, freespecies, stream.isrigid, stream.ffidx)
end
load(stream::StreamSimulationStep) = _load(stream.io, stream)

function Base.iterate(stream::StreamSimulationStep, _=nothing)
    isopen(stream.io) || return nothing
    return load(stream), nothing
end

Base.IteratorSize(::Type{StreamSimulationStep{N,T}}) where {N,T} = Base.SizeUnknown()
Base.eltype(::Type{StreamSimulationStep{N,T}}) where {N,T} = SimulationStep{N,T}
Base.isdone(stream::StreamSimulationStep, _=nothing) = !isopen(stream.io)
