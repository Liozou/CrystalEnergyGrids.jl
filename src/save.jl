using CodecZstd, TranscodingStreams

function _save_init(io::IO, step::SimulationStep)
    s = Serializer(io)
    Serialization.writeheader(s)
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


"""
    ProtoSimulationStep

Structure that contains the same information as a [`SimulationStep`](@ref) but organized to
be fast to create and read, not to be used in a simulation.

All fields documented for a [`SimulationStep`](@ref) are accessible on a
`ProtoSimulationStep` except `freespecies` and `posidx`, by assuming that
`freespecies` is empty and `posidx` is the reciprocal function to `atoms`.

Use `SimulationStep(x::ProtoSimulationStep)` to construct a `SimulationStep` from it.
"""
struct ProtoSimulationStep
    ff::ForceField
    charges::Vector{Te_au}
    mat::SMatrix{3,3,TÅ,9}
    positions::Vector{SVector{3,TÅ}}
    parallel::Bool
    atoms::Vector{NTuple{3,Int}}
    isrigid::BitVector
    ffidx::Vector{Vector{Int}}
end

"""
"""
function SimulationStep(x::ProtoSimulationStep)
    posidx = [Vector{Int}[] for _ in x.ffidx]
    not_in_freespecies = [BitSet() for _ in x.ffidx]
    for (l, (i,j,k)) in enumerate(x.atoms)
        I = posidx[i]
        if length(I) < j
            lenJ = length(x.ffidx[i])
            append!(I, Vector{Int}(undef, lenJ) for _ in j-length(I))
        end
        J = I[j]
        J[k] = l
        push!(not_in_freespecies[i], j)
    end
    freespecies = [setdiff(1:length(posidxi), not_in_i) for (posidxi, not_in_i) in zip(posidx, not_in_freespecies)]
    psystem = PeriodicSystem(; xpositions=x.positions,
                            ypositions=SVector{3,TÅ}[],
                            unitcell=x.mat,
                            parallel=x.parallel,
                            cutoff=x.ff.cutoff, output=0.0u"K")
    SimulationStep(x.ff, x.charges, psystem, x.atoms, posidx, freespecies, x.isrigid, x.ffidx)
end

"""
    StreamSimulationStep

Iterable reader of a "steps.stream" or "steps.zst" output file.

Each iteration returns a [`ProtoSimulationStep`](@ref).
"""
struct StreamSimulationStep
    io::Union{IOStream,TranscodingStream{ZstdDecompressor,IOStream}}
    ff::ForceField
    charges::Vector{Te_au}
    isrigid::BitVector
    ffidx::Vector{Vector{Int}}
    parallel::Bool
end

function _load_init(io::IO)
    s = Serializer(io)
    ff = deserialize(s)::ForceField
    charges = deserialize(s)::Vector{Te_au}
    isrigid = deserialize(s)::BitVector
    ffidx = deserialize(s)::Vector{Vector{Int}}
    parallel = deserialize(s)::Bool
    StreamSimulationStep(io, ff, charges, isrigid, ffidx, parallel)
end

function StreamSimulationStep(path)
    if splitext(path)[2] == ".zst"
        _load_init(ZstdDecompressorStream(open(path)))
    else
        _load_init(open(path))
    end
end


function _load(io::IO, stream::StreamSimulationStep)
    mat = read(io, SMatrix{3,3,TÅ,9})
    atoms = read_vector(io, Vector{Tuple{Int,Int,Int}})
    positions = read_vector(io, Vector{SVector{3,TÅ}})
    eof(io) && close(io)
    ProtoSimulationStep(stream.ff, stream.charges, mat, positions, stream.parallel, atoms, stream.isrigid, stream.ffidx)
end
load(stream::StreamSimulationStep) = _load(stream.io, stream)

function Base.iterate(stream::StreamSimulationStep, _=nothing)
    isopen(stream.io) || return nothing
    return load(stream), nothing
end

Base.IteratorSize(::Type{StreamSimulationStep}) = Base.SizeUnknown()
Base.eltype(::Type{StreamSimulationStep})= ProtoSimulationStep
Base.isdone(stream::StreamSimulationStep, _=nothing) = eof(stream.io)

"""
    stream_to_pdb(dir::AbstractString)

Takes either "steps.serial" or "steps.zst" for the given folder `dir` and
create an `output.pdb` file in the same folder containing the same information.
"""
function stream_to_pdb(dir::AbstractString)
    trajectory = joinpath(dir, "trajectory.pdb")
    if isfile(joinpath(dir, "trajectory.pdb"))
        error(lazy"$(joinpath(dir, \"trajectory.pdb\")) already exists.")
    end
    _serialfile = joinpath(dir, "steps.serial")
    stream = StreamSimulationStep(isfile(_serialfile) ? _serialfile : joinpath(dir, "steps.zst"))
    atomcounter = Counter3D()
    for (i, step) in enumerate(stream)
        lengths, angles = cell_parameters(step.mat)
        output_pdb(trajectory, step, lengths, angles, i, atomcounter)
    end
    nothing
end
