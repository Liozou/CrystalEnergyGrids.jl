module CrystalEnergyGrids

using StaticArrays
using LinearAlgebra: LinearAlgebra, norm, det, cross, dot
using Pkg: TOML

using StaticArrays
using AtomsBase
using Scratch: @get_scratch!

export EnergyGrid, parse_grid, interpolate_grid
export CrystalEnergySetup
export energy_point, energy_grid
# other exports in files

const MODULE_VERSION = VersionNumber(TOML.parsefile(joinpath(dirname(@__DIR__), "Project.toml"))["version"])
scratchspace::String = ""
function __init__()
    global scratchspace
    scratchspace = @get_scratch!("hnc-$(MODULE_VERSION.major).$(MODULE_VERSION.minor)")
end

include("constants.jl")
include("lebedev.jl")
include("utils.jl")
include("ChangePositionSystems.jl")
include("interactions.jl")
include("forcefields.jl")
include("coordinates.jl")
include("probes.jl")
include("energy.jl")
include("energy_minimizer.jl")
include("raspa.jl")
include("ewald.jl")
include("grids.jl")
include("montecarlo.jl")
include("output.jl")
include("simulation.jl")
include("parameterinputs.jl")
include("lammpstrj_parser.jl")
include("hypernettedchain.jl")
include("rdf.jl")
include("mdft.jl")
include("viewer.jl")

# const BUFFERS = [(Vector{Cdouble}(undef, 64),Vector{Cdouble}(undef, 64)) for _ in 1:Base.Threads.nthreads()]
# const GLOBAL_LOCK = Base.Threads.ReentrantLock()
# function acquire_buffers()
#     bs = lock(GLOBAL_LOCK) do
#         if isempty(BUFFERS)
#             (Vector{Cdouble}(undef, 64),Vector{Cdouble}(undef, 64))
#         else
#             pop!(BUFFERS)
#         end
#     end
#     return bs
# end
# function release_buffers(bs)
#     lock(GLOBAL_LOCK) do
#         push!(BUFFERS, bs)
#     end
# end

end
