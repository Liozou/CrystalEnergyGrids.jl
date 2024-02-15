module CrystalEnergyGrids

using Statistics: mean
using LinearAlgebra: LinearAlgebra, norm, det, cross, dot, mul!
using Pkg: TOML
using Random: rand, default_rng, TaskLocalRNG, randexp

using StaticArrays
using AtomsBase

export EnergyGrid, parse_grid, interpolate_grid
export CrystalEnergySetup
export energy_point, energy_grid
# other exports in files

using Scratch: @get_scratch!
const MODULE_VERSION = VersionNumber(TOML.parsefile(joinpath(dirname(@__DIR__), "Project.toml"))["version"])
scratchspace::String = ""
function __init__()
    global scratchspace
    scratchspace = @get_scratch!("CEG-$(MODULE_VERSION.major).$(MODULE_VERSION.minor)")
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
include("basins.jl")
include("raspa.jl")
include("ewald.jl")
include("grids.jl")
include("mcmoves.jl")
include("tailcorrection.jl")
include("montecarlo.jl")
include("hash.jl")
include("save.jl")
include("output.jl")
include("simulation.jl")
include("parameterinputs.jl")
include("averageclusters.jl")


end
