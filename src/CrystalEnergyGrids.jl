module CrystalEnergyGrids

using StaticArrays
using LinearAlgebra: norm, det
using StaticArrays
using OffsetArrays
using AtomsBase

export parse_grid, interpolate_grid # other exports in files

include("constants.jl")
include("raspa.jl")
include("ewald.jl")

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


struct RaspaGrid
    spacing::Cdouble
    dims::NTuple{3,Cint}
    size::NTuple{3,Cdouble}
    shift::NTuple{3,Cdouble}
    Δ::NTuple{3,Cdouble}
    unitcell::NTuple{3,Cdouble}
    num_unitcell::NTuple{3,Cint}
    ε_Ewald::Cfloat
    grid::Array{Cfloat,4}
end
function Base.show(io::IO, ::MIME"text/plain", rg::RaspaGrid)
    print(io, rg.ε_Ewald == Inf ? "VdW" : "Coulomb", " grid with ")
    join(io, rg.dims .+ 1, '×')
    print(io, " points for a ")
    if rg.num_unitcell != (1,1,1)
        join(io, rg.num_unitcell, '×')
        print(io, " supercell of the ")
    end
    print("unit cell of size ")
    join(io, rg.unitcell, " Å × ")
    print(io, " Å")
end

macro triplet(foo)
    esc(quote
        ($foo, $foo, $foo)
    end)
end

function parse_grid(file, iscoulomb)
    open(file) do io
        spacing = read(io, Cdouble)
        dims = @triplet read(io, Cint)
        size = @triplet read(io, Cdouble)
        shift = @triplet read(io, Cdouble)
        Δ = @triplet read(io, Cdouble)
        unitcell = @triplet read(io, Cdouble)
        num_unitcell = @triplet read(io, Cint)
        ε_Ewald = iscoulomb ? read(io, Cdouble) : Inf
        grid = Array{Cfloat, 4}(undef, dims[1]+1, dims[2]+1, dims[3]+1, 8)
        read!(io, grid)
        grid .*= ENERGY_TO_KELVIN
        RaspaGrid(spacing, dims, size, shift, Δ, unitcell, num_unitcell, ε_Ewald, grid)
    end
end

function interpolate_grid(g::RaspaGrid, point)
    shifted = @. (point - g.shift)/g.size*g.dims + 1
    p0 = floor.(Int, shifted)
    p1 = p0 .+ 1
    r = shifted .- p0
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    rx, ry, rz = r
    # X, a = acquire_buffers()
    X = MVector{64,Cfloat}(undef)
    for i in 1:8
        X[i*8-7] = g.grid[z0,y0,x0,i]
        X[i*8-3] = g.grid[z1,y0,x0,i]
        X[i*8-5] = g.grid[z0,y1,x0,i]
        X[i*8-1] = g.grid[z1,y1,x0,i]
        X[i*8-6] = g.grid[z0,y0,x1,i]
        X[i*8-2] = g.grid[z1,y0,x1,i]
        X[i*8-4] = g.grid[z0,y1,x1,i]
        X[i*8  ] = g.grid[z1,y1,x1,i]
    end
    if any(>(1e10), X)
        # release_buffers((X,a))
        return 1e8
    end
    a = COEFF * X

    ret = 0.0
    for k in 0:3, j in 0:3, i in 0:3
        ret += a[1+i+4*j+16*k]*(rx^i)*(ry^j)*(rz^k)
    end

    # release_buffers((X,a))
    return ret
end



end
