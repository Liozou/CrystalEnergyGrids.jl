module CrystalEnergyGrids

using StaticArrays
using LinearAlgebra: LinearAlgebra, norm, det, cross, dot
using Pkg: TOML

using StaticArrays
using OffsetArrays
using AtomsBase
using Scratch: @get_scratch!

export CrystalEnergyGrid, parse_grid, interpolate_grid
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
include("energy.jl")
include("raspa.jl")
include("ewald.jl")
include("lammpstrj_parser.jl")
include("hypernettedchain.jl")
include("rdf.jl")
include("mdft.jl")

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


struct CrystalEnergyGrid
    spacing::Cdouble
    dims::NTuple{3,Cint}
    size::NTuple{3,Cdouble}
    shift::NTuple{3,Cdouble}
    Δ::NTuple{3,Cdouble}
    unitcell::NTuple{3,Cdouble}
    num_unitcell::NTuple{3,Cint}
    ε_Ewald::Cfloat
    mat::SMatrix{3,3,Float64,9}
    invmat::SMatrix{3,3,Float64,9}
    grid::Array{Cfloat,4}
end
function CrystalEnergyGrid()
    CrystalEnergyGrid(0.0, (0,0,0), (0.0,0.0,0.0), (0.0,0.0,0.0), (0.0,0.0,0.0), (0.0,0.0,0.0), (0,0,0), -Inf, zero(SMatrix{3,3,Float64,9}), zero(SMatrix{3,3,Float64,9}), Array{Cfloat,4}(undef, 0, 0, 0, 0))
end
function Base.show(io::IO, ::MIME"text/plain", rg::CrystalEnergyGrid)
    if rg.ε_Ewald == -Inf
        print("Empty grid")
        return
    end
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

function parse_grid(file, iscoulomb, mat)
    invmat = inv(mat)
    open(file) do io
        spacing = read(io, Cdouble)
        dims = @triplet read(io, Cint)
        size = @triplet read(io, Cdouble)
        shift = @triplet read(io, Cdouble)
        Δ = @triplet read(io, Cdouble)
        unitcell = @triplet read(io, Cdouble)
        num_unitcell = @triplet read(io, Cint)
        ε_Ewald = iscoulomb ? read(io, Cdouble) : Inf
        grid = Array{Cfloat, 4}(undef, dims[3]+1, dims[2]+1, dims[1]+1, 8)
        read!(io, grid)
        grid .*= NoUnits(ENERGY_TO_KELVIN/u"K")
        CrystalEnergyGrid(spacing, dims, size, shift, Δ, unitcell, num_unitcell, ε_Ewald, mat, invmat, grid)
    end
end

function interpolate_grid(g::CrystalEnergyGrid, point)
    shifted = offsetpoint(point, g.mat, g.invmat, g.shift, g.size, g.dims)
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
    return ret*u"K"
end


struct CrystalEnergySetup{TFramework,TMolecule}
    framework::TFramework
    molecule::TMolecule
    coulomb::CrystalEnergyGrid
    grids::Vector{CrystalEnergyGrid}
    atomsidx::Vector{Int} # index of the grid corresponding to the atom
    ewald::Tuple{Float64,Int,Int,Int,SMatrix{3,3,Float64,9},Float64,Int,Float64,Float64,Vector{Float64},Float64,Vector{ComplexF64},Float64}
    block::Union{Nothing, BitArray{3}}
end

function energy_point(setup::CrystalEnergySetup, positions)
    num_atoms = length(setup.atomsidx)
    pos_strip::Vector{SVector{3,Float64}} = positions / u"Å"
    vdw = sum(interpolate_grid(setup.grids[setup.atomsidx[i]], pos_strip[i]) for i in 1:num_atoms)
    coulomb_direct = if setup.coulomb.ε_Ewald == -Inf
        0.0u"K"
    else
        sum((NoUnits(setup.molecule[i,:atomic_charge]/u"e_au"))*interpolate_grid(setup.coulomb, pos_strip[i]) for i in 1:num_atoms)
    end
    newmolecule = ChangePositionSystem(setup.molecule, positions)
    host_adsorbate_reciprocal, adsorbate_adsorbate_reciprocal = compute_ewald(setup.ewald, (newmolecule,))
    # @show vdw, coulomb_direct, host_adsorbate_reciprocal, adsorbate_adsorbate_reciprocal
    return (vdw, coulomb_direct + host_adsorbate_reciprocal + adsorbate_adsorbate_reciprocal)
end


function energy_grid(setup::CrystalEnergySetup, step, num_rotate=40)
    axeA, axeB, axeC = bounding_box(setup.framework)
    numA = floor(Int, norm(axeA) / step) + 1
    numB = floor(Int, norm(axeB) / step) + 1
    numC = floor(Int, norm(axeC) / step) + 1
    # grid = Array{Float64}(undef, numA, numB, numC)
    # mingrid = Array{Float64}(undef, numA, numB, numC)
    stepA = (axeA / norm(axeA)) * step
    stepB = (axeB / norm(axeB)) * step
    stepC = (axeC / norm(axeC)) * step
    gr0 = setup.coulomb.ε_Ewald == -Inf ? setup.grids[1] : setup.coulomb
    mat = gr0.mat
    invmat = gr0.invmat
    __pos = position(setup.molecule) / u"Å"
    rotpos::Vector{Vector{SVector{3,Float64}}} = if num_rotate == 0 || length(setup.molecule) == 1
        [[SVector{3}(p) for p in __pos]]
    # elseif num_rotate == -1
        # [__pos, __switch_yz.(__pos), __switch_xz.(__pos)]
    else
        rots, _ = get_rotation_matrices(setup.molecule, num_rotate)
        [[SVector{3}(r*p) for p in __pos] for r in rots]
    end
    allvals = Array{Float64}(undef, length(rotpos), numA, numB, numC)
    allvals .= NaN
    Base.Threads.@threads for idx in CartesianIndices((numA, numB, numC))
        iA, iB, iC = Tuple(idx)
        thisofs = (iA-1)*stepA + (iB-1)*stepB + (iC-1)*stepC
        if setup.block isa BitArray && num_rotate >= 0
            a0, b0, c0 = floor.(Int, offsetpoint(NoUnits.(thisofs/u"Å"), mat, invmat, gr0.shift, gr0.size, gr0.dims))
            a1 = a0 + 1; b1 = b0 + 1; c1 = c0 + 1;
            if setup.block[a0,b0,c0]+setup.block[a1,b0,c0]+setup.block[a0,b1,c0]+setup.block[a1,b1,c0]+setup.block[a0,b0,c1]+setup.block[a1,b0,c1]+setup.block[a0,b1,c1]+setup.block[a1,b1,c1] > 3
                # grid[iA,iB,iC] = Inf
                allvals[:,iA,iB,iC] .= 1e8
                continue
            end
        end
        # vals = 0.0
        # minval = Inf
        for (k, pos) in enumerate(rotpos)
            ofs = if num_rotate < 0
                rand()*numA*stepA + rand()*numB*stepB + rand()*numC*stepC
            else
                thisofs
            end
            newval = NoUnits(sum(energy_point(setup, [SVector{3}(ofs + p*u"Å") for p in pos]))/u"K")
            # vals += newval
            # vals == Inf && break
            # minval = min(minval, newval)
            allvals[k,iA,iB,iC] = newval
            # if newval < -10000
            #     @show iA, iB, iC
            #     @show ofs
            #     @show SVector{3}(wrap_atom(NoUnits.(ofs/u"Å") + p, mat, invmat)*u"Å" for p in pos)
            #     throw("!")
            # end
        end
        # grid[iA,iB,iC] = vals / length(rotpos)
        # mingrid[iA,iB,iC] = minval
    end
    allvals
end


end
