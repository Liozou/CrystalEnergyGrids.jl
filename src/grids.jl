using StaticArrays

"""
    CrystalEnergyGrid

Representation of an interpolable energy grid.

Use [`interpolate_grid`](@ref) to access the value of the grid at any point in space.
"""
struct CrystalEnergyGrid
    spacing::Cdouble
    dims::NTuple{3,Cint}
    size::NTuple{3,Cdouble}
    shift::NTuple{3,Cdouble}
    Δ::NTuple{3,Cdouble}
    unitcell::NTuple{3,Cdouble}
    num_unitcell::NTuple{3,Cint}
    ε_Ewald::Cfloat # Inf for VdW grid, -Inf for empty grid
    mat::SMatrix{3,3,Float64,9}
    invmat::SMatrix{3,3,Float64,9}
    higherorder::Bool # true if grid contains derivatives, false if only raw values
    grid::Array{Cfloat,4}
end
function CrystalEnergyGrid()
    CrystalEnergyGrid(0.0, (0,0,0), (0.0,0.0,0.0), (0.0,0.0,0.0), (0.0,0.0,0.0), (0.0,0.0,0.0), (0,0,0), -Inf, zero(SMatrix{3,3,Float64,9}), zero(SMatrix{3,3,Float64,9}), false, Array{Cfloat,4}(undef, 0, 0, 0, 0))
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

"""
    parse_grid(file, iscoulomb, mat)

Parse a .grid file and return the corresponding `CrystalEnergyGrid`.

Set `iscoulomb` if the file corresponds to an eletrostatics grid. `mat` is the unit cell
matrix of the framework.
"""
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
        CrystalEnergyGrid(spacing, dims, size, shift, Δ, unitcell, num_unitcell, ε_Ewald, mat, invmat, true, grid)
    end
end

"""
    interpolate_grid(g::CrystalEnergyGrid, point)

Interpolate grid `g` at the given `point`, which should be a triplet of coordinates (with
their corresponding unit).
"""
function interpolate_grid(g::CrystalEnergyGrid, point)
    shifted = offsetpoint(point, g.mat, g.invmat, g.shift, g.size, g.dims)
    p0 = floor.(Int, shifted)
    p1 = p0 .+ 1
    r = shifted .- p0
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    rx, ry, rz = r
    # X, a = acquire_buffers()

    ret = 0.0
    if g.higherorder # use derivatives
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
            return 1e8*u"K"
        end
        a = COEFF * X
        for k in 0:3, j in 0:3, i in 0:3
            ret += a[1+i+4*j+16*k]*(rx^i)*(ry^j)*(rz^k)
        end
    else # no derivatives
        mrx = 1 - rx
        mry = 1 - ry
        mrz = 1 - rz
        ret = begin
            g.grid[x0,y0,z0,1]*mrx*mry*mrz + g.grid[x1,y0,z0,1]*rx*mry*mrz +
            g.grid[x0,y1,z0,1]*mrx*ry*mrz  + g.grid[x0,y0,z1,1]*mrx*mry*rz +
            g.grid[x1,y1,z0,1]*rx*ry*mrz   + g.grid[x1,y0,z1,1]*rx*mry*rz  +
            g.grid[x0,y1,z1,1]*mrx*ry*rz   + g.grid[x1,y1,z1,1]*rx*ry*rz
        end
    end

    # release_buffers((X,a))
    return ret*u"K"
end

"""
    CrystalEnergySetup{TFramework,TMolecule}

Structure representing a system composed of a framework material (of type `TFramework`) and
a guest molecule (of type `TMolecule`).
Both `TFramework` and `TMolecule` should obey the `AbstractSystem` interface of `AtomsBase`.

See also [`setup_RASPA`](@ref) to build a `CrystalEnergySetup` from RASPA2 files.
"""
struct CrystalEnergySetup{TFramework,TMolecule}
    framework::TFramework
    molecule::TMolecule
    coulomb::CrystalEnergyGrid
    grids::Vector{CrystalEnergyGrid}
    atomsidx::Vector{Int} # index of the grid corresponding to the atom
    ewald::EwaldContext
    block::Union{Nothing, BitArray{3}}
end

"""
    energy_point(setup::CrystalEnergySetup, positions)

Compute the energy of a guest molecule whose atoms are at the given `positions` in the
framework. `positions` is a list of triplet of coordinates (with their unit).

Return a pair `(vdw, coulomb)` where `vdw` is the Van der Waals contribution to the energy
and `coulomb` is the electrostatic ones, both in K.

!!! warn
    No validation is used to ensure that the input `positions` are consistent with the
    shape of the molecule.
"""
function energy_point(setup::CrystalEnergySetup, positions)
    num_atoms = length(setup.atomsidx)
    vdw = sum(interpolate_grid(setup.grids[setup.atomsidx[i]], positions[i]) for i in 1:num_atoms)
    coulomb_direct = if setup.coulomb.ε_Ewald == -Inf
        0.0u"K"
    else
        sum((NoUnits(setup.molecule[i,:atomic_charge]/u"e_au"))*interpolate_grid(setup.coulomb, positions[i]) for i in 1:num_atoms)
    end
    newmolecule = ChangePositionSystem(setup.molecule, positions)
    host_adsorbate_reciprocal, adsorbate_adsorbate_reciprocal = compute_ewald(setup.ewald, (newmolecule,))
    # @show vdw, coulomb_direct, host_adsorbate_reciprocal, adsorbate_adsorbate_reciprocal
    return (vdw, coulomb_direct + host_adsorbate_reciprocal + adsorbate_adsorbate_reciprocal)
end


"""
    energy_grid(setup::CrystalEnergySetup, step, num_rotate=40)

Compute the energy on the the system given by `setup` on a regular grid with the given
`step` (with its unit).

If the guest molecule is not monoatomic, the first axe of the returned grid will represent
the rotation angle of the molecule, and its size will be greater or equal to `num_rotate`.
Otherwise (or if `num_rotate == 0`), the first axe can be dropped through`dropdims`.

A value of 1e8 in the grid indicates an inaccessible point due to blocking spheres.
"""
function energy_grid(setup::CrystalEnergySetup, step, num_rotate=40)
    axeA, axeB, axeC = bounding_box(setup.framework)
    numA = floor(Int, norm(axeA) / step) + 1
    numB = floor(Int, norm(axeB) / step) + 1
    numC = floor(Int, norm(axeC) / step) + 1
    stepA = (axeA / norm(axeA)) * step
    stepB = (axeB / norm(axeB)) * step
    stepC = (axeC / norm(axeC)) * step
    gr0 = setup.coulomb.ε_Ewald == -Inf ? setup.grids[1] : setup.coulomb
    mat = gr0.mat
    invmat = gr0.invmat
    __pos = position(setup.molecule) / u"Å"
    rotpos::Vector{Vector{SVector{3,Float64}}} = if num_rotate == 0 || length(setup.molecule) == 1
        [[SVector{3}(p) for p in __pos]]
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
            a0, b0, c0 = floor.(Int, offsetpoint(thisofs, mat, invmat, gr0.shift, gr0.size, gr0.dims))
            a1 = a0 + 1; b1 = b0 + 1; c1 = c0 + 1;
            if setup.block[a0,b0,c0]+setup.block[a1,b0,c0]+setup.block[a0,b1,c0]+setup.block[a1,b1,c0]+setup.block[a0,b0,c1]+setup.block[a1,b0,c1]+setup.block[a0,b1,c1]+setup.block[a1,b1,c1] > 3
                allvals[:,iA,iB,iC] .= 1e8
                continue
            end
        end
        for (k, pos) in enumerate(rotpos)
            ofs = if num_rotate < 0
                rand()*numA*stepA + rand()*numB*stepB + rand()*numC*stepC
            else
                thisofs
            end
            newval = NoUnits(sum(energy_point(setup, [SVector{3}(ofs + p*u"Å") for p in pos]))/u"K")
            allvals[k,iA,iB,iC] = newval
        end
    end
    allvals
end
