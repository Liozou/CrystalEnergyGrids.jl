using StaticArrays

"""
    EnergyGrid

Representation of an interpolable energy grid.

Use [`interpolate_grid`](@ref) to access the value of the grid at any point in space.
"""
struct EnergyGrid
    csetup::GridCoordinatesSetup
    num_unitcell::NTuple{3,Cint}
    ewald_precision::Cfloat # Inf for VdW grid, -Inf for empty grid
    higherorder::Bool # true if grid contains derivatives, false if only raw values
    grid::Array{Cfloat,4}
end
function EnergyGrid()
    EnergyGrid(GridCoordinatesSetup(), (0, 0, 0), -Inf, false, Array{Cfloat,4}(undef, 0, 0, 0, 0))
end
function Base.show(io::IO, ::MIME"text/plain", rg::EnergyGrid)
    if rg.ewald_precision == -Inf
        print("Empty grid")
        return
    end
    print(io, rg.ewald_precision == Inf ? "VdW" : "Coulomb", " grid with ")
    join(io, rg.csetup.dims .+ 1, '×')
    print(io, " points for a ")
    if rg.num_unitcell != (1,1,1)
        join(io, rg.num_unitcell, '×')
        print(io, " supercell of the ")
    end
    print("unit cell of size ")
    join(io, rg.csetup.unitcell, " × ")
end

macro triplet(foo)
    esc(quote
        ($foo, $foo, $foo)
    end)
end

"""
    parse_grid(file, iscoulomb, mat)

Parse a .grid file and return the corresponding `EnergyGrid`.

Set `iscoulomb` if the file corresponds to an eletrostatics grid. `mat` is the unit cell
matrix of the framework.
"""
function parse_grid(file, iscoulomb, mat=nothing)
    open(file) do io
        spacing = read(io, Cdouble)*u"Å"
        dims = @triplet read(io, Cint)
        size = @triplet read(io, Cdouble).*u"Å"
        shift = @triplet read(io, Cdouble).*u"Å"
        Δ = @triplet read(io, Cdouble).*u"Å"
        unitcell = @triplet read(io, Cdouble).*u"Å"
        num_unitcell = @triplet read(io, Cint)

        # @printf "ShiftGrid: %g %g %g\n" shift[1] shift[2] shift[3]
        # @printf "SizeGrid: %g %g %g\n" size[1] size[2] size[3]
        # @printf "Number of grid points: %d %d %d\n" dims[1] dims[2] dims[3]

        ewald_precision = iscoulomb ? read(io, Cdouble) : Inf
        grid = Array{Cfloat, 4}(undef, dims[3]+1, dims[2]+1, dims[1]+1, 8)
        read!(io, grid)
        grid .*= NoUnits(ENERGY_TO_KELVIN/u"K")
        newmat = if mat isa CellMatrix
            mat
        elseif mat isa AbstractMatrix
            CellMatrix(mat)
        else
            pos = position(io)
            seekend(io)
            if position(io) == pos
                error("Missing `mat` argument to `parse_grid` not provided by the grid.")
            end
            seek(io, pos)
            CellMatrix(read(io, SMatrix{3,3,Float64,9}))
        end
        EnergyGrid(GridCoordinatesSetup(newmat, spacing, dims, size, shift, unitcell, Δ), num_unitcell, ewald_precision, true, grid)
    end
end

macro write_no_unit(x)
    return esc(quote
        write(f, NoUnits(($x)[1]/u"Å"), NoUnits(($x)[2]/u"Å"), NoUnits(($x)[3]/u"Å"))
    end)
end

function _create_grid_common(file, framework::AbstractSystem{3}, spacing::TÅ, cutoff::TÅ)
    # isfile(file) && error(lazy"File $file already exists, please remove it if you want to overwrite it.")

    csetup = GridCoordinatesSetup(framework, spacing)
    num_unitcell = Cint.(find_supercell(framework, cutoff))

    # @printf "ShiftGrid: %g %g %g\n" shift[1] shift[2] shift[3]
    # @printf "SizeGrid: %g %g %g\n" size[1] size[2] size[3]
    # @printf "Number of grid points: %d %d %d\n" dims[1] dims[2] dims[3]

    open(file, "w") do f
        write(f, Cdouble(NoUnits(spacing/u"Å")))
        write(f, csetup.dims[1], csetup.dims[2], csetup.dims[3])
        @write_no_unit csetup.size
        @write_no_unit csetup.shift
        @write_no_unit csetup.Δ
        @write_no_unit csetup.unitcell
        write(f, num_unitcell[1], num_unitcell[2], num_unitcell[3])
    end

    return csetup, num_unitcell
end

function _set_gridpoint!(grid, i, j, k, Δ, λ, λ⁻¹e7, derivatives)
    value, ∂1, ∂2, ∂3 = derivatives
    if value > λ⁻¹e7
        value = 2*λ⁻¹e7
        ∂1 = clamp.(∂1, -λ⁻¹e7, λ⁻¹e7)
        ∂2 = zero(SVector{3,Float64})
        ∂3 = 0.0
    end
    grid[k+1,j+1,i+1,1] = value * λ
    grid[k+1,j+1,i+1,2] = ∂1[1]*Δ[1] * λ
    grid[k+1,j+1,i+1,3] = ∂1[2]*Δ[2] * λ
    grid[k+1,j+1,i+1,4] = ∂1[3]*Δ[3] * λ
    grid[k+1,j+1,i+1,5] = ∂2[1]*(Δ[1]*Δ[2]) * λ
    grid[k+1,j+1,i+1,6] = ∂2[2]*(Δ[1]*Δ[3]) * λ
    grid[k+1,j+1,i+1,7] = ∂2[3]*(Δ[2]*Δ[3]) * λ
    grid[k+1,j+1,i+1,8] = ∂3*(Δ[1]*Δ[2]*Δ[3]) * λ
    nothing
end

function create_grid_vdw(file, framework::AbstractSystem{3}, forcefield::ForceField, spacing::TÅ, atom::Symbol)
    cset, num_unitcell = _create_grid_common(file, framework, spacing, forcefield.cutoff)
    grid = Array{Cfloat,4}(undef, cset.dims[3]+1, cset.dims[2]+1, cset.dims[1]+1, 8)
    probe_vdw = ProbeSystem(framework, forcefield, atom)
    Δ = NoUnits.(cset.Δ/u"Å")
    λ⁻¹ = NoUnits(ENERGY_TO_KELVIN/u"K")
    λ = inv(λ⁻¹)
    @threads for i in 0:cset.dims[1]
        for j in 0:cset.dims[2], k in 0:cset.dims[3]
            pos = abc_to_xyz(cset, i, j, k)
            derivatives = compute_derivatives_vdw(probe_vdw, pos)
            _set_gridpoint!(grid, i, j, k, Δ, λ, λ⁻¹*1e7, derivatives)
        end
    end
    open(file, "a") do f
        write(f, grid)
        write(f, NoUnits.(cset.cell.mat./u"Å")) # not part of RASPA grids
    end
    grid
end

function create_grid_coulomb(file, framework::AbstractSystem{3}, forcefield::ForceField, spacing::TÅ, _ewald=nothing)
    cset, num_unitcell = _create_grid_common(file, framework, spacing, 12.0u"Å")
    ewald = if _ewald isa EwaldFramework
        _ewald
    else
        initialize_ewald(framework, num_unitcell)
    end
    grid = Array{Cfloat,4}(undef, cset.dims[3]+1, cset.dims[2]+1, cset.dims[1]+1, 8)
    probe_coulomb = ProbeSystem(framework, forcefield)
    Δ = NoUnits.(cset.Δ/u"Å")
    λ = COULOMBIC_CONVERSION_FACTOR
    λ⁻¹e7 = inv(λ)*1e7
    @threads for i in 0:cset.dims[1]
        for j in 0:cset.dims[2], k in 0:cset.dims[3]
            pos = abc_to_xyz(cset, i, j, k)
            derivatives = compute_derivatives_ewald(probe_coulomb, ewald, pos)
            _set_gridpoint!(grid, i, j, k, Δ, λ, λ⁻¹e7, derivatives)
        end
    end
    open(file, "a") do f
        write(f, ewald.precision)
        write(f, grid)
        write(f, NoUnits.(cset.cell.mat./u"Å")) # not part of RASPA grids
    end
    grid
end


function BlockFile(g::EnergyGrid)
    a, b, c = g.csetup.dims .+ 1
    block = falses(a, b, c)
    for i in 1:a-1, j in 1:b-1, k in 1:c-1
        if g.grid[k,j,i,1] > 5e6
            block[i  ,j  ,k  ] = true
            block[i+1,j  ,k  ] = true
            block[i  ,j+1,k  ] = true
            block[i+1,j+1,k  ] = true
            block[i  ,j  ,k+1] = true
            block[i+1,j  ,k+1] = true
            block[i  ,j+1,k+1] = true
            block[i+1,j+1,k+1] = true
        end
    end
    BlockFile(g.csetup, block)
end

"""
    interpolate_grid(g::EnergyGrid, point)

Interpolate grid `g` at the given `point`, which should be a triplet of coordinates (with
their corresponding unit).
"""
function interpolate_grid(g::EnergyGrid, point)
    shifted = offsetpoint(point, g.csetup)
    p0 = floor.(Int, shifted)
    # The following is p1 = p0 .+ 1 adapted to avoid BoundsError on numerical imprecision
    p1 = p0 .+ (p0[1] != size(g.grid, 3), p0[2] != size(g.grid, 2), p0[3] != size(g.grid, 1))
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
        if g.ewald_precision == Inf && any(>(5e6), @view X[1:8])
            # release_buffers((X,a))
            return 1e100*u"K"
        end
        a = COEFF * X
        rx2 = rx*rx; rxs = SVector{4,Float64}((1.0, rx, rx2, rx2*rx))
        ry2 = ry*ry; rys = SVector{4,Float64}((1.0, ry, ry2, ry2*ry))
        rz2 = rz*rz; rzs = SVector{4,Float64}((1.0, rz, rz2, rz2*rz))
        @inbounds for k in 0:3, j in 0:3, i in 0:3
            ret += a[1+i+4*j+16*k]*rxs[1+i]*rys[1+j]*rzs[1+k]
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
    coulomb::EnergyGrid
    grids::Vector{EnergyGrid}
    atomsidx::Vector{Int} # index of the grid corresponding to the atom
    ewald::EwaldFramework
    forcefield::ForceField
    block::BlockFile
end

"""
    energy_point(setup::CrystalEnergySetup, positions)

Compute the energy of a guest molecule whose atoms are at the given `positions` in the
framework. `positions` is a list of triplet of coordinates (with their unit).

Return a pair `(vdw, coulomb)` where `vdw` is the Van der Waals contribution to the energy
and `coulomb` is the electrostatic ones, both in K.

If one of the atom is on a blocking sphere, return `(1e100u"K", 0.0u"K")`

!!! warning
    No validation is used to ensure that the input `positions` are consistent with the
    shape of the molecule.
"""
function energy_point(setup::CrystalEnergySetup, positions)
    for pos in positions
        setup.block[pos] && return (1e100u"K", 0.0u"K")
    end
    num_atoms = length(setup.atomsidx)
    vdw = sum(interpolate_grid(setup.grids[setup.atomsidx[i]], positions[i]) for i in 1:num_atoms)
    setup.coulomb.ewald_precision == -Inf && return vdw, 0.0u"K"
    coulomb_direct = sum((Float64(setup.molecule[i,:atomic_charge]/u"e_au"))*interpolate_grid(setup.coulomb, positions[i]) for i in 1:num_atoms)
    newmolecule = ChangePositionSystem(setup.molecule, positions)
    coulomb_reciprocal = compute_ewald(setup.ewald, (newmolecule,))
    return (vdw, coulomb_direct + coulomb_reciprocal)
end


"""
    energy_grid(setup::CrystalEnergySetup, step, num_rotate=40)

Compute the energy on the the system given by `setup` on a regular grid with the given
`step` (with its unit).

If the guest molecule is not monoatomic, the first axe of the returned grid will represent
the rotation angle of the molecule, and its size will be greater or equal to `num_rotate`.
Otherwise (or if `num_rotate == 0`), the first axe can be dropped through`dropdims`.

A value of 1e100 in the grid indicates an inaccessible point due to blocking spheres.
"""
function energy_grid(setup::CrystalEnergySetup, step, num_rotate=40)
    axeA, axeB, axeC = bounding_box(setup.framework)
    numA = floor(Int, norm(axeA) / step) + 1
    numB = floor(Int, norm(axeB) / step) + 1
    numC = floor(Int, norm(axeC) / step) + 1
    stepA = (axeA / norm(axeA)) * step
    stepB = (axeB / norm(axeB)) * step
    stepC = (axeC / norm(axeC)) * step
    __pos = position(setup.molecule) / u"Å"
    rotpos::Vector{Vector{SVector{3,Float64}}} = if num_rotate == 0 || length(setup.molecule) == 1
        [[SVector{3}(p) for p in __pos]]
    else
        rots, _ = get_rotation_matrices(setup.molecule, num_rotate)
        [[SVector{3}(r*p) for p in __pos] for r in rots]
    end
    pre_printwarn = PRINT_CHARGE_WARNING[]
    # do one empty computation to trigger the warnings if need be
    if pre_printwarn
        energy_point(setup, [p*u"Å" for p in first(rotpos)])
        PRINT_CHARGE_WARNING[] = false
    end
    allvals = Array{Float64}(undef, length(rotpos), numA, numB, numC)
    allvals .= NaN
    Base.Threads.@threads for iABC in CartesianIndices((numA, numB, numC))
        iA, iB, iC = Tuple(iABC)
        thisofs = (iA-1)*stepA + (iB-1)*stepB + (iC-1)*stepC
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
    pre_printwarn && (PRINT_CHARGE_WARNING[] = pre_printwarn)
    allvals
end
