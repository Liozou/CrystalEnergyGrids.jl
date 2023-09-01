using StaticArrays

"""
    EnergyGrid

Representation of an interpolable energy grid.

Use [`interpolate_grid`](@ref) to access the value of the grid at any point in space.
"""
struct EnergyGrid
    spacing::Cdouble
    dims::NTuple{3,Cint}
    size::NTuple{3,Cdouble}
    shift::NTuple{3,Cdouble}
    Δ::NTuple{3,Cdouble}
    unitcell::NTuple{3,Cdouble}
    num_unitcell::NTuple{3,Cint}
    ewald_precision::Cfloat # Inf for VdW grid, -Inf for empty grid
    mat::SMatrix{3,3,Float64,9}
    invmat::SMatrix{3,3,Float64,9}
    higherorder::Bool # true if grid contains derivatives, false if only raw values
    grid::Array{Cfloat,4}
end
function EnergyGrid()
    EnergyGrid(0.0, (0,0,0), (0.0,0.0,0.0), (0.0,0.0,0.0), (0.0,0.0,0.0), (0.0,0.0,0.0), (0,0,0), -Inf, zero(SMatrix{3,3,Float64,9}), zero(SMatrix{3,3,Float64,9}), false, Array{Cfloat,4}(undef, 0, 0, 0, 0))
end
function Base.show(io::IO, ::MIME"text/plain", rg::EnergyGrid)
    if rg.ewald_precision == -Inf
        print("Empty grid")
        return
    end
    print(io, rg.ewald_precision == Inf ? "VdW" : "Coulomb", " grid with ")
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

Parse a .grid file and return the corresponding `EnergyGrid`.

Set `iscoulomb` if the file corresponds to an eletrostatics grid. `mat` is the unit cell
matrix of the framework.
"""
function parse_grid(file, iscoulomb, mat=nothing)
    open(file) do io
        spacing = read(io, Cdouble)
        dims = @triplet read(io, Cint)
        size = @triplet read(io, Cdouble)
        shift = @triplet read(io, Cdouble)
        Δ = @triplet read(io, Cdouble)
        unitcell = @triplet read(io, Cdouble)
        num_unitcell = @triplet read(io, Cint)

        # @printf "ShiftGrid: %g %g %g\n" shift[1] shift[2] shift[3]
        # @printf "SizeGrid: %g %g %g\n" size[1] size[2] size[3]
        # @printf "Number of grid points: %d %d %d\n" dims[1] dims[2] dims[3]

        ewald_precision = iscoulomb ? read(io, Cdouble) : Inf
        grid = Array{Cfloat, 4}(undef, dims[3]+1, dims[2]+1, dims[1]+1, 8)
        read!(io, grid)
        grid .*= NoUnits(ENERGY_TO_KELVIN/u"K")
        newmat = if mat isa AbstractMatrix
            mat
        else
            pos = position(io)
            seekend(io)
            if position(io) == pos
                error("Missing `mat` argument to `parse_grid` not provided by the grid.")
            end
            seek(io, pos)
            read(io, SMatrix{3,3,Float64,9})
        end
        EnergyGrid(spacing, dims, size, shift, Δ, unitcell, num_unitcell, ewald_precision, newmat, inv(newmat), true, grid)
    end
end

function _create_grid_common(file, framework::AbstractSystem{3}, spacing::Float64, cutoff::Float64)
    # isfile(file) && error(lazy"File $file already exists, please remove it if you want to overwrite it.")

    a, b, c = [NoUnits.(x./u"Å") for x in bounding_box(framework)]
    unitcell = norm(a), norm(b), norm(c)
    mat = stack3((a,b,c))
    invmat = inv(mat)
    size = SVector{3,Cdouble}(abs.(a)) + SVector{3,Cdouble}(abs.(b)) + SVector{3,Cdouble}(abs.(c))
    shift = SVector{3,Cdouble}(min.(a, 0.0)) + SVector{3,Cdouble}(min.(b, 0.0)) + SVector{3,Cdouble}(min.(c, 0.0))
    _dims = floor.(Cint, (size ./ spacing))
    dims = _dims + iseven.(_dims)
    Δ = size ./ dims
    num_unitcell = Cint.(find_supercell((a, b, c), cutoff))

    # @printf "ShiftGrid: %g %g %g\n" shift[1] shift[2] shift[3]
    # @printf "SizeGrid: %g %g %g\n" size[1] size[2] size[3]
    # @printf "Number of grid points: %d %d %d\n" dims[1] dims[2] dims[3]

    open(file, "w") do f
        write(f, Cdouble(spacing))
        write(f, dims[1], dims[2], dims[3])
        write(f, size[1], size[2], size[3])
        write(f, shift[1], shift[2], shift[3])
        write(f, Δ[1], Δ[2], Δ[3])
        write(f, unitcell[1], unitcell[2], unitcell[3])
        write(f, num_unitcell[1], num_unitcell[2], num_unitcell[3])
    end

    return size, shift, dims, Δ, num_unitcell, mat, invmat
end

function _set_gridpoint!(grid, i, j, k, Δ, λ, derivatives)
    value, ∂1, ∂2, ∂3 = derivatives
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

function create_grid_vdw(file, framework::AbstractSystem{3}, forcefield::ForceField, spacing::Float64, atom::Symbol)
    size, shift, dims, Δ, num_unitcell, mat, invmat = _create_grid_common(file, framework, spacing, NoUnits(forcefield.cutoff/u"Å"))
    grid = Array{Cfloat,4}(undef, dims[3]+1, dims[2]+1, dims[1]+1, 8)
    probe_vdw = ProbeSystem(framework, forcefield, atom)
    λ = inv(NoUnits(ENERGY_TO_KELVIN/u"K"))
    @threads for i in 0:dims[1]
        for j in 0:dims[2], k in 0:dims[3]
            pos = SVector{3,Float64}((i*size[1]/dims[1] + shift[1], j*size[2]/dims[2] + shift[2], k*size[3]/dims[3] + shift[3]))
            derivatives = compute_derivatives_vdw(probe_vdw, pos)
            _set_gridpoint!(grid, i, j, k, Δ, λ, derivatives)
        end
    end
    open(file, "a") do f
        write(f, grid)
        write(f, mat) # not part of RASPA grids
    end
    grid
end

function create_grid_coulomb(file, framework::AbstractSystem{3}, forcefield::ForceField, spacing::Float64, _ewald=nothing)
    size, shift, dims, Δ, num_unitcell, mat, invmat = _create_grid_common(file, framework, spacing, 12.0)
    ewald = if _ewald isa EwaldFramework
        _ewald
    else
        initialize_ewald(framework, num_unitcell)
    end
    grid = Array{Cfloat,4}(undef, dims[3]+1, dims[2]+1, dims[1]+1, 8)
    probe_coulomb = ProbeSystem(framework, forcefield)
    @threads for i in 0:dims[1]
        for j in 0:dims[2], k in 0:dims[3]
            pos = SVector{3,Float64}((i*size[1]/dims[1] + shift[1], j*size[2]/dims[2] + shift[2], k*size[3]/dims[3] + shift[3]))
            derivatives = compute_derivatives_ewald(probe_coulomb, ewald, pos)
            _set_gridpoint!(grid, i, j, k, Δ, COULOMBIC_CONVERSION_FACTOR, derivatives)
        end
    end
    open(file, "a") do f
        write(f, ewald.precision)
        write(f, grid)
        write(f, mat) # not part of RASPA grids
    end
    grid
end

"""
    interpolate_grid(g::EnergyGrid, point)

Interpolate grid `g` at the given `point`, which should be a triplet of coordinates (with
their corresponding unit).
"""
function interpolate_grid(g::EnergyGrid, point)
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
            return 1e100*u"K"
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
    coulomb::EnergyGrid
    grids::Vector{EnergyGrid}
    atomsidx::Vector{Int} # index of the grid corresponding to the atom
    ewald::EwaldFramework
    forcefield::ForceField
    block::Union{Nothing, BitArray{3}}
end

"""
    energy_point(setup::CrystalEnergySetup, positions)

Compute the energy of a guest molecule whose atoms are at the given `positions` in the
framework. `positions` is a list of triplet of coordinates (with their unit).

Return a pair `(vdw, coulomb)` where `vdw` is the Van der Waals contribution to the energy
and `coulomb` is the electrostatic ones, both in K.

!!! warning
    No validation is used to ensure that the input `positions` are consistent with the
    shape of the molecule.
"""
function energy_point(setup::CrystalEnergySetup, positions)
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
    gr0 = setup.coulomb.ewald_precision == -Inf ? setup.grids[1] : setup.coulomb
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
                allvals[:,iA,iB,iC] .= 1e100
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
