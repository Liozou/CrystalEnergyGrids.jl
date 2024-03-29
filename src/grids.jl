using StaticArrays

"""
    EnergyGrid

Representation of an interpolable energy grid.

Use [`interpolate_grid`](@ref) to access the value of the grid at any point in space.
"""
struct EnergyGrid
    csetup::GridCoordinatesSetup
    num_unitcell::NTuple{3,Cint}
    ewald_precision::Cfloat # Inf for VdW grid, -Inf for zero grid, NaN for invalid grid
    higherorder::Bool # true if grid contains derivatives, false if only raw values
    grid::Array{Cfloat,4}
end
function EnergyGrid(zero::Bool)
    EnergyGrid(GridCoordinatesSetup(), (0, 0, 0), ifelse(zero, -Inf, NaN), false, Array{Cfloat,4}(undef, 0, 0, 0, 0))
end
function Base.show(io::IO, ::MIME"text/plain", rg::EnergyGrid)
    if rg.ewald_precision == -Inf
        print("Zero grid")
        return
    elseif isnan(rg.ewald_precision)
        print("No grid")
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
        grid .*= GRID_TO_KELVIN
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
        write(io, NoUnits(($x)[1]/u"Å"), NoUnits(($x)[2]/u"Å"), NoUnits(($x)[3]/u"Å"))
    end)
end

function _setup_grid_common(framework::AbstractSystem{3}, spacing::TÅ, cutoff::TÅ)
    csetup = GridCoordinatesSetup(framework, spacing)
    num_unitcell = Cint.(find_supercell(framework, cutoff))
    csetup, num_unitcell
end

function _create_grid_common(io::IO, csetup::GridCoordinatesSetup, num_unitcell)
    write(io, Cdouble(ustrip(u"Å", csetup.spacing)))
    write(io, csetup.dims[1], csetup.dims[2], csetup.dims[3])
    @write_no_unit csetup.size
    @write_no_unit csetup.shift
    @write_no_unit csetup.Δ
    @write_no_unit csetup.unitcell
    write(io, Cint(num_unitcell[1]), Cint(num_unitcell[2]), Cint(num_unitcell[3]))
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
    cset, num_unitcell = _setup_grid_common(framework, spacing, forcefield.cutoff)
    grid = Array{Cfloat,4}(undef, cset.dims[3]+1, cset.dims[2]+1, cset.dims[1]+1, 8)
    probe_vdw = ProbeSystem(framework, forcefield, atom)
    Δ = NoUnits.(cset.Δ/u"Å")
    λ⁻¹ = GRID_TO_KELVIN
    λ = inv(λ⁻¹)
    @threads for i in 0:cset.dims[1]
        for j in 0:cset.dims[2], k in 0:cset.dims[3]
            pos = abc_to_xyz(cset, i, j, k)
            derivatives = compute_derivatives_vdw(probe_vdw, pos)
            _set_gridpoint!(grid, i, j, k, Δ, λ, λ⁻¹*1e7, derivatives)
        end
    end
    open(file, "w") do f
        _create_grid_common(f, cset, num_unitcell)
        write(f, grid)
        write(f, NoUnits.(cset.cell.mat./u"Å")) # not part of RASPA grids
    end
    grid
end

function create_grid_coulomb(file, framework::AbstractSystem{3}, forcefield::ForceField, spacing::TÅ, _ewald=nothing)
    cset, num_unitcell = _setup_grid_common(framework, spacing, 12.0u"Å")
    ewald = if _ewald isa EwaldFramework
        _ewald
    else
        initialize_ewald(framework, num_unitcell)
    end
    grid = Array{Cfloat,4}(undef, cset.dims[3]+1, cset.dims[2]+1, cset.dims[1]+1, 8)
    probe_coulomb = ProbeSystem(framework, forcefield)
    Δ = NoUnits.(cset.Δ/u"Å")
    λ = ustrip(u"K*Å/e_au^2", COULOMBIC_CONVERSION_FACTOR)/GRID_TO_KELVIN
    λ⁻¹e7 = inv(λ)*1e7
    @threads for i in 0:cset.dims[1]
        for j in 0:cset.dims[2], k in 0:cset.dims[3]
            pos = abc_to_xyz(cset, i, j, k)
            derivatives = compute_derivatives_ewald(probe_coulomb, ewald, pos)
            _set_gridpoint!(grid, i, j, k, Δ, λ, λ⁻¹e7, derivatives)
        end
    end
    open(file, "w") do f
        _create_grid_common(f, cset, num_unitcell)
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
    g.ewald_precision === -Inf32 && return 0.0u"K"
    isnan(g.ewald_precision) && error("Invalid grid cannot be interpolated!")
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
        X = SA[
            g.grid[z0,y0,x0,1], g.grid[z0,y0,x1,1], g.grid[z0,y1,x0,1], g.grid[z0,y1,x1,1],
            g.grid[z1,y0,x0,1], g.grid[z1,y0,x1,1], g.grid[z1,y1,x0,1], g.grid[z1,y1,x1,1],
            g.grid[z0,y0,x0,2], g.grid[z0,y0,x1,2], g.grid[z0,y1,x0,2], g.grid[z0,y1,x1,2],
            g.grid[z1,y0,x0,2], g.grid[z1,y0,x1,2], g.grid[z1,y1,x0,2], g.grid[z1,y1,x1,2],
            g.grid[z0,y0,x0,3], g.grid[z0,y0,x1,3], g.grid[z0,y1,x0,3], g.grid[z0,y1,x1,3],
            g.grid[z1,y0,x0,3], g.grid[z1,y0,x1,3], g.grid[z1,y1,x0,3], g.grid[z1,y1,x1,3],
            g.grid[z0,y0,x0,4], g.grid[z0,y0,x1,4], g.grid[z0,y1,x0,4], g.grid[z0,y1,x1,4],
            g.grid[z1,y0,x0,4], g.grid[z1,y0,x1,4], g.grid[z1,y1,x0,4], g.grid[z1,y1,x1,4],
            g.grid[z0,y0,x0,5], g.grid[z0,y0,x1,5], g.grid[z0,y1,x0,5], g.grid[z0,y1,x1,5],
            g.grid[z1,y0,x0,5], g.grid[z1,y0,x1,5], g.grid[z1,y1,x0,5], g.grid[z1,y1,x1,5],
            g.grid[z0,y0,x0,6], g.grid[z0,y0,x1,6], g.grid[z0,y1,x0,6], g.grid[z0,y1,x1,6],
            g.grid[z1,y0,x0,6], g.grid[z1,y0,x1,6], g.grid[z1,y1,x0,6], g.grid[z1,y1,x1,6],
            g.grid[z0,y0,x0,7], g.grid[z0,y0,x1,7], g.grid[z0,y1,x0,7], g.grid[z0,y1,x1,7],
            g.grid[z1,y0,x0,7], g.grid[z1,y0,x1,7], g.grid[z1,y1,x0,7], g.grid[z1,y1,x1,7],
            g.grid[z0,y0,x0,8], g.grid[z0,y0,x1,8], g.grid[z0,y1,x0,8], g.grid[z0,y1,x1,8],
            g.grid[z1,y0,x0,8], g.grid[z1,y0,x1,8], g.grid[z1,y1,x0,8], g.grid[z1,y1,x1,8],
        ]
        if g.ewald_precision == Inf && any(>(5e6), @view X[1:8])
            # release_buffers((X,a))
            return 1e100*u"K"
        end
        a::Vector{Float64} = get!(task_local_storage(), :CEG_buffer_interpolation) do
            Vector{Float64}(undef, 64)
        end
        @inbounds mul!(a, COEFF, X)
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
            g.grid[x0,y0,z0,1]*mrx*mry*mrz + g.grid[x1,y0,z0,1]* rx*mry*mrz +
            g.grid[x0,y1,z0,1]*mrx* ry*mrz + g.grid[x0,y0,z1,1]*mrx*mry* rz +
            g.grid[x1,y1,z0,1]* rx* ry*mrz + g.grid[x1,y0,z1,1]* rx*mry* rz +
            g.grid[x0,y1,z1,1]*mrx* ry* rz + g.grid[x1,y1,z1,1]* rx* ry* rz
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
    charges::Vector{Float64}
    grids::Vector{EnergyGrid}
    atomsidx::Vector{Int} # index of the grid corresponding to the atom
    ewald::EwaldFramework
    forcefield::ForceField
    block::BlockFile
end

"""
    energy_point(setup::CrystalEnergySetup, positions, ctx=nothing)

Compute the energy of a guest molecule whose atoms are at the given `positions` in the
framework. `positions` is a list of triplet of coordinates (with their unit).

Return a pair `(vdw, coulomb)` where `vdw` is the Van der Waals contribution to the energy
and `coulomb` is the electrostatic ones, both in K.

If one of the atom is on a blocking sphere, return `(1e100u"K", 0.0u"K")`

!!! warning
    No validation is used to ensure that the input `positions` are consistent with the
    shape of the molecule.
"""
function energy_point(setup::CrystalEnergySetup, positions, ctx=nothing)
    for pos in positions
        setup.block[pos] && return (1e100u"K", 0.0u"K")
    end
    num_atoms = length(setup.atomsidx)
    vdw = sum(interpolate_grid(setup.grids[setup.atomsidx[i]], positions[i]) for i in 1:num_atoms)
    setup.coulomb.ewald_precision == -Inf && return vdw, 0.0u"K"
    coulomb_direct = sum(setup.charges[i]*interpolate_grid(setup.coulomb, positions[i]) for i in 1:num_atoms)
    coulomb_reciprocal = if ctx isa Nothing
        newmolecule = ChangePositionSystem(setup.molecule, positions)
        compute_ewald(setup.ewald, ((newmolecule,),))
    else
        move_one_system!(ctx, 1, positions)
        compute_ewald(ctx)
    end
    return (vdw, coulomb_direct + coulomb_reciprocal)
end


"""
    energy_grid(setup::Union{CrystalEnergySetup,MonteCarloSetup}, step, num_rotate=40)

Compute the energy on the the system given by `setup` on a regular grid with the given
`step` (with its unit).

If the guest molecule is not monoatomic, the first axe of the returned grid will represent
the rotation angle of the molecule, and its size will be greater or equal to `num_rotate`.
Otherwise (or if `num_rotate == 0`), the first axe can be dropped through`dropdims`.

A value of 1e100 in the grid indicates an inaccessible point due to blocking spheres.

If `setup` is not a [`CrystalEnergyGris`](@ref) but a [`MonteCarloSetup`](@ref), the guest
molecule will be given by the first species kind of `setup`. Existing molecules will be
retained while building the grid.
"""
function energy_grid(setup, step, num_rotate=40)
    setup::Union{CrystalEnergySetup,MonteCarloSetup}
    molecule = setup isa CrystalEnergySetup ? setup.molecule : NanoSystem(setup, 1)
    store = string(hash_string(setup), '-', ustrip(u"Å", step), '-', lebedev_num(molecule, num_rotate))
    path = joinpath(scratchspace, store)
    if isfile(path)
        printstyled("Retrieved energy grid at ", path, ".\n"; color=:cyan)
        return deserialize(path)::Array{Float64,4}
    end

    __pos = position(molecule) / u"Å"
    per_task = if setup isa MonteCarloSetup
        mc0 = MonteCarloSetup(setup; parallel=false) # use a copy
        add_one_system!(mc0, 1) # the newly added molecule will be the probe
        @assert mc0.ewald.last[] ≤ 0 # 0 for neutral species, -1 otherwise
        mc_tasks = [(MonteCarloSetup(mc0), Inf*u"K", NaN*u"K") for _ in 1:(nthreads()-1)]
        push!(mc_tasks, (mc0, Inf*u"K", NaN*u"K"))
        mc_tasks
    else
        ctx0 = EwaldContext(setup.ewald, ((molecule,),))
        ctxs = EwaldContext[deepcopy(ctx0) for _ in 1:(nthreads()-1)]
        push!(ctxs, ctx0)
        pre_printwarn = PRINT_CHARGE_WARNING[]
        # do one empty computation to trigger the warnings if need be
        if pre_printwarn
            energy_point(setup, [SVector{3}(p)*u"Å" for p in __pos])
            PRINT_CHARGE_WARNING[] = false
        end
        ctxs
    end

    printstyled("Creating energy grid at ", path, " ... "; color=:yellow)
    axeA, axeB, axeC = setup isa CrystalEnergySetup ? bounding_box(setup.framework) : eachcol(setup.speciesblocks[1].csetup.cell.mat)
    numA = floor(Int, norm(axeA) / step) + 1
    numB = floor(Int, norm(axeB) / step) + 1
    numC = floor(Int, norm(axeC) / step) + 1
    stepA = axeA / numA
    stepB = axeB / numB
    stepC = axeC / numC
    rotpos::Vector{Vector{SVector{3,Float64}}} = if num_rotate == 0 || length(molecule) == 1
        [[SVector{3}(p) for p in __pos]]
    else
        rots, _ = get_rotation_matrices(molecule, num_rotate)
        [[SVector{3}(r*p) for p in __pos] for r in rots]
    end
    allvals = Array{Float64}(undef, length(rotpos), numA, numB, numC)
    allvals .= NaN

    @loadbalance 0.4 for iABC in CartesianIndices((numA, numB, numC))
        iA, iB, iC = Tuple(iABC)
        thisofs = (iA-1)*stepA + (iB-1)*stepB + (iC-1)*stepC
        bufferpos = Vector{SVector{3,TÅ}}(undef, length(__pos))
        if setup isa MonteCarloSetup
            mc, current, before = per_task[taskid]
        else
            ctx = per_task[taskid]
        end
        for (k, pos) in enumerate(rotpos)
            ofs = if num_rotate < 0
                rand()*numA*stepA + rand()*numB*stepB + rand()*numC*stepC
            else
                thisofs
            end
            bufferpos .= SVector{3,TÅ}.((ofs,) .+ pos.*u"Å")
            newval = ustrip(u"K", if setup isa MonteCarloSetup
                _newval, current, before = energy_point_mc!(mc, bufferpos, current, before)
                per_task[taskid] = (mc, current, before)
                _newval
            else
                sum(energy_point(setup, bufferpos, ctx))
            end)
            allvals[k,iA,iB,iC] = newval
        end
    end
    setup isa CrystalEnergySetup && pre_printwarn && (PRINT_CHARGE_WARNING[] = pre_printwarn)
    serialize(path, allvals)
    printstyled("Done.\n"; color=:cyan)
    return allvals
end
