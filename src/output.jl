using Printf

export output_pdb, output_restart

SimulationStep(mc::MonteCarloSetup) = SimulationStep(mc.step, :output)

function output_pdb(path, o::Union{ProtoSimulationStep,SimulationStep}, (a, b, c), (α, β, γ), modelidx, atomcounter)
    invmat = inv(ustrip.(u"Å", o.mat))*u"Å^-1"
    open(path, "a") do io
        @printf io "MODEL     %4d\n" modelidx
        @printf io "CRYST1%9g%9g%9g%7g%7g%7g\n" NoUnits(a/u"Å") NoUnits(b/u"Å") NoUnits(c/u"Å") α β γ
        for (l, opos) in enumerate(o.positions)
            i, j, k = o.atoms[l]
            ix = o.ffidx[i][k]
            abc = invmat * opos
            pos = o.mat*(abc .- floor.(abc))/u"Å"
            symb = String(o.ff.symbols[ix])
            molid, atomid = atomcounter[i, j, k]
            @printf io "ATOM  %-6d%4.4s MOL  %-8d%8.4lf%8.4lf%8.4lf  1.00  0.00          %2.2s  \n" atomid symb molid pos[1] pos[2] pos[3] symb
        end
        @printf io "ENDMDL\n"
        nothing
    end
end

function output_pdb(path, sh::SiteHopping, (a, b, c), (α, β, γ), modelidx, atomcounter)
    open(path, "a") do io
        @printf io "MODEL     %4d\n" modelidx
        @printf io "CRYST1%9g%9g%9g%7g%7g%7g\n" NoUnits(a/u"Å") NoUnits(b/u"Å") NoUnits(c/u"Å") α β γ
        for i in sh.population, (k, pos) in enumerate(sh.sites[i])
            symb = String(sh.atomnames[k])
            molid, atomid = atomcounter[1, i, k]
            @printf io "ATOM  %-6d%4.4s MOL  %-8d%8.4lf%8.4lf%8.4lf  1.00  0.00          %2.2s  \n" atomid symb molid NoUnits(pos[1]/u"Å") NoUnits(pos[2]/u"Å") NoUnits(pos[3]/u"Å") symb
        end
        @printf io "ENDMDL\n"
        nothing
    end
end


function output_pdb(path, o::Union{ProtoSimulationStep,SimulationStep,SiteHopping})
    if splitext(path)[2] != ".pdb"
        path = path*".pdb"
    end
    lengths, angles = cell_parameters(o.mat)
    atomcounter = Counter3D()
    output_pdb(path, o, lengths, angles, 0, atomcounter)
end

"""
    output_pdb(path, mc::MonteCarloSetup, o=SimulationStep(mc), i=0)

Output a .pdb at the given `path` representing the positions of the atoms for the given
output `o`, which corresponds to the `i`-th simulation step of the input `mc`.
If not provided, the output `o` corresponds to the current status of `mc`.

The output is appended to the file at the given `path`, if any.

See also [`output_restart`](@ref).
"""
function output_pdb(path, mc::MonteCarloSetup, o=SimulationStep(mc), i=0, atomcounter=Counter3D())
    lengths, angles = cell_parameters(o.mat)
    output_pdb(path, o, lengths, angles, i, atomcounter)
end


function output_restart(path, o::SimulationStep, (a, b, c), (α, β, γ), molnames)
    positions = [Vector{SVector{3,TÅ}}[] for _ in o.ffidx]
    invmat = inv(ustrip.(u"Å", o.mat))*u"Å^-1"
    for (l, pos) in enumerate(o.positions)
        i, j, k = o.atoms[l]
        posi = positions[i]
        length(posi) < j && resize!(posi, j)
        isassigned(posi, j) || (posi[j] = SVector{3,TÅ}[])
        poss = posi[j]
        length(poss) < k && resize!(poss, k)
        poss[k] = pos
    end
    for posi in positions
        hasassigned = [j for j in 1:length(posi) if isassigned(posi, j)]
        keepat!(posi, hasassigned)
    end
    open(path, "w") do io
        println(io, """Cell info:
========================================================================
number-of-unit-cells: 1 1 1""")
        for s in ("unit-cell-vector-", "cell-vector-")
            for (i, x) in enumerate(('a', 'b', 'c'))
                @printf io "%s%c:%19.12f%19.12f%19.12f\n" s x NoUnits(o.mat[1,i]/u"Å") NoUnits(o.mat[2,i]/u"Å") NoUnits(o.mat[3,i]/u"Å")
            end
            println(io)
        end
        @printf io "cell-lengths:%19.12f%19.12f%19.12f\n" NoUnits(a/u"Å") NoUnits(b/u"Å") NoUnits(c/u"Å")
        @printf io "cell-angles:%19.12f%19.12f%19.12f\n\n" α β γ
        println(io, """

Maximum changes for MC-moves:
========================================================================
Maximum-volume-change: 0.000500
Maximum-Gibbs-volume-change: 0.000500
Maximum-box-shape-change: 0.100000 0.100000 0.100000, 0.100000 0.100000 0.100000, 0.100000 0.100000 0.100000


Acceptance targets for MC-moves:
========================================================================
Target-volume-change: 0.500000
Target-box-shape-change: 0.500000
Target-Gibbs-volume-change: 0.500000


Components: """, length(o.ffidx), " (Adsorbates ", sum(length, positions), """, Cations 0)
========================================================================""")
        for (i, name) in enumerate(molnames)
            im1 = i-1
            println(io, "Component ", im1, " (", name, ')')
            println(io, "	Fractional-molecule-id component ", im1, ": -1")
            print(io, "	Lambda-factors component ", im1, ": ")
            join(io, " 0.000000" for _ in o.ffidx[i])
            println(io)
            println(io, "	Number-of-biasing-factors component ", im1, ": 21")
            println(io, "	Biasing-factors component ", im1, ":  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000")
            println(io, "	Maximum-CF-Lambda-change component ", im1, ": 0.500000")
            println(io, "	Maximum-CBCF-Lambda-change component ", im1, ": 0.500000\n")
            println(io, "	Maximum-translation-change component ", im1, ": 0.246897,0.246897,0.246897")
            println(io, "	Maximum-translation-in-plane-change component ", im1, ": 0.000000,0.000000,0.000000")
            println(io, "	Maximum-rotation-change component ", im1, ": 0.123456 0.123456 0.123456")
        end
        println(io, "\nReactions: 0\n")
        for (i, name) in enumerate(molnames)
            posi = positions[i]
            num = length(posi)
            ffidxi = o.ffidx[i]
            m = length(ffidxi)
            println(io, "Component: ", i-1, "     Adsorbate    ", num, " molecules of ", name)
            println(io, "------------------------------------------------------------------------")
            for (j, poss) in enumerate(posi), (k, opos) in enumerate(poss)
                abc = invmat * opos
                pos = NoUnits.(o.mat*(abc .- floor.(abc))./u"Å")
                @printf io "Adsorbate-atom-position: %d %d%19.12f%19.12f%19.12f\n" (j-1) (k-1) pos[1] pos[2] pos[3]
            end
            for s in ("velocity:", "force:   "), j in 0:num-1, k in 0:m-1
                @printf io "Adsorbate-atom-%s %d %d     0.000000000000     0.000000000000     0.000000000000\n" s j k
            end
            for j in 0:num-1, (k, ix) in enumerate(ffidxi)
                @printf io "Adsorbate-atom-charge:   %d %d%19.12f\n" j (k-1) NoUnits(o.charges[ix]/u"e_au")
            end
            for j in 0:num-1, k in 0:m-1
                @printf io "Adsorbate-atom-scaling:  %d %d     1.000000000000\n" j k
            end
            for j in 0:num-1, k in 0:m-1
                @printf io "Adsorbate-atom-fixed:    %d %d       0 0 0\n" j k
            end
            println(io)
        end
        println(io, '\n')
        nothing
    end
end

"""
    function output_restart(path, mc::MonteCarloSetup)
    function output_restart(path, step::SimulationStep(mc))

Output a restart file compatible with RASPA at the given `path` for the given output `step`,
which represents a simulation step of the input `mc`.
If not provided, the output `step` corresponds to the current status of `mc`.

See also [`output_pdb`](@ref).
"""
function output_restart(path, step::SimulationStep)
    lengths, angles = cell_parameters(step.mat)
    molnames = [identify_molecule([step.ff.symbols[ix] for ix in ffidxi]) for ffidxi in step.ffidx]
    output_restart(path, step, lengths, angles, molnames)
end
output_restart(path, mc::MonteCarloSetup) = output_restart(path, SimulationStep(mc))

function output_handler(outdir::AbstractString, outtype::Vector{Symbol}, step::SimulationStep, initial::Bool)
    taskref = Ref{Task}()
    if isempty(outdir) || !any(x -> (y = String(x); !endswith(y, "energies") && (!initial ⊻ startswith(y, "initial"))), outtype)
        Channel{SimulationStep}(1; taskref) do channel
            while !isempty(take!(channel).ffidx) end # skip until isempty(o.ffidx)
        end
    else
        lengths, angles = cell_parameters(step.mat)
        atomcounter = Counter3D()
        prefix = initial ? :initial_ : Symbol("")
        stream_path = Symbol(prefix, :stream) in outtype ? joinpath(outdir, string(prefix, "steps.stream")) : ""
        isempty(stream_path) || save_init(stream_path, step)
        zst_path = Symbol(prefix, :zst) in outtype ? joinpath(outdir, string(prefix, "steps.zst")) : ""
        isempty(zst_path) || save_init(zst_path, step)
        pdb_path = Symbol(prefix, :pdb) in outtype ? joinpath(outdir, string(prefix, "trajectory.pdb")) : ""
        Channel{SimulationStep}(1; taskref) do channel
            i = 0
            while true
                i += 1
                o = take!(channel)
                isempty(o.ffidx) && break
                isempty(zst_path) || save(zst_path, o)
                isempty(stream_path) || save(stream_path, o)
                isempty(pdb_path) || output_pdb(pdb_path, o, lengths, angles, i, atomcounter)
            end
            nothing
        end
    end, taskref
end

function output_cube(path, grid::Array{Float64,4}, framework, T=300)
    output_cube(path, meanBoltzmann(grid, T), framework)
end

function output_cube(path, grid::Array{Float64,3}, framework)
    cif = ispath(framework) ? framework : joinpath(RASPADIR[], "structures", "cif", framework*".cif")
    system = load_system(AtomsIO.ChemfilesParser(), cif)
    box = AtomsBase.bounding_box(system)./u"bohr"
    atoms = [(a, NoUnits.(p./u"bohr")) for (a, p) in zip(AtomsBase.atomic_number(system), AtomsBase.position(system))]
    open(path, "w") do io
        println(io, "CPMD CUBE FILE")
        println(io, "exported by CrystalEnergyGrids.jl")
        @printf io "%5d %12.6g %12.6g %12.6g\n" length(atoms) 0.0 0.0 0.0
        for (i, l) in enumerate(box)
            a = size(grid,i)
            @printf io "%5d %12.6g %12.6g %12.6g\n" size(grid, i) NoUnits(l[1])/a NoUnits(l[2])/a NoUnits(l[3])/a
        end
        for (a, pos) in atoms
            @printf io "%5d 0.0 %12.6g %12.6g %12.6g\n" a pos[1] pos[2] pos[3]
        end
        counter = 0
        for i1 in axes(grid, 1), i2 in axes(grid, 2), i3 in axes(grid, 3)
            counter += 1
            @printf io "%12.6g" exp(-grid[i1,i2,i3]/300)
            if counter%6==0
                println(io)
            else
                print(io, ' ')
            end
        end
        println(io)
    end
end

# Utility to export the detected local clusters. Visualize in VMD with Points, coloring on ResID
function output_basins_pdb(path, points, nodes, mat)
    (a, b, c), (α, β, γ) = cell_parameters(mat)
    na, nb, nc = size(nodes)
    ϵ = min(ustrip(a)/na, ustrip(b)/nb, ustrip(c)/nc)/2
    counter = 0
    open(path, "w") do io
        @printf io "MODEL     %4d\n" 0
        @printf io "CRYST1%9g%9g%9g%7g%7g%7g\n" ustrip(a) ustrip(b) ustrip(c) α β γ
        for (i,j,k) in points
            px, py, pz = mat * SVector{3,Float64}(i/na, j/nb, k/nc)
            basins = nodes[mod1(i, na),mod1(j, nb),mod1(k, nc)]
            for basin in basins
                counter = mod(counter+1, 100000)
                @printf io "ATOM  %-6d%4.4s MOL  %-8d%8.4lf%8.4lf%8.4lf  1.00  0.00          %2.2s  \n" counter :X basin (px+ϵ*(rand()-0.5)) (py+ϵ*(rand()-0.5)) (pz+ϵ*(rand()-0.5)) :Y
            end
        end
    end
end


function output_density(path, density, mat; mode=:vtf, factor=1000)
    mode === :vtf || mode === :pdb || throw(ArgumentError(lazy"Only mode :vtf and :pdb are supported, not $vtf"))
    _mat = unit(eltype(mat)) === NoUnits ? mat : ustrip.(u"Å", mat)
    (a, b, c), (α, β, γ) = cell_parameters(_mat)
    na, nb, nc = size(density)
    ϵ = min(ustrip(a)/na, ustrip(b)/nb, ustrip(c)/nc)/2
    counter = 0
    _path = endswith(path, string('.', mode)) ? path : string(path, '.', mode)
    nums = [min(30, floor(Int, factor*density[i,j,k]*(1+randexp()))) for k in 1:nc, j in 1:nb, i in 1:na]
    tot = sum(nums)
    @show tot
    open(_path, "w") do io
        if mode === :vtf
            @printf io "atom default radius 0 name X\n"
            @printf io "atom 1:%d resid 1 radius 0 name X\n" tot
            @printf io "unitcell %9g%9g%9g%7g%7g%7g\n" ustrip(a) ustrip(b) ustrip(c) α β γ
            println(io, "timestep")
        elseif mode === :pdb
            @printf io "MODEL     %4d\n" 0
            @printf io "CRYST1%9g%9g%9g%7g%7g%7g\n" ustrip(a) ustrip(b) ustrip(c) α β γ
        end
        idx = 0
        for k in 1:nc, j in 1:nb, i in 1:na
            idx += 1
            num = nums[idx]
            num == 0 && continue
            px, py, pz = ustrip.(_mat * SVector{3,Float64}(i/na, j/nb, k/nc))
            for _ in 1:num
                counter = mod(counter+1, 100000)
                x = px+ϵ*(2rand()-1)
                y = py+ϵ*(2rand()-1)
                z = pz+ϵ*(2rand()-1)
                if mode === :vtf
                    println(io, x, ' ', y, ' ', z)
                elseif mode === :pdb
                    @printf io "ATOM  %-6d%4.4s MOL  %-8d%8.4lf%8.4lf%8.4lf  1.00  0.00          %2.2s  \n" counter :X 1 x y z :Y
                end
            end
        end
    end
end

struct MismatchedUnitCells <: Exception
    c1::SMatrix{3,3,TÅ,9}
    c2::SMatrix{3,3,TÅ,9}
end
function Base.show(io::IO, e::MismatchedUnitCells)
    print(io, "Cannot jointly export cif of a simulation step and a framework with different unit cells: ", ustrip.(u"Å", e.c1), " ≠ ", ustrip.(u"Å", e.c2))
end

"""
    output_cif(path, step::SimulationStep[, framework])

Output `step` as a CIF file at the given `path`.

Only the content of the `SimulationStep` is exported, not the underlying framework.
The optional argument `framework` can be set to an `AtomsBase`-compliant system to export
its content along that of the `step`. It can also be the path to a cif file representing
that system.
"""
function output_cif(path, step::SimulationStep, framework::AbstractSystem=RASPASystem(SA[SA[0.,0.,0.]u"Å",SA[0.,0.,0.]u"Å",SA[0.,0.,0.]u"Å"],SVector{3,TÅ}[],Symbol[],Int[],typeof(1.0u"u")[],Te_au[],false))
    if !all(iszero, bounding_box(framework))
        step.mat ≈ stack(bounding_box(framework)) || throw(MismatchedUnitCells(step.mat, SMatrix{3,3,TÅ,9}(stack(bounding_box(framework)))))
    end

    name = basename(path)
    if splitext(path)[2] != ".cif"
        path = path*".cif"
    end
    mkpath(dirname(path))

    names = Vector{Symbol}(undef, length(step.ff.symbols))
    for (s, idx) in step.ff.sdict
        names[idx] = s
    end
    invmat = inv(ustrip.(u"Å", step.mat))


    open(path, "w") do io
        println(io, "data_", name)
        println(io)
        lengths, angles = cell_parameters(step.mat)
        for (symb, val) in zip((:a, :b, :c), lengths)
            @printf io "_cell_length_%s   %.8g\n" symb ustrip(u"Å", val)
        end
        for (symb, val) in zip((:alpha, :beta, :gamma), angles)
            @printf io "_cell_angle_%s   %.8g\n" symb val
        end
        println(io, """
        _symmetry_space_group_name_H-M	'P 1'
        _symmetry_Int_Tables_number	1
        _symmetry_cell_setting	triclinic

        loop_
        _symmetry_equiv_pos_as_xyz
        x,y,z

        loop_
        _atom_site_label
        _atom_site_type_symbol
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z""")
        counter = 0
        for (ppos, name, num) in zip(position(framework), atomic_symbol(framework), atomic_number(framework))
            counter += 1
            pos = invmat*ustrip.(u"Å", ppos)
            pos -= floor.(pos)
            @printf io "%s\t%s\t%12g\t%12g\t%12g\n" name ATOMS_PER_NUM[num] pos[1] pos[2] pos[3]
        end
        for (i, posidxi) in enumerate(step.posidx)
            ffidxi = step.ffidx[i]
            for (j, posidxij) in enumerate(posidxi), (k, l) in enumerate(posidxij)
                counter += 1
                step.atoms[l] == (i,j,k) || continue
                ix = ffidxi[k]
                name = names[ix]
                pos = invmat*ustrip.(u"Å", step.positions[l])
                pos -= floor.(pos)
                symb = step.ff.symbols[ix]
                @printf io "%s_%i\t%s\t%12g\t%12g\t%12g\n" name counter symb pos[1] pos[2] pos[3]
            end
        end
    end
    nothing
end
function output_cif(path, step::SimulationStep, framework::AbstractString)
    cif = ispath(framework) ? framework : joinpath(RASPADIR[], "structures", "cif", framework*".cif")
    system = load_system(AtomsIO.ChemfilesParser(), cif)
    output_cif(path, step, system)
end

function output_cif(path, framework::AbstractSystem{3})
    mat = stack3(bounding_box(framework))*u"Å"
    a, b, c = perpendicular_lengths(mat)
    zeroff = ForceField(Matrix{Union{InteractionRule,InteractionRuleSum}}(undef, 0, 0), IdDict{Symbol,Int}(), Symbol[], min(a, b, c)/3, "")
    zerostep = SimulationStep(ProtoSimulationStep(zeroff, Te_au[], mat,SVector{3,TÅ}[], false, NTuple{3,Int}[], falses(0), Vector{Int}[]))
    output_cif(path, zerostep, framework)
end
function output_cif(path, framework::AbstractString)
    cif = ispath(framework) ? framework : joinpath(RASPADIR[], "structures", "cif", framework*".cif")
    system = load_system(AtomsIO.ChemfilesParser(), cif)
    output_cif(path, system)
end
