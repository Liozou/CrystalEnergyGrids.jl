using Printf

export output_pdb, output_restart

SimulationStep(mc::MonteCarloSetup) = SimulationStep(mc.step, :output)

function output_pdb(path, o::Union{ProtoSimulationStep,SimulationStep}, (a, b, c), (α, β, γ), i, atomcounter)
    invmat = inv(ustrip.(u"Å", o.mat))*u"Å^-1"
    open(path, "a") do io
        @printf io "MODEL     %4d\n" i
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

function output_handler(outdir::AbstractString, outtype::Vector{Symbol}, step::SimulationStep)
    taskref = Ref{Task}()
    if isempty(outdir) || length(outtype)  == (:energies in outtype)
        Channel{SimulationStep}(Inf; taskref) do channel
            while !isempty(take!(channel).ffidx) end # skip until isempty(o.ffidx)
        end
    else
        lengths, angles = cell_parameters(step.mat)
        atomcounter = Counter3D()
        stream_path = :stream in outtype ? joinpath(outdir, "steps.stream") : ""
        isempty(stream_path) || save_init(stream_path, step)
        zst_path = :zst in outtype ? joinpath(outdir, "steps.zst") : ""
        isempty(zst_path) || save_init(zst_path, step)
        pdb_path = :pdb in outtype ? joinpath(outdir, "positions.pdb") : ""
        Channel{SimulationStep}(Inf; taskref) do channel
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
