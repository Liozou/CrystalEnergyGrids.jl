using Printf

export output_pdb, output_restart

SimulationStep(mc::MonteCarloSetup) = SimulationStep(mc.step, :output)

function output_pdb(path, o::SimulationStep, (a, b, c), (α, β, γ), i, atomcounter)
    invmat = inv(ustrip.(u"Å", o.mat))*u"Å^-1"
    open(path, "a") do io
        @printf io "MODEL %4d\n" i
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
    open(path, "w") do io
        println(io, """Cell info:
========================================================================
number-of-unit-cells: 1 1 1""")
        for s in ("unit-cell-vector-", "cell-vector-")
            for (i, x) in enumerate(('a', 'b', 'c'))
                @printf io "%19.12f%19.12f%19.12f\n" NoUnits(o.mat[1,i]/u"Å") NoUnits(o.mat[2,i]/u"Å") NoUnits(o.mat[3,i]/u"Å")
            end
        end
        posi = positions[end]
        println(io, "------------------------------------------------------------------------")
        for (j, poss) in enumerate(posi), (k, opos) in enumerate(poss)
            abc = invmat * opos
            pos = NoUnits.(o.mat*(abc .- floor.(abc))./u"Å")
            @printf io "Adsorbate-atom-position: %d %d%19.12f%19.12f%19.12f\n" (j-1) (k-1) pos[1] pos[2] pos[3]
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

function pdb_output_handler(path, mat::SMatrix{3,3,TÅ,9})
    taskref = Ref{Task}()
    if isempty(path)
        Channel{SimulationStep}(Inf; taskref) do channel
            while !isempty(take!(channel).ffidx) end # skip until isempty(o.ffidx)
        end
    else
        lengths, angles = cell_parameters(mat)
        atomcounter = Counter3D()
        Channel{SimulationStep}(Inf; taskref) do channel
            i = 0
            while true
                i += 1
                o = take!(channel)
                isempty(o.ffidx) && break
                output_pdb(path, o, lengths, angles, i, atomcounter)
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
