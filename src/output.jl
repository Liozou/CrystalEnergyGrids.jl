using Printf, Serialization

export output_pdb, output_restart

struct OutputSimulationStep
    cell::CellMatrix
    ff::ForceField
    positions::Vector{SVector{3,typeof(1.0u"Å")}}
    nummol::Vector{Int}
    ffidx::Vector{Vector{Int}}
end

function OutputSimulationStep(mc::MonteCarloSetup)
    positions = Vector{SVector{3,typeof(1.0u"Å")}}(undef, mc.numatoms[])
    nummol = Vector{Int}(undef, length(mc.step.positions))
    k = 0
    for (i, positioni) in enumerate(mc.step.positions)
        nummol[i] = length(positioni)
        for poss in positioni
            for pos in poss
                k += 1
                positions[k] = pos
            end
        end
    end
    OutputSimulationStep(mc.step.cell, mc.step.ff, positions, nummol, mc.step.ffidx)
end

# Signal that the channel should be closed
function OutputSimulationStep(mc::MonteCarloSetup, ::Nothing)
    OutputSimulationStep(mc.step.cell, mc.step.ff, SVector{3,typeof(1.0u"Å")}[], Int[], Vector{Int}[])
end

function output_pdb(path, o::OutputSimulationStep, (a, b, c), (α, β, γ), i)
    open(path, "a") do io
        @printf io "MODEL %4d\n" i
        @printf io "CRYST1%9g%9g%9g%7g%7g%7g\n" NoUnits(a/u"Å") NoUnits(b/u"Å") NoUnits(c/u"Å") α β γ
        k = 0
        for (ffidxi, nummoli) in zip(o.ffidx, o.nummol), _ in 1:nummoli, ix in ffidxi
            symb = String(o.ff.symbols[ix])
            k += 1
            abc = o.cell.invmat * o.positions[k]
            pos = o.cell.mat*(abc .- floor.(abc))/u"Å"
            @printf io "ATOM  %5d %4.4s MOL          %8.4lf%8.4lf%8.4lf  1.00  0.00          %2.2s  \n" k symb pos[1] pos[2] pos[3] symb
        end
        @printf io "ENDMDL\n"
        nothing
    end
end

"""
    output_pdb(path, mc::MonteCarloSetup, o=OutputSimulationStep(mc), i=0)

Output a .pdb at the given `path` representing the positions of the atoms for the given
output `o`, which corresponds to the `i`-th simulation step of the input `mc`.
If not provided, the output `o` corresponds to the current status of `mc`.

The output is appended to the file at the given `path`, if any.

See also [`output_restart`](@ref).
"""
function output_pdb(path, mc::MonteCarloSetup, o=OutputSimulationStep(mc), i=0)
    lengths, angles = cell_parameters(mc.step.cell.mat)
    output_pdb(path, o, lengths, angles, i)
end


function output_restart(path, o::OutputSimulationStep, (a, b, c), (α, β, γ), molnames, charges)
    open(path, "w") do io
        println(io, """Cell info:
========================================================================
number-of-unit-cells: 1 1 1""")
        for s in ("unit-cell-vector-", "cell-vector-")
            for (i, x) in enumerate(('a', 'b', 'c'))
                @printf io "%s%c:%19.12f%19.12f%19.12f\n" s x NoUnits(o.cell.mat[1,i]/u"Å") NoUnits(o.cell.mat[2,i]/u"Å") NoUnits(o.cell.mat[3,i]/u"Å")
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


Components: """, length(o.ffidx), " (Adsorbates ", sum(o.nummol), """, Cations 0)
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
        t = 0
        for (i, (name, num)) in enumerate(zip(molnames, o.nummol))
            ffidxi = o.ffidx[i]
            m = length(ffidxi)
            println(io, "Component: ", i-1, "     Adsorbate    ", num, " molecules of ", name)
            println(io, "------------------------------------------------------------------------")
            for j in 0:num-1, k in 0:m-1
                t += 1
                abc = o.cell.invmat * o.positions[t]
                pos = NoUnits.(o.cell.mat*(abc .- floor.(abc))./u"Å")
                @printf io "Adsorbate-atom-position: %d %d%19.12f%19.12f%19.12f\n" j k pos[1] pos[2] pos[3]
            end
            for s in ("velocity:", "force:   "), j in 0:num-1, k in 0:m-1
                @printf io "Adsorbate-atom-%s %d %d     0.000000000000     0.000000000000     0.000000000000\n" s j k
            end
            for j in 0:num-1, (k, ix) in enumerate(ffidxi)
                @printf io "Adsorbate-atom-charge:   %d %d%19.12f\n" j (k-1) NoUnits(charges[ix]/u"e_au")
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
    function output_restart(path, mc::MonteCarloSetup, o=OutputSimulationStep(mc))

Output a restart file compatible with RASPA at the given `path` for the given output `o`,
which represents a simulation step of the input `mc`.
If not provided, the output `o` corresponds to the current status of `mc`.

See also [`output_pdb`](@ref).
"""
function output_restart(path, mc::MonteCarloSetup, o=OutputSimulationStep(mc))
    lengths, angles = cell_parameters(mc.step.cell.mat)
    molnames = [identify_molecule([mc.step.ff.symbols[ix] for ix in ffidxi]) for ffidxi in mc.step.ffidx]
    output_restart(path, o, lengths, angles, molnames, mc.step.charges)
end

function pdb_output_handler(path, cell::CellMatrix)
    taskref = Ref{Task}()
    if isempty(path)
        Channel{OutputSimulationStep}(1; taskref) do channel
            for o in channel
                isempty(o.ffidx) && break
            end
        end
    else
        lengths, angles = cell_parameters(cell.mat)
        Channel{OutputSimulationStep}(100; taskref) do channel
            for (i, o) in enumerate(channel)
                isempty(o.ffidx) && break
                output_pdb(path, o, lengths, angles, i)
            end
            nothing
        end
    end, taskref
end
