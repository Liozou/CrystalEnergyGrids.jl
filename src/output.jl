using Printf

struct OutputSimulationStep
    cell::CellMatrix
    ff::ForceField
    positions::Vector{SVector{3,typeof(1.0u"Å")}}
    nummol::Vector{Int}
    idx::Vector{Vector{Int}}
end

function OutputSimulationStep(mc::MonteCarloSimulation)
    positions = Vector{SVector{3,typeof(1.0u"Å")}}(undef, mc.numatoms[])
    nummol = Vector{Int}(undef, length(mc.positions))
    k = 0
    for (i, positioni) in enumerate(mc.positions)
        nummol[i] = length(positioni)
        for poss in positioni
            for pos in poss
                k += 1
                positions[k] = pos
            end
        end
    end
    OutputSimulationStep(mc.cell, mc.ff, positions, nummol, mc.idx)
end

# Signal that the channel should be closed
function OutputSimulationStep(mc::MonteCarloSimulation, ::Nothing)
    OutputSimulationStep(mc.cell, mc.ff, SVector{3,typeof(1.0u"Å")}[], Int[], Vector{Int}[])
end

function output_pdb(file, o::OutputSimulationStep, (a, b, c), (α, β, γ), i)
    open(file, "a") do io
        @printf io "MODEL %4d\n" i
        @printf io "CRYST1%9g%9g%9g%7g%7g%7g\n" NoUnits(a/u"Å") NoUnits(b/u"Å") NoUnits(c/u"Å") α β γ
        k = 0
        for (idxi, nummoli) in zip(o.idx, o.nummol), _ in 1:nummoli, ix in idxi
            symb = String(o.ff.symbols[ix])
            k += 1
            abc = o.cell.invmat * o.positions[k]
            pos = o.cell.mat*(abc .- floor.(abc))/u"Å"
            @printf io "ATOM  %5d %4.4s MOL          " k symb
            @printf io "%8.4lf%8.4lf%8.4lf" pos[1] pos[2] pos[3]
            @printf io "  1.00  0.00          %2.2s  \n" symb
        end
        @printf io "ENDMDL\n"
    end
end

function pdb_output_handler(file, cell::CellMatrix)
    lengths, angles = cell_parameters(cell.mat)
    Channel{OutputSimulationStep}(100) do channel
        for (i, o) in enumerate(channel)
            isempty(o.idx) && break
            output_pdb(file, o, lengths, angles, i)
        end
        nothing
    end
end