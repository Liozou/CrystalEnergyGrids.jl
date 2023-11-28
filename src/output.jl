using Printf

export output_pdb, output_restart

SimulationStep(mc::MonteCarloSetup) = SimulationStep(mc.step, :output)

function output_pdb(path, o::SimulationStep, (a, b, c), (α, β, γ), i, atomcounter)
    invmat = inv(ustrip.(u"Å", o.mat))*u"Å^-1"
    open(path, "a") do io
        @printf io "MODEL %4d\n" i
        @printf io "CRYST1%9g%9g%9g%7g%7g%7g\n" NoUnits(a/u"Å") NoUnits(b/u"Å") NoUnits(c/u"Å") α β γ
        for (l, opos) in enumerate(o.positions)
            abc = invmat * opos
            pos = o.mat*(abc .- floor.(abc))/u"Å"
            symb = "Na"
            molid, atomid = atomcounter[1, l, 1]
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
    positions = [Vector{SVector{3,TÅ}}[] for _ in o.positions]
    invmat = inv(ustrip.(u"Å", o.mat))*u"Å^-1"
    for (l, pos) in enumerate(o.positions)
        posi = positions[1]
        length(posi) < l && resize!(posi, l)
        isassigned(posi, l) || (posi[l] = SVector{3,TÅ}[])
        poss = posi[l]
        length(poss) < 1 && resize!(poss, 1)
        poss[1] = pos
    end
    open(path, "w") do io
        for s in ("unit-cell-vector-", "cell-vector-")
            for (i, x) in enumerate(('a', 'b', 'c'))
                @printf io "%19.12f%19.12f%19.12f\n" NoUnits(o.mat[1,i]/u"Å") NoUnits(o.mat[2,i]/u"Å") NoUnits(o.mat[3,i]/u"Å")
            end
        end
        posi = positions[end]
        for (j, poss) in enumerate(posi), (k, opos) in enumerate(poss)
            abc = invmat * opos
            pos = NoUnits.(o.mat*(abc .- floor.(abc))./u"Å")
            @printf io "Adsorbate-atom-position: %d %d%19.12f%19.12f%19.12f\n" (j-1) (k-1) pos[1] pos[2] pos[3]
        end
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
    molnames = ["Na" for _ in step.positions]
    output_restart(path, step, lengths, angles, molnames)
end
output_restart(path, mc::MonteCarloSetup) = output_restart(path, SimulationStep(mc))

function pdb_output_handler(path, mat::SMatrix{3,3,TÅ,9})
    taskref = Ref{Task}()
    if isempty(path)
        Channel{SimulationStep}(Inf; taskref) do channel
            while true
                take!(channel)
            end
        end
    else
        lengths, angles = cell_parameters(mat)
        atomcounter = Counter3D()
        Channel{SimulationStep}(Inf; taskref) do channel
            i = 0
            while true
                i += 1
                o = take!(channel)
                output_pdb(path, o, lengths, angles, i, atomcounter)
            end
            nothing
        end
    end, taskref
end
