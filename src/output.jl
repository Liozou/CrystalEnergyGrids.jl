using Printf



function output_restart(path, o::SimulationStep, (a, b, c), (α, β, γ), molnames)
    positions = [Vector{SVector{3,TÅ}}[] for _ in o.positions]
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
output_restart(path, mc::MonteCarloSetup) = output_restart(path, deepcopy(mc.step))
