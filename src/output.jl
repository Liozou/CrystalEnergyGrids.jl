using Printf

@noinline opaque_job(@nospecialize(x)) = Base.inferencebarrier(x)+1
@noinline function opaque_busywait()
    for i in 1:100
        opaque_job(i)
    end
end


function output_restart(path, o::SimulationStep, (a, b, c), (α, β, γ), molnames)
    positions = [[[rand(SVector{3,Float64})]] for _ in 1:o.len]
    invmat = inv(o.mat)
    open(path, "w") do io
        for s in ("unit-cell-vector-", "cell-vector-")
            for (i, x) in enumerate(('a', 'b', 'c'))
                @printf io "%19.12f%19.12f%19.12f\n" o.mat[1,i] o.mat[2,i] o.mat[3,i]
            end
        end
        posi = positions[end]
        for (j, poss) in enumerate(posi), (k, opos) in enumerate(poss)
            abc = invmat * opos
            pos = o.mat*abc
            @printf io "Adsorbate-atom-position: %19.12f%19.12f%19.12f\n" pos[1] pos[2] pos[3]
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
