"""
    TailCorrection

Internal representation of the tail correction computation for truncated force field
interactions.
"""
struct TailCorrection
    value::Base.RefValue{TK}
    framework::Vector{TK}
    cross::Matrix{TK}
    numspecies::Vector{Int}
end

Base.getindex(tc::TailCorrection) = tc.value[]

function TailCorrection(ff::ForceField, ffidx::Vector{Vector{Int}}, framework_atoms::Vector{Int}, λ::Float64, numspecies)
    n = length(framework_atoms)
    @assert n == length(ff.sdict)
    allatoms = Int[i for (i,x) in enumerate(framework_atoms) if x > 0]
    for ffidxi in ffidx
        append!(allatoms, ffidxi)
    end
    sort!(allatoms); unique!(allatoms)
    atom_pairs = fill(NaN*u"K", n, n)
    for (idxi, i) in enumerate(allatoms)
        atom_pairs[i,i] = tailcorrection(ff[i,i], ff.cutoff)*λ
        for idxj in (idxi+1):length(allatoms)
            j = allatoms[idxj]
            atom_pairs[j,i] = atom_pairs[i,j] = tailcorrection(ff[i,j], ff.cutoff)*λ
            # λ is 2π / V where V is the volume of the unit cell in Å^3
        end
    end
    m = length(numspecies)
    @assert m == length(ffidx)
    cross = fill(NaN*u"K", m, m)
    framework = Vector{TK}(undef, m)

    # initial value is given by the empty framework
    value = 0.0u"K"
    for (ki, xi) in enumerate(framework_atoms)
        xi == 0 && continue
        for (kj, xj) in enumerate(framework_atoms)
            xj == 0 && continue
            value += xi*xj*atom_pairs[ki,kj]
        end
    end

    for (iA, ffidxiA) in enumerate(ffidx)
        f = 0.0u"K"
        for (k2, x2) in enumerate(framework_atoms)
            x2 == 0 && continue
            for k1 in ffidxiA
                f += atom_pairs[k1,k2]*x2
            end
        end
        framework[iA] = 2f
        numA = numspecies[iA]
        value += 2f*numA
        for iB in iA:m
            c = 0.0u"K"
            ffidxiB = ffidx[iB]
            for kA in ffidxiA
                for kB in ffidxiB
                    c += atom_pairs[kA,kB]
                end
            end
            cross[iA,iB] = cross[iB,iA] = c
            numB = numspecies[iB]
            value += (2 - (iA==iB))*c*numA*numB
        end
    end
    TailCorrection(Ref(value), framework, cross, numspecies)
end

function Base.copy(tc::TailCorrection)
    TailCorrection(Ref(tc.value[]), tc.framework, tc.cross, copy(tc.numspecies))
end

function modify_species_dryrun(tc::TailCorrection, i, num)
    diff = tc.framework[i]
    for (j, other) in enumerate(tc.numspecies)
        if j == i
            diff += (num + 2*other)*tc.cross[i,i]
        else
            diff += 2other*tc.cross[j,i]
        end
    end
    diff * num
end

function modify_species!(tc::TailCorrection, i, num)
    diffnum = modify_species_dryrun(tc, i, num)
    tc.numspecies[i] += num
    tc.value[] += diffnum
end
