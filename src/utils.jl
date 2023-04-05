# General utils

function wrap_atom(point, mat, invmat)
    abc = invmat * (point isa SVector ? point : SVector{3,Float64}(point))
    # return mat * (abc .- floor.(abc))
    x, y, z = mat * (abc .- floor.(abc))
    ε = 1e-14
    return SVector{3,Float64}(clamp(x, ε, mat[1,1]-ε), clamp(y, ε, mat[2,2]-ε), clamp(z, ε, mat[3,3]-ε))
end

function offsetpoint(point, mat, invmat, shift, size, dims)
    newpoint = wrap_atom(point, mat, invmat)
    @. (newpoint - shift)/size*dims + 1
end

function inverse_abc_offsetpoint(ipoint, invmat, shift, size, dims)
    abc = invmat*(@. (ipoint - 1)/dims*size + shift)
    abc .- floor.(abc)
end

# The following are copied from PeriodicGraphEmbeddings.jl

function cell_parameters(mat::AbstractMatrix)
    _a, _b, _c = eachcol(mat)
    a = norm(_a)
    b = norm(_b)
    c = norm(_c)
    α = acosd(_b'_c/(b*c))
    β = acosd(_c'_a/(c*a))
    γ = acosd(_a'_b/(a*b))
    (a, b, c), (α, β, γ)
end

function prepare_periodic_distance_computations(mat)
    (a, b, c), (α, β, γ) = cell_parameters(mat)
    ortho = all(x -> isapprox(Float16(x), 90; rtol=0.02), (α, β, γ))
    _a, _b, _c = eachcol(mat)
    safemin = min(Float64(dot(cross(_b, _c), _a)/(b*c)),
                  Float64(dot(cross(_c, _a), _b)/(a*c)),
                  Float64(dot(cross(_a, _b), _c)/(a*b)))/2
    # safemin is the half-distance between opposite planes of the unit cell
    return MVector{3,Float64}(undef), ortho, safemin
end

function periodic_distance!(buffer, mat, ortho, safemin)
    @simd for i in 1:3
        diff = buffer[i] + 0.5
        buffer[i] = diff - floor(diff) - 0.5
    end
    ref = norm(mat*buffer)
    (ortho || ref ≤ safemin) && return ref
    @inbounds for i in 1:3
        buffer[i] += 1
        newnorm = norm(mat*buffer)
        newnorm < ref && return newnorm # in a reduced lattice, there should be at most one
        buffer[i] -= 2
        newnorm = norm(mat*buffer)
        newnorm < ref && return newnorm
        buffer[i] += 1
    end
    return ref
end

function fraction_sites(egrid)
    bins = zeros(100)
    for val in egrid
        -10000 < val < 0 || continue
        bins[1+floor(Int, -val/100)] += 1
    end
    s = sum(bins)
    @show s, s/length(egrid)
    bins ./= s
    bins
end
