using Artifacts
using StaticArrays

const _rotmatrix = SMatrix{3,3,Float64,9}([-0.17963068200890037 -0.21953827352603253 -0.9589242746631385; -0.9599246581752935 0.25228151379218244 0.12206018362173197; 0.21512198564550156 0.9424208106021077 -0.2560577025515984])

const lebedev_path = artifact"lebedev"
const lebedev_sizes = sort!([parse(Int, i) for i in readdir(lebedev_path)])

export get_rotation_matrices

struct LebedevGrid
    size::Int32
    degree::Int32
    weights::Vector{Float64}
    points::Vector{SVector{3,Float64}}
end

function read_lebedev_grid(file, islinearsymmetric)
    open(file) do f
        s = read(f, Int32)
        d = read(f, Int32)
        w = Vector{Float64}(undef, s)
        read!(f, w)
        p = Vector{SVector{3,Float64}}(undef, s)
        read!(f, p)
        delete = Vector{Int}(undef, s)
        j = 1
        for i in 1:s
            p[i] = _rotmatrix*p[i]
            if islinearsymmetric
                w[i] = 2*w[i]
                if p[i][1] < 0.0
                    delete[j] = i
                    j += 1
                end
            end
        end
        @assert !islinearsymmetric || 2*(j-1) == s
        resize!(delete, j-1)
        deleteat!(p, delete)
        deleteat!(w, delete)
        LebedevGrid(length(p), d, w, p)
    end
end

function lebedev_num(size::Integer, islinearsymmetric::Bool)
    idx = searchsortedfirst(lebedev_sizes, size*(1+islinearsymmetric))
    lebedev_sizes[idx-(idx>length(lebedev_sizes))]
end

function lebedev_num(mol, num::Integer)
    length(mol) == 1 || num ≤ 1 && return 1
    islin = is_zaxis_linear(mol)
    issym = is_zaxis_linear_symmetric(mol, islin)
    lebedev_num(islin ? num : div(num, 5), issym)
end

function get_lebedev(size, islinearsymmetric)
    kept = lebedev_num(size, islinearsymmetric)
    ret = read_lebedev_grid(joinpath(lebedev_path, string(kept)), islinearsymmetric)
    # ret.size != size && @info "Using a Lebedev grid of size $(ret.size)"
    ret
end

function get_lebedev_direct(size)
    idx = searchsortedfirst(lebedev_sizes, size)
    islinearsymmetric = false
    realsize = size
    if get(lebedev_sizes, idx, -1) != size
        realsize = 2*size
        idx = searchsortedfirst(lebedev_sizes, realsize)
        islinearsymmetric = true
    end
    if get(lebedev_sizes, idx, -1) != realsize
        size == 1 && return LebedevGrid(one(Int32), one(Int32), [4π], [zero(SVector{3,Float64})])
        error(lazy"Size $size does not correspond to a registered lebedev grid")
    end
    read_lebedev_grid(joinpath(lebedev_path, string(realsize)), islinearsymmetric)
end

"""
    get_rotation_matrices(mol::AbstractSystem, num, strict=false)

Return a tuple `(rots, weights)` where:
- `rots` is a list of 3×3 rotation matrices
- `weights` is the list of their corresponding weight

With these, it is possible to approximate ``∫f(θ,ϕ,ψ) dθ dϕ dψ``, where `f` is a function
of the orientation of system `mol` along the three Euler angles `θ`, `φ` and `ψ`, by:
```julia
sum(weight * f(ChangePositionSystem(mol, (rot,) .* position(mol))) for (rot, weight) in zip(rots, weights))
```

!!! warning
    In order to compute the average of `f`, the result should be divided by 4π since
    `sum(weights) == 4π`

`num` is the approximate precision of the lebedev grid. The actual number may be much lower
if the molecule has symmetries, which will keep the same precision.
Use `strict=true` to make sure that the lebedev points is exactly `num`.
"""
function get_rotation_matrices(mol::AbstractSystem, num, strict=false)
    if length(mol) == 1
        strict && num != 1 && error("Cannot return more than one rotation matrix for a monoatomic molecule")
        return [one(SMatrix{3,3,Float64,9})], [4π]
    end
    islin = is_zaxis_linear(mol)
    issym = is_zaxis_linear_symmetric(mol, islin)
    lebedev = strict ? get_lebedev_direct(num) : get_lebedev(islin ? num : div(num, 5), issym)
    zrots = [SMatrix{3,3,Float64,9}([cospi(2i/5) -sinpi(2i/5) 0; sinpi(2i/5) cospi(2i/5) 0; 0 0 1]) for i in 0:(4-4islin)]
    rots = SMatrix{3,3,Float64,9}[]
    for zrot in zrots, point in lebedev.points
        push!(rots, SMatrix{3,3,Float64,9}(hcat([1, 0, 0], [0, 1, 0], point))*zrot)
    end
    return rots, (islin ? lebedev.weights : repeat(lebedev.weights, 5)./5)
end

function lebedev_average(l, weights)
    tot = 0.0
    for (x, w) in zip(l, weights)
        tot += w*x
    end
    tot / (4π)
end
