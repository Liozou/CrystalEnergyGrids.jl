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

function get_lebedev(size, islinearsymmetric)
    idx = searchsortedfirst(lebedev_sizes, size*(1+islinearsymmetric))
    kept = lebedev_sizes[idx-(idx>length(lebedev_sizes))]
    ret = read_lebedev_grid(joinpath(lebedev_path, string(kept)), islinearsymmetric)
    # ret.size != size && @info "Using a Lebedev grid of size $(ret.size)"
    ret
end

"""
    get_rotation_matrices(mol::AbstractSystem, num)

Return a tuple `(rots, weights)` where:
- `rots` is a list of 3×3 rotation matrices
- `weights` is the list of their corresponding weight

With these, it is possible to approximate ``∫f(θ,ϕ,ψ) dθ dϕ dψ``, where `f` is a function
of the orientation of system `mol` along the three Euler angles `θ`, `φ` and `ψ`, by:
```julia
sum(weight * f(ChangePositionSystem(mol, (rot,) .* position(mol))) for (rot, weight) in zip(rots, weights))
```

!!! note
    In order to compute the average of `f`, the result should be divided by 4π since
    `sum(weights) == 4π`
"""
function get_rotation_matrices(mol::AbstractSystem, num)
    islin = is_zaxis_linear(mol)
    issym = is_zaxis_linear_symmetric(mol, islin)
    lebedev = get_lebedev(islin ? num : div(num, 5), issym)
    zrots = [SMatrix{3,3,Float64,9}([cospi(2i/5) -sinpi(2i/5) 0; sinpi(2i/5) cospi(2i/5) 0; 0 0 1]) for i in 0:(4-4islin)]
    rots = SMatrix{3,3,Float64,9}[]
    for zrot in zrots, point in lebedev.points
        push!(rots, SMatrix{3,3,Float64,9}(hcat([1, 0, 0], [0, 1, 0], point))*zrot)
    end
    rots, (islin ? lebedev.weights : repeat(lebedev.weights, 5)./5)
end
