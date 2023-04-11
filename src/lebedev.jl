using Artifacts
using StaticArrays

const _rotmatrix = SMatrix{3,3,Float64,9}([-0.17963068200890037 -0.21953827352603253 -0.9589242746631385; -0.9599246581752935 0.25228151379218244 0.12206018362173197; 0.21512198564550156 0.9424208106021077 -0.2560577025515984])

const lebedev_path = artifact"lebedev"
const lebedev_sizes = sort!([parse(Int, i) for i in readdir(lebedev_path)])

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
