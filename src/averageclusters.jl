function pdb_properties(trajectory::AbstractString)
    open(trajectory) do io
        readline(io)
        l0 = readline(io)
        a, b, c, α, β, γ = parse.(Float64, l0[rge] for rge in (7:15, 16:24, 25:33, 34:40, 41:47, 48:54))::Vector{Float64}
        cosα = cosd(α); cosβ = cosd(β); sinγ, cosγ = sincosd(γ)
        ω = sqrt(1 - cosα^2 - cosβ^2 - cosγ^2 + 2*cosα*cosβ*cosγ)
        mat = SMatrix{3,3,Float64,9}(a, 0, 0, b*cosγ, b*sinγ, 0, c*cosβ, c*(cosα - cosβ*cosγ)/sinγ, c*ω/sinγ)
        l = readline(io)
        counter = 0
        while !startswith(l, "ENDMDL")
            l = readline(io)
            counter += 1
        end
        mat*u"Å", counter
    end
end

function bin_trajectory(trajectory, energies::Vector{TK}; refT::TK=300.0u"K", step=0.15u"Å")
    mat, numatoms = pdb_properties(trajectory)
    numatoms > 0 || error(lazy"Invalid or empty PDB file at $trajectory: could not find any atom")
    invmat = inv(NoUnits.(mat./u"Å"))
    na, nb, nc = floor.(Int, NoUnits.(cell_lengths(mat) ./ step))
    bins = zeros(Float64, na, nb, nc)
    mine::TK = minimum(energies)
    lines = eachline(trajectory)
    state = iterate(lines)
    modelcounter = 1
    while !isnothing(state) && startswith(state[1], "MODEL")
        energy = energies[modelcounter]
        iterate(lines)
        for _ in 1:numatoms
            l = first(iterate(lines))
            pos = SVector{3,Float64}(parse(Float64, @view l[31:38]), parse(Float64, @view l[39:46]), parse(Float64, @view l[47:54]))
            y = invmat*pos
            i, j, k = y .- floor.(y)
            val = exp((mine-energy)/refT)
            bins[ceil(Int, i*na)+(i==0), ceil(Int, j*nb)+(j==0), ceil(Int, k*nc)+(k==0)] += val
        end
        lend = first(iterate(lines))
        @assert startswith(lend, "ENDMDL")
        modelcounter += 1
        state = iterate(lines)
    end
    λ = sum(exp((mine-e)/refT) for e in energies; init=0.0)
    bins ./= λ
    bins
end

function _add_cyclic((x, y, z), (i, j, k), (a, b, c), λ=true)
    ma, mb, mc = a÷2, b÷2, c÷2
    u = x ≤ ma && i > ma ? (i-a) : x > ma && i ≤ ma ? (i+a) : i
    v = y ≤ mb && j > mb ? (j-b) : y > mb && j ≤ mb ? (j+b) : j
    w = z ≤ mc && k > mc ? (k-c) : z > mc && k ≤ mc ? (k+c) : k
    SVector{3,Float64}(λ*u, λ*v, λ*w)
end

function average_clusters(bins::Array{Float64,3}; merge_manhattan=3, threshold=1e-5)
    basinsets = local_basins(bins, -Inf; lt=(>), merge_manhattan, skip=iszero)
    n = length(basinsets)
    @assert !any(isempty, basinsets)
    points, nodes = decompose_basins(bins, basinsets)
    probabilities = zeros(n)
    center = [zero(SVector{3,Float64}) for _ in 1:n]
    numnodes = zeros(Float64, n)
    for (i,j,k) in points
        belongs = nodes[i,j,k]
        @assert !isempty(belongs)
        λ = inv(length(belongs))
        for b in belongs
            probabilities[b] += bins[i,j,k]*λ
            numnodes[b] += λ
            center[b] += _add_cyclic(center[b], SVector{3,Float64}(i,j,k), size(bins), λ)
        end
    end
    ret = [(p, c/λ) for (p, c, λ) in zip(probabilities, center, numnodes)]
    filter!((>(threshold))∘first, ret)
    ret
end

function average_clusters(trajectory, energies::Vector{TK}; merge_manhattan=3,
                          refT::TK=300.0u"K", step=0.15u"Å", threshold=1e-5)
    average_clusters(bin_trajectory(trajectory, energies; refT, step); merge_manhattan, threshold)
end
