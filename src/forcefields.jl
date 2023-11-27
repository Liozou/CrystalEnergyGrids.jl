# Definition of a force field

export ForceField

struct ForceField
    interactions::Matrix{Vector{Float64}}
    sdict::IdDict{Symbol,Int}
    symbols::Vector{Symbol}
    cutoff::TÅ
    name::String
end
ForceField() = ForceField(Matrix{Vector{Float64}}(undef, 0, 0), IdDict{Symbol,Int}(), Symbol[], 12.0u"Å", "")
