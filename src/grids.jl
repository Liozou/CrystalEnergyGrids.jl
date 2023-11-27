using StaticArrays

"""
    EnergyGrid

Representation of an interpolable energy grid.

Use [`interpolate_grid`](@ref) to access the value of the grid at any point in space.
"""
struct EnergyGrid
    num_unitcell::NTuple{3,Cint}
    ewald_precision::Cfloat # Inf for VdW grid, -Inf for empty grid
    higherorder::Bool # true if grid contains derivatives, false if only raw values
    grid::Array{Cfloat,4}
end
function EnergyGrid()
    EnergyGrid((0, 0, 0), -Inf, false, Array{Cfloat,4}(undef, 0, 0, 0, 0))
end
