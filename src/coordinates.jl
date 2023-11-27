"""
    GridCoordinatesSetup

Structure representing the coordinates of points in a 3D regular grid.
The following attributes are part of the API:
- `cell`: the [`CellMatrix`] of the unit cell on which the grid is mapped.
- `spacing`: the required spacing between consecutive points. The actual spacing is given
  by `Δ`.
- `dims`: the number of points along each dimension minus one.
- `size`: furthest point from the origin (after shifting).
- `shift`: shift vector to give all points in the unit cell positives coordinates.
- `unitcell`: length of each axis of the unit cell.
- `Δ`: spacing between two consecutive grid points, along each dimension.
"""
struct GridCoordinatesSetup
    cell::CellMatrix
    spacing::TÅ
    dims::SVector{3,Cint}
    size::SVector{3,TÅ}
    shift::SVector{3,TÅ}
    unitcell::SVector{3,TÅ}
    Δ::SVector{3,TÅ}
end
GridCoordinatesSetup() = GridCoordinatesSetup(CellMatrix(),
                                              0.0u"Å",
                                              SVector{3}((0, 0, 0)),
                                              SVector{3}((0.0u"Å", 0.0u"Å", 0.0u"Å")),
                                              SVector{3}((0.0u"Å", 0.0u"Å", 0.0u"Å")),
                                              SVector{3}((0.0u"Å", 0.0u"Å", 0.0u"Å")),
                                              SVector{3}((0.0u"Å", 0.0u"Å", 0.0u"Å")))

function GridCoordinatesSetup(cell::CellMatrix, spacing::TÅ)
    a, b, c = eachcol(cell.mat)
    size = SVector{3,TÅ}(abs.(a)) + SVector{3,TÅ}(abs.(b)) + SVector{3,TÅ}(abs.(c))
    shift = SVector{3,TÅ}(min.(a, 0.0u"Å")) + SVector{3,TÅ}(min.(b, 0.0u"Å")) + SVector{3,TÅ}(min.(c, 0.0u"Å"))
    _dims = floor.(Cint, (size ./ spacing))
    dims = _dims + iseven.(_dims)
    unitcell = SVector{3,TÅ}((norm(a), norm(b), norm(c)))
    Δ = size ./ dims
    GridCoordinatesSetup(cell, spacing, dims, size, shift, unitcell, Δ)
end

function GridCoordinatesSetup(framework::AbstractSystem{3}, spacing::TÅ)
    GridCoordinatesSetup(CellMatrix(framework), spacing)
end


function wrap_atom(point, cell::CellMatrix)
    abc = cell.invmat * (point isa SVector ? point : SVector{3}(point))
    cell.mat * (abc .- floor.(abc))
end

function offsetpoint(point, csetup::GridCoordinatesSetup)
    newpoint = wrap_atom(point, csetup.cell)
    @. (newpoint - csetup.shift)*csetup.dims/csetup.size + 1
end

function inverse_offsetpoint(ipoint, csetup::GridCoordinatesSetup)
    @. (ipoint - 1)*csetup.Δ + csetup.shift
end

function abc_to_xyz(cset::GridCoordinatesSetup, i, j, k)
    SVector{3,TÅ}((i*cset.size[1]/cset.dims[1] + cset.shift[1],
                                j*cset.size[2]/cset.dims[2] + cset.shift[2],
                                k*cset.size[3]/cset.dims[3] + cset.shift[3]))
end
