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
    spacing::typeof(1.0u"Å")
    dims::SVector{3,Cint}
    size::SVector{3,typeof(1.0u"Å")}
    shift::SVector{3,typeof(1.0u"Å")}
    unitcell::SVector{3,typeof(1.0u"Å")}
    Δ::SVector{3,typeof(1.0u"Å")}
end
GridCoordinatesSetup() = GridCoordinatesSetup(CellMatrix(),
                                              0.0u"Å",
                                              SVector{3}((0, 0, 0)),
                                              SVector{3}((0.0u"Å", 0.0u"Å", 0.0u"Å")),
                                              SVector{3}((0.0u"Å", 0.0u"Å", 0.0u"Å")),
                                              SVector{3}((0.0u"Å", 0.0u"Å", 0.0u"Å")),
                                              SVector{3}((0.0u"Å", 0.0u"Å", 0.0u"Å")))

function GridCoordinatesSetup(framework::AbstractSystem{3}, spacing::typeof(1.0u"Å"))
    a, b, c = bounding_box(framework)
    cell = CellMatrix(framework)
    size = SVector{3,typeof(1.0u"Å")}(abs.(a)) + SVector{3,typeof(1.0u"Å")}(abs.(b)) + SVector{3,typeof(1.0u"Å")}(abs.(c))
    shift = SVector{3,typeof(1.0u"Å")}(min.(a, 0.0u"Å")) + SVector{3,typeof(1.0u"Å")}(min.(b, 0.0u"Å")) + SVector{3,typeof(1.0u"Å")}(min.(c, 0.0u"Å"))
    _dims = floor.(Cint, (size ./ spacing))
    dims = _dims + iseven.(_dims)
    unitcell = SVector{3,typeof(1.0u"Å")}((norm(a), norm(b), norm(c)))
    Δ = size ./ dims
    GridCoordinatesSetup(cell, spacing, dims, size, shift, unitcell, Δ)
end

function wrap_atom(point, cell::CellMatrix)
    abc = cell.invmat * (point isa SVector ? point : SVector{3}(point))
    cell.mat * (abc .- floor.(abc))
end

function offsetpoint(point, csetup::GridCoordinatesSetup)
    newpoint = wrap_atom(point, csetup.cell)
    @. (newpoint - csetup.shift)*csetup.dims/csetup.size + 1
end

function inverse_offsetpoint(ipoint, cset::GridCoordinatesSetup)
    @. (ipoint - 1)*cset.size/cset.dims + cset.shift
end

function abc_to_xyz(cset::GridCoordinatesSetup, i, j, k)
    SVector{3,typeof(1.0u"Å")}((i*cset.size[1]/cset.dims[1] + cset.shift[1],
                                j*cset.size[2]/cset.dims[2] + cset.shift[2],
                                k*cset.size[3]/cset.dims[3] + cset.shift[3]))
end

"""
    BlockFile

Structure representing the content of a .block file projected on a grid.
"""
struct BlockFile
    csetup::GridCoordinatesSetup
    block::BitArray{3}
    empty::Bool
end
BlockFile() = BlockFile(GridCoordinatesSetup(), BitArray{3}(undef, 0, 0, 0), true)
function BlockFile(csetup::GridCoordinatesSetup, block::BitArray{3})
    if any(block)
        BlockFile(csetup, block, false)
    else
        BlockFile()
    end
end
Base.isempty(x::BlockFile) = x.empty
function Base.getindex(x::BlockFile, pos::SVector{3,typeof(1.0u"Å")})
    x.empty && return false
    a, b, c = round.(Int, offsetpoint(pos, x.csetup))
    x.block[a,b,c]
end
Base.getindex(x::BlockFile, i, j, k) = x[SVector{3}(i, j, k)]

"""
    parse_blockfile(file, csetup)

Parse the .block input `file` into a grid of dimensions specified by `csetup::GridCoordinatesSetup`.

Return a [`BlockFile`](@ref)
"""
function parse_blockfile(file, csetup)
    lines = readlines(file)
    num = parse(Int, popfirst!(lines))
    num == 0 && return BlockFile()
    a, b, c = csetup.dims .+ 1
    # δmax = Int.(cld.(csetup.dims, 2))
    block = falses(a, b, c)
    mat = NoUnits.(csetup.cell.mat./u"Å")
    invmat = NoUnits.(csetup.cell.invmat.*u"Å")
    buffer, ortho, safemin = prepare_periodic_distance_computations(mat)
    safemin2 = safemin^2
    buffer2 = MVector{3,Float64}(undef)
    for (idx, l) in enumerate(lines)
        idx > num && ((@assert all(isspace, l)); continue)
        sl = split(l)
        radius = parse(Float64, pop!(sl))
        _center = csetup.cell.mat*SVector{3,Float64}(parse.(Float64, sl))
        ofscenter = offsetpoint(_center, csetup)
        ic, jc, kc = round.(Int, ofscenter)
        center = inverse_offsetpoint(ofscenter, csetup)
        # δi, δj, δk = min.(1 .+ ceil.(Int, radius*u"Å" ./ csetup.Δ), δmax)
        radius2 = radius^2
        for i in 1:a, j in 1:b, k in 1:c
        # for _i in ic-δi:ic+δi, _j in jc-δj:jc+δj, _k in kc-δk:kc+δk
        #     i = mod1(_i, a); j = mod1(_j, b); k = mod1(_k, c)
        # Using this approach is incorrect because a line in the block makes a line in the
        # crystal aligned with the cartesian referential, not the fractional one.
        # In other words, if the blocking sphere straddles a periodic boundary, the two
        # neighbour points in the crystal may not be neighbours in the block.
            block[i,j,k] && continue
            buffer .= NoUnits.((center .- inverse_offsetpoint(SVector{3,Int}(i,j,k), csetup)) ./ u"Å")
            if periodic_distance2_fromcartesian!(buffer, mat, invmat, ortho, safemin2, buffer2) < radius2
                block[i,j,k] = true
            end
        end
        if !block[ic,jc,kc]
            @warn "Blocking sphere of $radius Å around position $center is too small to be used with grid step $(csetup.spacing)"
        end
    end
    return BlockFile(csetup, block)
end
