const _TGERG = typeof(Clapeyron.GERG2008(["water"]))
const _TPR = typeof(Clapeyron.PR(["water"]))

struct GCMCData
    volumes::Vector{typeof(1.0u"Å^3")} # volumes[i] is the volume accessible to kind  i,
    # i.e. the volume of the framework minus that of the blocking spheres
    model::Vector{_TGERG}
    model0::Vector{_TPR}
    Π::Int # number of unit cells in the supercell
end

function Base.:(==)(g1::GCMCData, g2::GCMCData)
    (g1.volumes == g2.volumes && g1.Π == g2.Π) || return false
    n = length(g1.model)
    (n == length(g2.model) == length(g1.model0) == length(g2.model0)) || return false
    for i in 1:n, P in (100u"Pa", 1e6u"Pa"), T in (300u"K", 500u"K")
        φ1 = only(Clapeyron.fugacity_coefficient(g1.model[i], P, T; phase=:stable, vol0=Clapeyron.volume(g1.model0[i], P, T)))
        φ2 = only(Clapeyron.fugacity_coefficient(g2.model[i], P, T; phase=:stable, vol0=Clapeyron.volume(g2.model0[i], P, T)))
        φ1 ≈ φ2 || return false
    end
    return true
end

function GCMCData(ff::ForceField, ffidx, speciesblocks::Vector{BlockFile}, unitcell)
    n = length(ffidx)
    Π = prod(find_supercell(unitcell, ff.cutoff))
    cellvolume = det(unitcell)*Π
    volumes = Vector{typeof(1.0u"Å^3")}(undef, n)
    model = Vector{_TGERG}(undef, n)
    model0 = Vector{_TPR}(undef, n)
    for i in 1:length(ffidx)
        gasname = identify_molecule([ff.symbols[k] for k in ffidx[i]])
        gaskey = get(GAS_NAMES, gasname, "")
        isempty(gaskey) && continue
        block = speciesblocks[i]
        accessible_fraction = block.empty ? 1.0 : (1.0 - count(block.block)/length(block.block))
        volumes[i] = accessible_fraction * cellvolume
        model[i] = Clapeyron.GERG2008([gaskey])
        model0[i] = Clapeyron.PR([gaskey])
    end
    GCMCData(volumes, model, model0, Π)
end

const GAS_NAMES = Dict{String,String}(
    "Ar"    => "argon",
    "CH4"   => "methane",
    "CH4O"  => "methanol",
    "CH3OH" => "methanol",
    "CO"    => "carbon monoxide",
    "CO2"   => "carbon dioxide",
    "H2"    => "hydrogen",
    "H4Si"  => "silane",
    "SiH4"  => "silane",
    "H2O"   => "water",
    "He"    => "helium", # missing for PCSAFT
    "Kr"    => "krypton",
    "N2"    => "nitrogen",
    "NO"    => "nitrogen monoxide", # missing for PCSAFT
    "NO2"   => "nitrogen dioxide", # missing for PCSAFT
    "Ne"    => "neon", # missing for PCSAFT
    "O2"    => "oxygen",
    "O2S"   => "sulfur dioxide",
    "SO2"   => "sulfur dioxide",
    "Xe"    => "xenon",
)



# internal struct used to convey information about an attempted swap move
struct SwapInformation
    φPV_div_k::TK
    i::Int
    isswap::Bool
    isinsertion::Bool
end

function compute_accept_move_swap(diff::TK, T, mc, swapinfo::SwapInformation)
    expterm = exp(-diff/T)
    prefix = swapinfo.φPV_div_k/T
    N = length(mc.step.posidx[swapinfo.i])
    rand(mc.rng) < if swapinfo.isinsertion # if insertion
        # @show prefix/(N+1) * expterm
        prefix/(N+1) * expterm
    else # deletion
        # @show N/prefix * expterm
        N/prefix * expterm
    end
end
