using Test
using CrystalEnergyGrids
import CrystalEnergyGrids as CEG

using StaticArrays
using Unitful, UnitfulAtomic
using AtomsBase

using Serialization

# using Aqua
# Aqua.test_all(CrystalEnergyGrids; ambiguities=false)

const TESTDIR = joinpath(dirname(dirname(pathof(CEG))), "test")
setdir_RASPA!(joinpath(TESTDIR, "raspa"))

const GRIDDIR = joinpath(TESTDIR, "savegrids")
if isdir(GRIDDIR)
    rm(GRIDDIR; recursive=true)
end
mkdir(GRIDDIR)

@testset "CrystalEnergyGrids" begin
    setupArCHA = setup_RASPA("CHA_1.4_3b4eeb96_Na_11812", "BoulfelfelSholl2021", "Ar", "TraPPE"; blockfile=nothing);
    egridArCHA = dropdims(energy_grid(setupArCHA, 0.3u"Å"); dims=1)
    serialize(joinpath(TESTDIR, "savegrids", "Ar_CHA_1.4_3b4eeb96_Na_11812"), egridArCHA)

    vdw, coulomb = energy_point(setupArCHA, [SVector{3}([0.0, 0.0, 0.0])*u"Å"])
    @test iszero(coulomb)
    @test !iszero(vdw) && vdw == egridArCHA[1,1,1]u"K"
    # vdw, coulomb = energy_point(setupArCHA, [SVector{3}([0.2, 0.0, 0.0])*u"Å"])
    # @test iszero(coulomb)
    # @test !iszero(vdw) && vdw == egridArCHA[2,1,1]*u"K"

    minimaAr = sort!(CEG.local_minima(egridArCHA))
    @test minimaAr == [CartesianIndex(51, 50, 91)]
    @test egridArCHA[only(minimaAr)] == minimum(egridArCHA)
    @test egridArCHA[only(minimaAr)] ≈ -1854.9285197062766


    setupNaCHA = setup_RASPA("CHA_1.4_3b4eeb96", "BoulfelfelSholl2021", "Na", "TraPPE");
    egridNaCHA = dropdims(energy_grid(setupNaCHA, 0.3u"Å"); dims=1)
    serialize(joinpath(TESTDIR, "savegrids", "Na_CHA_1.4_3b4eeb96"), egridNaCHA)
    vdw, coulomb = energy_point(setupNaCHA, [SVector{3}([0.0, 0.0, 0.0])*u"Å"])
    @test vdw ≈ -11083.13758653269u"K"
    @test coulomb ≈ -1.8645348919601506e6u"K"
    @test vdw + coulomb == egridNaCHA[1,1,1]u"K"

    minimaNa = sort!(CEG.local_minima(egridNaCHA, 0.0))
    @test minimaNa == [CartesianIndex(29, 60, 60)]
    @test egridNaCHA[only(minimaNa)] ≈ -1.941489105927073e6
    @test egridNaCHA[only(minimaNa)] == minimum(egridNaCHA)

    ewald = setupNaCHA.ewald
    co2 = CEG.load_molecule_RASPA("CO2", "TraPPE", "BoulfelfelSholl2021");
    reciprocal = compute_ewald(ewald, [CEG.ChangePositionSystem(co2, [[1.3, 2.9, 1.149]u"Å", [1.3, 2.9, 0.0]u"Å", [1.3, 2.9, -1.149]u"Å"]), CEG.ChangePositionSystem(co2, [[4.7, 10.1, 1.149]u"Å", [4.7, 10.1, 0.0]u"Å", [4.7, 10.1, -1.149]u"Å"])])
    @test reciprocal ≈ 1177.5215489122043u"K" + 7.4133039820109055u"K"

    # ar = CEG.load_molecule_RASPA("Ar", "TraPPE", "BoulfelfelSholl2021")
    # boulfelfelsholl2021 = CEG.parse_forcefield_RASPA("BoulfelfelSholl2021")
end


@testset "IncrementalEwaldContext" begin
    setupNaCHA = setup_RASPA("CHA_1.4_3b4eeb96", "BoulfelfelSholl2021", "Na", "TraPPE");
    co2 = CEG.load_molecule_RASPA("CO2", "TraPPE", "BoulfelfelSholl2021");
    pos1 = [SVector{3}((1.0, 2.5, 1.7))*u"Å"]
    ctx = CEG.EwaldContext(setupNaCHA.ewald,[
        CEG.ChangePositionSystem(setupNaCHA.molecule, pos1),
        CEG.ChangePositionSystem(setupNaCHA.molecule, [SVector{3}((6.2, 5.1, 3.0))*u"Å"]),
        co2
    ])
    iec = CEG.IncrementalEwaldContext(ctx)
    reciprocal = compute_ewald(ctx)
    @test reciprocal == compute_ewald(iec)
    # /!\ call to single_contribution_ewald must be after call to compute_ewald
    @test reciprocal ≈ compute_ewald(iec, 1) + CEG.single_contribution_ewald(iec, 1, pos1)

end


rm(GRIDDIR; recursive=true)
