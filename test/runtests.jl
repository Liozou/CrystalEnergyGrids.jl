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
    setupArCHA = setup_RASPA("CHA_1.4_3b4eeb96_Na_11812", "BoulfelfelSholl2021", "Ar", "TraPPE", 0.15, (1,1,1), nothing);
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
    @test egridArCHA[only(minimaAr)] ≈ -1854.646947681577


    setupNaCHA = setup_RASPA("CHA_1.4_3b4eeb96", "BoulfelfelSholl2021", "Na", "TraPPE", 0.15, (1,1,1), nothing);
    egridNaCHA = dropdims(energy_grid(setupNaCHA, 0.3u"Å"); dims=1)
    serialize(joinpath(TESTDIR, "savegrids", "Na_CHA_1.4_3b4eeb96"), egridNaCHA)
    vdw, coulomb = energy_point(setupNaCHA, [SVector{3}([0.0, 0.0, 0.0])*u"Å"])
    @test vdw ≈ -11083.130753181682u"K"
    @test coulomb ≈ -1.850940229601566e6u"K"
    @test vdw + coulomb == egridNaCHA[1,1,1]u"K"

    minimaNa = sort!(CEG.local_minima(egridNaCHA))
    @test minimaNa == [CartesianIndex(72, 64, 87)]
    @test egridNaCHA[only(minimaNa)] ≈ -4.629398710521718e13
    @test egridNaCHA[only(minimaNa)] == minimum(egridNaCHA)

    ewald = setupNaCHA.ewald
    co2 = CEG.load_molecule_RASPA("CO2", "TraPPE", "BoulfelfelSholl2021");
    host, adsorbate = compute_ewald(ewald, [CEG.ChangePositionSystem(co2, [[1.3, 2.9, 1.149]u"Å", [1.3, 2.9, 0.0]u"Å", [1.3, 2.9, -1.149]u"Å"]), CEG.ChangePositionSystem(co2, [[4.7, 10.1, 1.149]u"Å", [4.7, 10.1, 0.0]u"Å", [4.7, 10.1, -1.149]u"Å"])])
    @test host ≈ 1177.521548912203u"K"
    @test adsorbate ≈ 7.4133039820109055u"K"

    # ar = CEG.load_molecule_RASPA("Ar", "TraPPE", "BoulfelfelSholl2021")
    # boulfelfelsholl2021 = CEG.parse_forcefield_RASPA("BoulfelfelSholl2021")
end

rm(GRIDDIR; recursive=true)
