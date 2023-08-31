using Test
using CrystalEnergyGrids
import CrystalEnergyGrids as CEG

using StaticArrays
using Unitful, UnitfulAtomic
using AtomsBase

using Serialization

# using Aqua
# Aqua.test_all(CrystalEnergyGrids; ambiguities=false)

const TESTDIR = joinpath(dirname(dirname(pathof(CEG))), "test"); setdir_RASPA!(joinpath(TESTDIR, "raspa"))

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
    @test egridArCHA[only(minimaAr)] ≈ -1854.9285197062766 rtol=0.001

    setupNaCHA = setup_RASPA("CHA_1.4_3b4eeb96", "BoulfelfelSholl2021", "Na", "TraPPE");
    egridNaCHA = dropdims(energy_grid(setupNaCHA, 0.3u"Å"); dims=1)
    serialize(joinpath(TESTDIR, "savegrids", "Na_CHA_1.4_3b4eeb96"), egridNaCHA)
    vdw, coulomb = energy_point(setupNaCHA, [SVector{3}([0.0, 0.0, 0.0])*u"Å"])
    @test vdw ≈ -11083.13758653269u"K" rtol=0.001
    @test coulomb ≈ -1.8509402225092095e6u"K" rtol=0.001
    @test vdw + coulomb == egridNaCHA[1,1,1]u"K"

    minimaNa = sort!(CEG.local_minima(egridNaCHA, 0.0))
    @test minimaNa == [CartesianIndex(29, 60, 60)]
    @test egridNaCHA[only(minimaNa)] ≈ -1.9278944364761321e6 rtol=0.001
    @test egridNaCHA[only(minimaNa)] == minimum(egridNaCHA)

    ewald = setupNaCHA.ewald
    co2 = CEG.load_molecule_RASPA("CO2", "TraPPE", "BoulfelfelSholl2021");
    reciprocal = compute_ewald(ewald, [CEG.ChangePositionSystem(co2, [[1.3, 2.9, 1.149]u"Å", [1.3, 2.9, 0.0]u"Å", [1.3, 2.9, -1.149]u"Å"]), CEG.ChangePositionSystem(co2, [[4.7, 10.1, 1.149]u"Å", [4.7, 10.1, 0.0]u"Å", [4.7, 10.1, -1.149]u"Å"])])
    @test reciprocal ≈ 1177.5215489122043u"K" + 7.4133039820109055u"K" rtol=0.001

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

@testset "MonteCarloSimulation" begin
    co2 = CEG.load_molecule_RASPA("CO2", "TraPPE", "BoulfelfelSholl2021");
    na = CEG.load_molecule_RASPA("Na", "TraPPE", "BoulfelfelSholl2021");
    ar = CEG.load_molecule_RASPA("Ar", "TraPPE", "BoulfelfelSholl2021");
    pos1 = [SVector{3}([0.0, 2.0, 4.0]u"Å")]
    pos2 = [SVector{3}([1.0, 3.0, 1.0]u"Å")]

    mcAr1, _ = setup_montecarlo("CHA_1.4_3b4eeb96", "BoulfelfelSholl2021",
                                [co2, CEG.ChangePositionSystem(ar, pos1)]);
    mcAr2, _ = setup_montecarlo("CHA_1.4_3b4eeb96", "BoulfelfelSholl2021",
                                [co2, CEG.ChangePositionSystem(ar, pos2)]);
    baseAr1 = CEG.baseline_energy(mcAr1)
    baseAr2 = CEG.baseline_energy(mcAr2)
    mov1 = CEG.movement_energy(mcAr1, (2,1))
    @test baseAr1.er.reciprocal == baseAr2.er.reciprocal
    @test baseAr1.er.framework.direct == baseAr2.er.framework.direct
    @test baseAr1.er.inter != baseAr2.er.inter
    @test baseAr1.er.framework.vdw != baseAr2.er.framework.vdw
    @test mov1.reciprocal == mov1.framework.direct == 0.0u"K"
    @test mov1 == CEG.movement_energy(mcAr1, (2,1), pos1)
    @test Float64(baseAr2) ≈ Float64(baseAr1 - mov1 + CEG.movement_energy(mcAr1, (2,1), pos2))

    mcNa1, _ = setup_montecarlo("CHA_1.4_3b4eeb96", "BoulfelfelSholl2021",
                                [co2, CEG.ChangePositionSystem(na, pos1)]);
    mcNa2, _ = setup_montecarlo("CHA_1.4_3b4eeb96", "BoulfelfelSholl2021",
                                [co2, CEG.ChangePositionSystem(na, pos2)]);
    baseNa1 = CEG.baseline_energy(mcNa1)
    movNa1_1 = CEG.movement_energy(mcNa1, (2,1))
    movNa1_2 = CEG.movement_energy(mcNa1, (2,1), pos2)
    @test movNa1_1 == CEG.movement_energy(mcNa1, (2,1), pos1)
    @test Float64(CEG.baseline_energy(mcNa2)) ≈ Float64(baseNa1 - movNa1_1 + movNa1_2)
end

@testset "Triclinic" begin
    # triclinic input
    setupArCIT7 = setup_RASPA("CIT7", "BoulfelfelSholl2021", "Ar", "TraPPE"; blockfile=false);
    @test only(setupArCIT7.grids).num_unitcell == (2, 3, 3)
    vdw, coulomb = energy_point(setupArCIT7, [SVector{3}([-6.2437738548095165, 22.4072046579578092, 5.7557837121224642]u"Å")])
    @test iszero(coulomb)
    @test vdw ≈ -1043.35893781u"K" rtol=0.001
    cit7 = setupArCIT7.framework
    ΠA, ΠB, ΠC = setupArCIT7.grids[1].num_unitcell
    Π = ΠA*ΠB*ΠC
    n = length(cit7)
    _positions = Vector{SVector{3,typeof(1.0u"Å")}}(undef, n*Π)
    _symbols = Vector{Symbol}(undef, n*Π)
    mat = setupArCIT7.grids[1].mat
    axeA, axeB, axeC = eachcol(mat)
    iofs = 0
    for πa in 1:ΠA, πb in 1:ΠB, πc in 1:ΠC
        posofs = ((πa-1) .* axeA .+ (πb-1) .* axeB .+ (πc-1) .* axeC)*u"Å"
        for (i, (at, pos)) in enumerate(zip(cit7.atomic_symbol, cit7.position))
            _positions[iofs+i] = pos .+ posofs
            _symbols[iofs+i] = Symbol(split(String(at), '_')[1])
        end
        iofs += n
    end
    _numbers = repeat(cit7.atomic_number, Π)
    _masses = repeat(cit7.atomic_mass, Π)
    _charges = repeat(cit7.atomic_charge, Π)
    superCIT7 = RASPASystem(AtomsBase.bounding_box(setupArCIT7.framework) .* (ΠA, ΠB, ΠC), _positions, _symbols, _numbers, _masses, _charges, false)
    supermat = [axeA .* ΠA;; axeB .* ΠB;; axeC .* ΠC]
    ar = CEG.load_molecule_RASPA("Ar", "TraPPE", "BoulfelfelSholl2021");
    molAr = CEG.ChangePositionSystem(ar, [SVector{3}([5.86207, 16.23616, 17.52779]u"Å")]);
    mcArCIT7, _ = setup_montecarlo("CIT7", "BoulfelfelSholl2021", [molAr]);
    mcVoidArCIT7, _ = setup_montecarlo(supermat, "BoulfelfelSholl2021", [molAr, superCIT7]);
    CEG.baseline_energy(mcArCIT7); movArCIT7 = CEG.movement_energy(mcArCIT7, (1,1))
    CEG.baseline_energy(mcVoidArCIT7); movVoidArCIT7 = CEG.movement_energy(mcVoidArCIT7, (1,1))
    @test Float64(movArCIT7) ≈ Float64(movVoidArCIT7) rtol=0.005

    molAr1 = CEG.ChangePositionSystem(ar, [SVector{3}([-7.7365250811304911,31.5070011601372251,1.5285305931479920]u"Å")]);
    molAr2 = CEG.ChangePositionSystem(ar, [SVector{3}([10.7586599791867421,-2.3259182727570948,20.5642722996513001]u"Å")]);
    mcAr2, _ = setup_montecarlo("CIT7", "BoulfelfelSholl2021", [molAr1, molAr2]);
    baseAr2 = Float64(CEG.baseline_energy(mcAr2))
    movAr2_1 = Float64(CEG.movement_energy(mcAr2, (1,1)))
    CEG.baseline_energy(mcAr2); movAr2_2 = Float64(CEG.movement_energy(mcAr2, (1,2)))
    @test baseAr2 ≈ -1789.77383582 rtol=0.001
    @test baseAr2 - movAr2_1 - movAr2_2 ≈ 0.28797384 rtol=0.001


    na = CEG.load_molecule_RASPA("Na", "TraPPE", "BoulfelfelSholl2021");
    molNaSm1 = CEG.ChangePositionSystem(na, [SVector{3}([4.935357501688667, 23.53557287745349, 25.71480449842175]u"Å")]);
    mcNaSm1, _ = setup_montecarlo("SuperCIT7m1", "BoulfelfelSholl2021", [molNaSm1]);
    @test Float64(CEG.baseline_energy(mcNaSm1)) ≈ -8996.975999017683 rtol=0.001

    # test one atom in supercell
    molNaMini = CEG.ChangePositionSystem(na, [SVector{3}([-1.401612509676063, 14.86235802394228, 15.37932058231622]u"Å")]);
    mcNaMini, _ = setup_montecarlo("Mini", "BoulfelfelSholl2021", [molNaMini]);
    mcNaMiniRef, _ = setup_montecarlo("MiniRef", "BoulfelfelSholl2021", [molNaMini]);
    baseMini = Float64(CEG.baseline_energy(mcNaMini))
    @test baseMini ≈ Float64(CEG.baseline_energy(mcNaMiniRef))
    @test baseMini ≈ -248304.58180794 rtol=0.001

    molNaPetit = CEG.ChangePositionSystem(na, [SVector{3}([18.77838182689036, 14.73031622108175, 3.669308624409278]u"Å")]);
    mcNaPetit, _ = setup_montecarlo("Petit", "BoulfelfelSholl2021", [molNaPetit]);
    mcNaPetitRef, _ = setup_montecarlo("PetitRef", "BoulfelfelSholl2021", [molNaPetit]);
    basePetit = Float64(CEG.baseline_energy(mcNaPetit))
    @test basePetit ≈ Float64(CEG.baseline_energy(mcNaPetitRef)) rtol=0.0001
    @test basePetit ≈ 262204.85720076 rtol=0.001

    molNaSolo = CEG.ChangePositionSystem(na, [SVector{3}([-5.485237153390736, 18.66002404042065, 24.29891456309461]u"Å")]);
    mcNaSolo, _ = setup_montecarlo("CIT7", "BoulfelfelSholl2021", [molNaSolo]);
    @test mcNaSolo.tailcorrection[] ≈ -70.44772635984882u"K" # account for supercell


    # molNa1 = CEG.ChangePositionSystem(na, [SVector{3}([0.9813979229493714, 10.43406962195019, 23.05994195152032]u"Å")]);
    # molNa2 = CEG.ChangePositionSystem(na, [SVector{3}([1.369507482331668, 19.16986027450343, 7.145326794850396]u"Å")]);
    # mcNa2NoVDW, _ = setup_montecarlo("CIT7", "BoulfelfelSholl2021", [molNa1, molNa2]);

    # co2 = CEG.load_molecule_RASPA("CO2", "TraPPE", "BoulfelfelSholl2021");
    # molNa0 = CEG.ChangePositionSystem(na, [SVector{3}([-6.901, 14.575, 22.728]u"Å")]);
    # molCO21 = CEG.ChangePositionSystem(co2, [SVector{3}([-7.464, 8.180, 20.644]u"Å"),
    #                                          SVector{3}([-8.467, 8.482, 21.117]u"Å"),
    #                                          SVector{3}([-9.469, 8.784, 21.589]u"Å"),]);
    # molCO22 = CEG.ChangePositionSystem(co2, [SVector{3}([12.666, 18.725, 24.722]u"Å"),
    #                                          SVector{3}([12.058, 17.758, 24.847]u"Å"),
    #                                          SVector{3}([11.450, 16.792, 24.972]u"Å"),]);
    # mcNaCO2, _ = setup_montecarlo("CIT7", "BoulfelfelSholl2021", [molNa0, molCO21, molCO22]);

end

rm(GRIDDIR; recursive=true)
rm(joinpath(TESTDIR, "raspa", "grids"); recursive=true)
