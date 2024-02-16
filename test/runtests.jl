using Test
using CrystalEnergyGrids
import CrystalEnergyGrids as CEG

prev_warning::Bool = CEG.PRINT_CHARGE_WARNING[]
CEG.PRINT_CHARGE_WARNING[] = false

using StaticArrays
using Unitful, UnitfulAtomic
using AtomsBase
using StableRNGs: StableRNG

using LinearAlgebra: det
using Serialization

# using Aqua
# Aqua.test_all(CrystalEnergyGrids; ambiguities=false)

const TESTDIR = joinpath(dirname(dirname(pathof(CEG))), "test"); setdir_RASPA!(joinpath(TESTDIR, "raspa"))
mkpath(joinpath(TESTDIR, "savegrids"))
# rm(joinpath(TESTDIR, "raspa", "grids"); recursive=true)

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
    @test minimaAr == [CartesianIndex(56, 88, 49), CartesianIndex(52, 50, 91), CartesianIndex(51, 51, 91)]
    @test egridArCHA[minimaAr[2]] == minimum(egridArCHA)
    @test egridArCHA[minimaAr[2]] ≈ -1841.0165092850448 rtol=0.001

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
    reciprocal = compute_ewald(ewald, ([CEG.ChangePositionSystem(co2, [[1.3, 2.9, 1.149]u"Å", [1.3, 2.9, 0.0]u"Å", [1.3, 2.9, -1.149]u"Å"]), CEG.ChangePositionSystem(co2, [[4.7, 10.1, 1.149]u"Å", [4.7, 10.1, 0.0]u"Å", [4.7, 10.1, -1.149]u"Å"])],))
    @test reciprocal ≈ 1177.5215489122043u"K" + 7.4133039820109055u"K" rtol=0.001

    # ar = CEG.load_molecule_RASPA("Ar", "TraPPE", "BoulfelfelSholl2021")
    # boulfelfelsholl2021 = CEG.parse_forcefield_RASPA("BoulfelfelSholl2021")
end

@testset "IncrementalEwaldContext" begin
    setupNaCHA = setup_RASPA("CHA_1.4_3b4eeb96", "BoulfelfelSholl2021", "Na", "TraPPE");
    co2 = CEG.load_molecule_RASPA("CO2", "TraPPE", "BoulfelfelSholl2021");
    pos1 = [SVector{3}((1.0, 2.5, 1.7))*u"Å"]
    pos2 = [SVector{3}((6.2, 5.1, 3.0))*u"Å"]
    ctx = CEG.EwaldContext(setupNaCHA.ewald,[
        [
            CEG.ChangePositionSystem(setupNaCHA.molecule, pos1),
            CEG.ChangePositionSystem(setupNaCHA.molecule, pos2),
        ],
        [co2]
    ])
    iec = CEG.IncrementalEwaldContext(ctx)
    reciprocal = compute_ewald(ctx)
    @test reciprocal == compute_ewald(iec)
    # /!\ call to single_contribution_ewald must be after call to compute_ewald
    @test reciprocal ≈ compute_ewald(iec.ctx, 1) + CEG.single_contribution_ewald(iec, 1, pos1)

    oldctx = deepcopy(ctx)
    @test CEG.remove_one_system!(ctx, 3) == 3
    @test CEG.add_one_system!(ctx, 2, position(co2)) == 3
    @test ctx == oldctx

    @test CEG.remove_one_system!(ctx, 2) == 3
    @test CEG.add_one_system!(ctx, 1, pos2) == 3
    @test compute_ewald(ctx) ≈ reciprocal
    @test ctx != oldctx

    tmpctx = deepcopy(ctx)
    @test CEG.remove_one_system!(ctx, 2) == 3
    @test CEG.add_one_system!(ctx, 2, position(co2)) == 3
    @test ctx == oldctx

    @test CEG.remove_one_system!(ctx, 2) == 3
    @test CEG.remove_one_system!(ctx, 1) == 2
    @test CEG.add_one_system!(ctx, 1, pos2) == 2
    @test CEG.add_one_system!(ctx, 1, pos1) == 3
    @test compute_ewald(ctx) ≈ reciprocal
    @test ctx != oldctx && ctx != tmpctx

    pos3 = [SVector{3}((4.1, 3.7, 2.2))*u"Å"]
    CEG.move_one_system!(ctx, 3, pos3)
    reciprocal3 = compute_ewald(ctx)
    @test CEG.remove_one_system!(ctx, 3) == 3
    @test CEG.add_one_system!(ctx, 1, pos3) == 3
    @test compute_ewald(ctx) ≈ reciprocal3

    co2_2 = CEG.ChangePositionSystem(co2, SVector{3}.([[5.491645446274333, 8.057854365959964, 8.669190836544463]u"Å",
                                                       [6.335120278303245, 7.462084936052019, 9.172986424179925]u"Å",
                                                       [7.178595110332157, 6.866315506144074, 9.676782011815387]u"Å"]));
    @test CEG.add_one_system!(ctx, 2, position(co2_2)) == 4
    ctx_twoco2 = CEG.EwaldContext(setupNaCHA.ewald,[
        [co2_2, co2],
        [
            CEG.ChangePositionSystem(setupNaCHA.molecule, pos3),
            CEG.ChangePositionSystem(setupNaCHA.molecule, pos2),
        ]
    ])
    @test compute_ewald(ctx) ≈ compute_ewald(ctx_twoco2)
end

@testset "MonteCarloSetup" begin
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
    mat = NoUnits.(setupArCIT7.grids[1].csetup.cell.mat./u"Å")
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
    mcArCIT7, _ = setup_montecarlo("CIT7", "BoulfelfelSholl2021", [molAr]; blockfiles=[false]);
    mcVoidArCIT7, _ = setup_montecarlo(supermat, "BoulfelfelSholl2021", [molAr, superCIT7]);
    CEG.baseline_energy(mcArCIT7); movArCIT7 = CEG.movement_energy(mcArCIT7, (1,1))
    CEG.baseline_energy(mcVoidArCIT7); movVoidArCIT7 = CEG.movement_energy(mcVoidArCIT7, (1,1))
    @test Float64(movArCIT7) ≈ Float64(movVoidArCIT7) rtol=0.005

    molAr1 = CEG.ChangePositionSystem(ar, [SVector{3}([-7.7365250811304911,31.5070011601372251,1.5285305931479920]u"Å")]);
    molAr2 = CEG.ChangePositionSystem(ar, [SVector{3}([10.7586599791867421,-2.3259182727570948,20.5642722996513001]u"Å")]);
    mcAr2, _ = setup_montecarlo("CIT7", "BoulfelfelSholl2021", [molAr1, molAr2]; blockfiles=[false, false]);
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
    @test baseMini ≈ Float64(CEG.baseline_energy(mcNaMiniRef)) rtol=0.001
    @test baseMini ≈ -248304.58180794 rtol=0.001

    molNaPetit = CEG.ChangePositionSystem(na, [SVector{3}([18.77838182689036, 14.73031622108175, 3.669308624409278]u"Å")]);
    mcNaPetit, _ = setup_montecarlo("Petit", "BoulfelfelSholl2021", [molNaPetit]);
    mcNaPetitRef, _ = setup_montecarlo("PetitRef", "BoulfelfelSholl2021", [molNaPetit]);
    basePetit = Float64(CEG.baseline_energy(mcNaPetit))
    @test basePetit ≈ Float64(CEG.baseline_energy(mcNaPetitRef)) rtol=0.0001
    @test basePetit ≈ 262204.85720076 rtol=0.001

    molNaSolo = CEG.ChangePositionSystem(na, [SVector{3}([-4.728415488310421, 32.03533696753957, 2.943765448968882]u"Å")]);
    mcNaSolo, _ = setup_montecarlo("CIT7", "BoulfelfelSholl2021", [molNaSolo]);
    @test mcNaSolo.tailcorrection[] ≈ -70.44772635984882u"K" # account for supercell
    baseSolo = Float64(CEG.baseline_energy(mcNaSolo))
    @test baseSolo ≈ -21375.116833457894 rtol=0.001
    posSoloNext = [SVector{3}([-5.036,31.876,3.117]u"Å")]
    molNaSoloNext = CEG.ChangePositionSystem(na, posSoloNext);
    mcNaSoloNext, _ = setup_montecarlo("CIT7", "BoulfelfelSholl2021", [molNaSoloNext]);
    @test mcNaSolo.tailcorrection[] == mcNaSoloNext.tailcorrection[]
    baseSoloNext = Float64(CEG.baseline_energy(mcNaSoloNext))
    @test baseSoloNext ≈ -21795.8765195143 rtol=0.001
    @test baseSoloNext ≈ baseSolo + Float64(CEG.movement_energy(mcNaSolo, (1,1), posSoloNext) - CEG.movement_energy(mcNaSolo, (1,1)))

    molNaDuo1 = CEG.ChangePositionSystem(na, [SVector{3}([1.001641978413878, 8.638263743446769, 17.23737576131632]u"Å")]);
    molNaDuo2 = CEG.ChangePositionSystem(na, [SVector{3}([4.163424680441308, 16.12704796355876, 17.20994387427006]u"Å")]);
    mcNaDuo, _ = setup_montecarlo("CIT7", "BoulfelfelSholl2021", [molNaDuo1, molNaDuo2]);
    @test Float64(CEG.baseline_energy(mcNaDuo)) ≈ -25957.746610866408 rtol=0.001

    co2 = CEG.load_molecule_RASPA("CO2", "TraPPE", "BoulfelfelSholl2021");
    molNaTrio = CEG.ChangePositionSystem(na, [SVector{3}([3.019388765467742, 0.8997706038543032, 26.11901621898599]u"Å")]);
    molCO2_1 = CEG.ChangePositionSystem(co2, SVector{3}.([[11.93940309885289, 8.48657378465003, 2.135736631609201]u"Å",
                                                         [11.10485516124311, 7.710040763525694, 1.991767166323031]u"Å",
                                                         [10.27030722363334, 6.933507742401357, 1.84779770103686]u"Å"]));
    molCO2_2 = CEG.ChangePositionSystem(co2, SVector{3}.([[5.491645446274333, 8.057854365959964, 8.669190836544463]u"Å",
                                                         [6.335120278303245, 7.462084936052019, 9.172986424179925]u"Å",
                                                         [7.178595110332157, 6.866315506144074, 9.676782011815387]u"Å"]));
    mcTrio, _ = setup_montecarlo("CIT7", "BoulfelfelSholl2021", [molNaTrio, molCO2_1, molCO2_2]; blockfiles=[false, false, false]);
    baseTrio = Float64(CEG.baseline_energy(mcTrio))
    @test baseTrio ≈ -28329.113561030445 rtol=0.001
    newposNaTrio = [SVector{3}([0.9, 0.1, 1.5]u"Å")]
    diffNa = Float64(CEG.movement_energy(mcTrio, (1,1), newposNaTrio) - CEG.movement_energy(mcTrio, (1,1)))
    @test diffNa ≈ 5440.529635958557 rtol=0.001
    mcTrio_diffNa, _ = setup_montecarlo("CIT7", "BoulfelfelSholl2021", [CEG.ChangePositionSystem(na, newposNaTrio), molCO2_1, molCO2_2]; blockfiles=[false, false, false]);
    @test baseTrio + diffNa ≈ Float64(CEG.baseline_energy(mcTrio_diffNa))
    newposCO2_2 = SVector{3}.([[1.3, 2.9, 1.149]u"Å", [1.3, 2.9, 0.0]u"Å", [1.3, 2.9, -1.149]u"Å"])
    diffCO2 = Float64(CEG.movement_energy(mcTrio, (2,2), newposCO2_2) - CEG.movement_energy(mcTrio, (2,2)))
    mcTrio_diffCO2, _ = setup_montecarlo("CIT7", "BoulfelfelSholl2021", [molNaTrio, molCO2_1, CEG.ChangePositionSystem(molCO2_1, newposCO2_2)]; blockfiles=[false, false, false]);
    @test baseTrio + diffCO2 ≈ Float64(CEG.baseline_energy(mcTrio_diffCO2)) rtol=0.0001

    # blocking sphere straddling a periodic boundary
    setupArCIT7block = setup_RASPA("CIT7block", "BoulfelfelSholl2021", "Ar", "TraPPE");
    @test energy_point(setupArCIT7block, [SVector{3}([12.5, 0.8, 0.3]u"Å")]) == (1e100u"K", 0.0u"K")
    @test energy_point(setupArCIT7block, [SVector{3}([20.175808361078516, 10.58027451750961, 9.185277912816744]u"Å")]) == (1e100u"K", 0.0u"K") # same in another cell
end

@testset "TailCorrection" begin
    co2 = CEG.load_molecule_RASPA("CO2", "TraPPE", "BoulfelfelSholl2021");
    na = CEG.load_molecule_RASPA("Na", "TraPPE", "BoulfelfelSholl2021");
    mc, _ = setup_montecarlo("CIT7", "BoulfelfelSholl2021", [(na, 2), (co2, 3)]; blockfiles=[false,false]);

    ff = mc.step.ff
    ffidx = mc.step.ffidx
    @test getindex.((ff.symbols,), ffidx) == [[:Na], [:O, :C, :O]]
    framework_atoms = [0, 720, 0, 0, 360, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    λ = ustrip(u"Å^-3", 2π/det(mc.step.mat))

    tc0 = CEG.TailCorrection(ff, ffidx, framework_atoms, λ, [0,0])
    @test tc0[] == 0.0u"K"
    CEG.modify_species!(tc0, 1, 2)
    CEG.modify_species!(tc0, 2, 3)
    @test tc0[] ≈ mc.tailcorrection[] ≈ -322.46442841047406u"K"

    tc1 = CEG.TailCorrection(ff, ffidx, framework_atoms, λ, [2,3])
    @test tc1[] ≈ mc.tailcorrection[]

    tc2 = CEG.TailCorrection(ff, ffidx, framework_atoms, λ, [1,5])
    CEG.modify_species!(tc2, 2, -1)
    CEG.modify_species!(tc2, 1, 1)
    CEG.modify_species!(tc2, 2, -1)
    @test tc2[] ≈ mc.tailcorrection[]
end

@testset "Deletion and addition" begin
    ar = CEG.load_molecule_RASPA("Ar", "TraPPE", "BoulfelfelSholl2021");
    co2 = CEG.load_molecule_RASPA("CO2", "TraPPE", "BoulfelfelSholl2021");
    molCO2_1 = CEG.ChangePositionSystem(co2, SVector{3}.([[4.221014336721, 2.235273775272, 5.118873160667]u"Å",
                                                         [4.915895823098, 3.103469549355, 4.829776606287]u"Å",
                                                         [5.610777309474, 3.971665323438, 4.540680051907]u"Å"]));
    molCO2_2 = CEG.ChangePositionSystem(co2, SVector{3}.([[13.732840754800, 5.744459327798, 8.119705823930]u"Å",
                                                         [13.418141549103, 6.398573841906, 7.229032139371]u"Å",
                                                         [13.103442343406, 7.052688356015, 6.338358454812]u"Å"]));
    molCO2_3 = CEG.ChangePositionSystem(co2, SVector{3}.([[21.473519181736, 19.719777629492, 13.809504713725]u"Å",
                                                         [21.888511731126, 20.756701257738, 13.539742646090]u"Å",
                                                         [22.303504280517, 21.793624885985, 13.269980578454]u"Å"]));
    pos1 = position(molCO2_1)
    pos2 = position(molCO2_2)
    posX = SVector{3}.([[5.488965064161, 14.335087694715, 16.087580364999]u"Å",
                        [4.887655471102, 13.508430918264, 15.562922050243]u"Å",
                        [4.286345878043, 12.681774141812, 15.038263735487]u"Å"])

    mcRef, _ = setup_montecarlo("CHA_1.4_3b4eeb96_Na_11812", "BoulfelfelSholl2021", [molCO2_1, molCO2_3, ar]);
    baseRef = baseline_energy(mcRef)
    moveRef = movement_energy(mcRef, (1,1), pos2)
    CEG.update_mc!(mcRef, (1, 1), pos2)
    othermoveRef = movement_energy(mcRef, (1,2), posX)

    for (mols, i, j) in (([ar, molCO2_1, molCO2_2, molCO2_3], 2, 2),
                         ([ar, molCO2_2, molCO2_1, molCO2_3], 2, 1),
                         ([molCO2_1, molCO2_2, molCO2_3, ar], 1, 2),
                         ([molCO2_2, molCO2_1, molCO2_3, ar], 1, 1))
        @info "Testing $i $j"

        mcA, _ = setup_montecarlo("CHA_1.4_3b4eeb96_Na_11812", "BoulfelfelSholl2021", mols);
        mcB = deepcopy(mcA)
        CEG.remove_one_system!(mcA, i, j)
        @test baseline_energy(mcA) ≈ baseRef
        @test movement_energy(mcA, (i, 3-j), pos2) ≈ moveRef
        CEG.update_mc!(mcA, (i, 3-j), pos2)
        @test movement_energy(mcA, (i, j), posX) ≈ othermoveRef
        @test CEG.add_one_system!(mcA, i, pos1) == 3
        CEG.remove_one_system!(mcA, i, 3-j)
        @test movement_energy(mcA, (i, 3-j), pos2) ≈ moveRef

        baseline_energy(mcB)
        CEG.remove_one_system!(mcB, i, j)
        @test movement_energy(mcB, (i, 3-j), pos2) ≈ moveRef
        CEG.update_mc!(mcB, (i, 3-j), pos2)
        @test movement_energy(mcB, (i, j), posX) ≈ othermoveRef
        @test CEG.add_one_system!(mcB, i, pos1) == 3
        CEG.remove_one_system!(mcB, i, 3-j)
        @test movement_energy(mcB, (i, 3-j), pos2) ≈ moveRef
    end
end

@testset "Random walk" begin
    na = CEG.load_molecule_RASPA("Na", "TraPPE", "BoulfelfelSholl2021");
    co2 = CEG.load_molecule_RASPA("CO2", "TraPPE", "BoulfelfelSholl2021");
    molNaTrio = CEG.ChangePositionSystem(na, [SVector{3}([3.019388765467742, 0.8997706038543032, 26.11901621898599]u"Å")]);
    molCO2_1 = CEG.ChangePositionSystem(co2, SVector{3}.([[11.93940309885289, 8.48657378465003, 2.135736631609201]u"Å",
                                                         [11.10485516124311, 7.710040763525694, 1.991767166323031]u"Å",
                                                         [10.27030722363334, 6.933507742401357, 1.84779770103686]u"Å"]));
    molCO2_2 = CEG.ChangePositionSystem(co2, SVector{3}.([[5.491645446274333, 8.057854365959964, 8.669190836544463]u"Å",
                                                         [6.335120278303245, 7.462084936052019, 9.172986424179925]u"Å",
                                                         [7.178595110332157, 6.866315506144074, 9.676782011815387]u"Å"]));
    mcTrio, _ = setup_montecarlo("CIT7", "BoulfelfelSholl2021", [molNaTrio, molCO2_1, molCO2_2]; blockfiles=[false, false, false]);
    reportsTrio = run_montecarlo!(mcTrio, SimulationSetup(300u"K", 4000))
    shadow, _ = setup_montecarlo("CIT7", "BoulfelfelSholl2021", [CEG.ChangePositionSystem(na, [mcTrio.step.positions[1]]), CEG.ChangePositionSystem(co2, mcTrio.step.positions[2:4]), CEG.ChangePositionSystem(co2, mcTrio.step.positions[5:7])]; blockfiles=[false,false,false]);
    baseTrio = Float64(baseline_energy(mcTrio))
    @test baseTrio ≈ Float64(baseline_energy(shadow))
    @test baseTrio ≈ Float64(reportsTrio[end])

    #= Update protocol:
    In the tests below, if the results fails because the simulation protocol was updated,
    the values must be updated as well. To do so, run the Monte-Carlo simulation, then
    export a restart file and launch RASPA from the restart file to extract the comparison
    values.
    For instance, to regenerate mc2, run in julia:

    ```julia
    framework2 = "CHA_1.4_3b4eeb96"
    mc2, _ = setup_montecarlo("CHA_1.4_3b4eeb96", "BoulfelfelSholl2021", [(na, 135),]; rng=StableRNG(2));
    reports2 = run_montecarlo!(mc2, SimulationSetup(300.0u"K", 1000; outdir="", printevery=1000))
    mkpath("/tmp/mc2/RestartInitial/System_0/")
    CEG.output_restart("/tmp/mc2/RestartInitial/System_0/restart_$(framework2)_1.1.1_300.000000_2e+08", mc2)
    ```

    then, from the "/tmp/mc2" directory, create a "simulation.input" file containing the
    following (after putting the appropriate molecule name and framework)

    ```
    SimulationType                MonteCarlo
    NumberOfCycles                0
    NumberOfInitializationCycles  0
    PrintEvery                    1
    RestartFile                   yes

    Forcefield                    BoulfelfelSholl2021

    Component 0 MoleculeName                     Na
                MoleculeDefinition               TraPPE
                ExtraFrameworkMolecule           no

    Framework 0
    ChargeMethod Ewald
    EwaldPrecision 1e-06
    FrameworkName CHA_1.4_3b4eeb96
    CutOffVDW 12.0
    UnitCells 1 1 1
    ExternalPressure 2.0e8
    RemoveAtomNumberCodeFromLabel yes
    Movies yes
    WriteMoviesEvery 1
    ExternalTemperature 300
    ```

    then run it with:

    ```bash
    $RASPA_DIR/bin/simulate -i simulation.input
    ```

    and finally fetch the result with:

    ```bash
    grep -m1 -A70 "Host/Adsorbate energy:" Output/System_0/*.data
    ```
    =#

    framework2 = "CHA_1.4_3b4eeb96"
    mc2, _ = setup_montecarlo(framework2, "BoulfelfelSholl2021", [(na, 135),]; rng=StableRNG(2));
    reports2 = run_montecarlo!(mc2, SimulationSetup(300.0u"K", 1000; outdir="", printevery=1000))
    @test length(reports2) == 2
    @test reports2[end] ≈ CEG.BaselineEnergyReport(-718484.79640543, -6214557.02749673, 626259.86449476, -248250821.39479709+121425571.41820429, -8575.23090402) rtol=0.0001

    ar = CEG.load_molecule_RASPA("Ar", "TraPPE", "BoulfelfelSholl2021");
    framework3 = "CHA_1.4_3b4eeb96_Na_11812"
    mc3, _ = setup_montecarlo(framework3, "BoulfelfelSholl2021", [(ar, 30)]; rng=StableRNG(3));
    reports3 = run_montecarlo!(mc3, CEG.SimulationSetup(300.0u"K", 1000; outdir="", printevery=0))
    @test reports3[end] ≈ CEG.BaselineEnergyReport(-42195.81574644, 0.0, -563.48165902, 0.0, -8575.23090402) rtol=0.0001

    framework4 = "FAU_1.4_03f7b3ba"
    mc4, _ = setup_montecarlo(framework4, "BoulfelfelSholl2021", [(na, -1)]; rng=StableRNG(4));
    reports4 = run_montecarlo!(mc4, SimulationSetup(300.0u"K", 1000, outdir="", printevery=0))
    @test reports4[end] ≈ CEG.BaselineEnergyReport(-422811.87279224, -3825797.13918730, 390589.13414715, -101231726.22144973+49027959.61839646, -4782.44322153) rtol=0.0003
end

@testset "Restart" begin
    co2 = CEG.load_molecule_RASPA("CO2", "TraPPE", "BoulfelfelSholl2021");
    path = joinpath(@__DIR__, "CHA_1.4_3b4eeb96_Na_11812.restart")
    mcRestart, _ = setup_montecarlo("CHA_1.4_3b4eeb96_Na_11812", "BoulfelfelSholl2021", [(co2,1)]; restart=path)
    @test Float64(baseline_energy(mcRestart)) ≈ -14128.888030042883 rtol=0.001
    CEG.output_restart(path*".COPY", mcRestart)
    @test read(path) == read(path*".COPY")
    rm(path*".COPY")
end

CEG.PRINT_CHARGE_WARNING[] = prev_warning
