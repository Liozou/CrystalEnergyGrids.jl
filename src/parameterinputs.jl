using Base.Threads

## Temperature input
abstract type TemperatureFunction <: Function end

struct TRamp <: TemperatureFunction
    start::TK
    Δ::TK

    function TRamp(start::TK, Δ::TK, ::Nothing)
        new(start, Δ)
    end
end
function TRamp(start, stop)
    TRamp(convert(TK, start), convert(TK, stop-start), nothing)
end
function (tf::TRamp)(k, n)
    tf.start + tf.Δ*(k-1)/(n-1)
end

struct TPoly <: TemperatureFunction
    ramp::TRamp
    γ::Float64
end
function (tf::TPoly)(k, n)
    tf.ramp.start + tf.ramp.Δ*((k-1)/(n-1))^tf.γ
end

struct TSteps <: TemperatureFunction
    ramp::TRamp
    nstepsp1::Int
end
function (tf::TSteps)(k, n)
    tf.ramp.start + tf.Δ*floor(tf.nstepsp1*(k-1)/(n-1))/tf.nstepsp1
end

struct TParts{N,T} <: TemperatureFunction
    fractions::NTuple{N,Float64} # sorted between 0 and 1, the last element is always 1
    parts::T
end
function (tf::TParts{N})(k, n) where N
    ρ = (k-1)/(n-1)
    for i in 1:N
        frac = tf.fractions[i]
        if frac ≥ ρ
            num_start = i==1 ? -1 : floor(Int, n*tf.fractions[i-1])-1
            num_stop = floor(Int, n*frac)
            num_current = floor(Int, n*ρ)
            return tf.parts[i](num_current - num_start, num_stop - num_start)
        end
    end
    @assert ("Please report this error."; false)
    0.0u"K"
end

function TAnneal(bottom::TK, top::TK, wait::Float64)
    TParts(((1.0 - wait)/2, (1.0 + wait)/2, 1.0), (
        TRamp(bottom, top),
        Returns(top),
        TRamp(top, bottom)
    ))
end
function TAnneal(sequence, wait::Number)
    N = length(sequence)
    TParts(ntuple(i -> begin
                        j, o = fldmod(i, 2)
                        (j + o*(1-wait))/(N-1)
                       end, 2*(N-1)),
           ntuple(i -> begin
                        j, o = fldmod1(i, 2)
                        seqjp1 = sequence[j+1]
                        o == 1 ? TRamp(sequence[j], seqjp1) : Returns(seqjp1)
                       end, 2*(N-1))
    )
end


## Record inputs
abstract type RecordFunction <: Function end

mutable struct RMinimumEnergy{T} <: RecordFunction
    mink::Int
    mine::BaselineEnergyReport
    minpos::SimulationStep{T}

    RMinimumEnergy{T}() where {T} = new{T}(BaselineEnergyReport(0, Inf*u"K", Inf*u"K", Inf*u"K", Inf*u"K", Inf*u"K"))
    function RMinimumEnergy(mink::Int, mine::BaselineEnergyReport, minpos::SimulationStep{T}) where T
        new{T}(mink, mine, minpos)
    end
end
function (record::RMinimumEnergy)(o::SimulationStep, e::BaselineEnergyReport, k::Int,
                                  mc::MonteCarloSetup, simu::SimulationSetup)
    if Float64(e) < Float64(record.mine)
        record.mink = k
        record.mine = e
        record.minpos = o
    end
    if k == simu.ncycles && !isempty(simu.outdir)
        open(Base.Fix2(println, record.mink), joinpath(simu.outdir, "min_k"), "w")
        output_pdb(joinpath(simu.outdir, "min_energy.pdb"), mc, record.minpos)
        output_restart(joinpath(simu.outdir, "min_energy.restart"), record.minpos)
    end
    nothing
end

_empty_elastic_matrix!(x::ElasticMatrix) = resize!(x, size(x, 1), 0)

struct ShootingStarMinimizer{T} <: RecordFunction
    temperature::TK
    every::Int
    length::Int
    positions::Vector{SimulationStep{T}}
    energies::Vector{BaselineEnergyReport}
    ninit::Int
    outdir::String
    printevery::Int
    printeveryinit::Int
    outtype::Vector{Symbol}
    lb::LoadBalancer{Tuple{Int,MonteCarloSetup{T},SimulationSetup{RMinimumEnergy{T}}}}
end
function ShootingStarMinimizer(; length::Int=100, every::Int=1, outdir="", printevery::Integer=0, printeveryinit::Integer=printevery, outtype::AbstractVector{Symbol}=[:energies, :zst], temperature=300u"K", ninit=20_000, numthreads=nthreads())
    T = typeof_psystem(Val(3))
    positions = Vector{SimulationStep{T}}(undef, 0)
    energies = Vector{BaselineEnergyReport}(undef, 0)
    lb = LoadBalancer{Tuple{Int,MonteCarloSetup{T},SimulationSetup{RMinimumEnergy{T}}}}(numthreads) do (ik, newmc, newsimu), _
        let ik=ik, newmc=newmc, newsimu=newsimu, positions=positions, energies=energies
            run_montecarlo!(newmc, newsimu)
            positions[ik] = newsimu.record.minpos
            energies[ik] = newsimu.record.mine

            # GC maintenance: force freeing the buffers to avoid memory leak
            _empty_elastic_matrix!(newmc.ewald.sums)
            foreach(_empty_elastic_matrix!, newmc.ewald.tmpEiks)
            empty!(newmc.ewald.tmpsums)
            GC.gc()
        end
    end
    ShootingStarMinimizer(temperature, every, length, positions, energies, ninit, outdir, printevery, printeveryinit, outtype, lb)
end
function initialize_record!(star::T, simu::SimulationSetup{T}) where {T <: ShootingStarMinimizer}
    n = simu.ncycles ÷ star.every
    resize!(star.positions, n)
    resize!(star.energies, n)
end

function (star::ShootingStarMinimizer)(o::SimulationStep, e::BaselineEnergyReport, k::Int, mc::MonteCarloSetup, _)
    k ≤ 0 && return
    ik, r = divrem(k, star.every)
    r == 0 || return
    newmc = MonteCarloSetup(mc, o; parallel=false)
    recordminimum = RMinimumEnergy(0, e, o)
    outdir = isempty(star.outdir) ? "" : joinpath(star.outdir, string(ik))
    newsimu = SimulationSetup(star.temperature, star.length; outdir, star.printevery, star.printeveryinit, star.outtype, star.ninit, record=recordminimum)
    put!(star.lb, (ik, newmc, newsimu))
    nothing
end

Base.fetch(x::ShootingStarMinimizer) = wait(x.lb)


function reconstitute_trace(path::AbstractString, skip, keep)
    numdirs = length(readdir(path)) - 1
    energies = [deserialize(joinpath(path, string(i), "energies.serial")) for i in 1:numdirs if begin
        x = ispath(joinpath(path, string(i), "energies.serial"))
        x || @warn "Subfolder $i does not contain energies.serial"
        x
    end]
    numdirs = length(energies)
    trace = Vector{Vector{Float64}}(undef, numdirs)
    for (i, energy) in enumerate(energies)
        n = length(energy)
        start = skip isa Integer ? skip : floor(Int, skip*n)
        start += (start == 0)
        stop = keep isa Integer ? start + keep - 1 : start + round(Int, keep*(n-start+1))
        stop -= (stop == length(energy) + 1)
        trace[i] = Float64.(@view energy[start:stop])
    end
    trace
end


## GCMC

export run_gcmc!

struct GCMCRecorder{T}
    Ns::Vector{Float64}
    subrecord::T

    function GCMCRecorder(simu::SimulationSetup{T}) where {T}
        new{T}(Int[], simu.record)
    end
    GCMCRecorder() = new{Returns{Nothing}}(Int[], Returns(nothing))
end

function (rec::GCMCRecorder{T})(o, e, k, mc, simu) where T
    if k > 0
        rec.Ns[k] = length(o.posidx[1])/mc.gcmcdata.Π
    end
    if !(T <: Returns{Nothing})
        rec.subrecord(o, e, k, mc, simu)
    end
end
function initialize_record!(rec::T, simu::SimulationSetup{T}) where {T <: GCMCRecorder}
    resize!(rec.Ns, simu.ncycles)
    fill!(rec.Ns, -1)
    # no need to initialize_record!(rec.subrecord) since that was already done before creating GCMCRecorder
    nothing
end

function run_gcmc!(mc::MonteCarloSetup, simu::SimulationSetup)
    rec = GCMCRecorder(simu)
    mcmoves1 = mc.mcmoves[1]
    if iszero(mcmoves1[:swap])
        iswap = findfirst(==(:swap), mcmovenames)
        newmcmoves1 = MCMoves(ntuple(i -> i < iswap ? mcmoves1.cumulatives[i]*2/3 : mcmoves1.cumulatives[i], Val(length(mcmovenames)-1))) # make swap probability 1/3
        mc = MonteCarloSetup(mc; mcmoves=[newmcmoves1; mc.mcmoves[2:end]])
    end
    printstyled("Running GCMC on species ", identify_molecule(mc.step.ff.symbols[mc.step.ffidx[1]]), " at pressure ", simu.pressure, '\n'; color=:green)
    newsimu = SimulationSetup(simu.temperatures, simu.pressure, simu.ncycles, simu.ninit, simu.outdir, simu.printevery, simu.printeveryinit, simu.outtype, rec)
    ret = run_montecarlo!(mc, newsimu)
    ret, rec.Ns
end

function make_isotherm(mc::MonteCarloSetup, simu::SimulationSetup, pressures)
    n = length(pressures)
    isotherm = Vector{Float64}(undef, n)
    energies = Vector{BaselineEnergyReport}(undef, n)
    @assert allequal(simu.temperatures)
    @info "The number of cycles is multiplied by 1 + 1e3/sqrt(ustrip(u\"Pa\", pressure))"
    @threads for i in 1:n
        newmc = MonteCarloSetup(mc; parallel=false)
        pressure = pressures[i]
        P = pressure isa Quantity ? uconvert(u"Pa", pressure) : pressure*u"Pa"
        ncycles = ceil(Int, simu.ncycles*(1+1e3/sqrt(ustrip(u"Pa", P))))
        temperatures = fill(simu.temperatures[1], ncycles+simu.ninit)
        newsimu = SimulationSetup(temperatures, P, ncycles, simu.ninit, simu.outdir, simu.printevery, simu.printeveryinit, simu.outtype, simu.record)
        es, Ns = run_gcmc!(newmc, newsimu)
        isotherm[i] = mean(Ns)
        energies[i] = mean(es)
        println("Finished $pressure Pa")
    end
    isotherm, energies
end
