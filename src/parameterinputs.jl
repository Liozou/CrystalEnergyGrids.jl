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

mutable struct RMinimumEnergy{N,T} <: RecordFunction
    mine::BaselineEnergyReport
    minpos::SimulationStep{N,T}

    RMinimumEnergy{N,T}() where {N,T} = new{N,T}(BaselineEnergyReport(Inf*u"K", Inf*u"K", Inf*u"K", Inf*u"K", Inf*u"K"))
    function RMinimumEnergy(mine::BaselineEnergyReport, minpos::SimulationStep{N,T}) where {N,T}
        new{N,T}(mine, minpos)
    end
end
function (record::RMinimumEnergy)(o::SimulationStep, e::BaselineEnergyReport, k::Int,
                                  mc::MonteCarloSetup, simu::SimulationSetup)
    if Float64(e) < Float64(record.mine)
        record.mine = e
        record.minpos = o
    end
    if k == simu.ncycles && !isempty(simu.outdir)
        output_pdb(joinpath(simu.outdir, "min_energy.pdb"), mc, record.minpos)
        output_restart(joinpath(simu.outdir, "min_energy.restart"), record.minpos)
    end
    nothing
end


struct ShootingStarMinimizer{N,T} <: RecordFunction
    every::Int
    length::Int
    positions::Vector{SimulationStep{N,T}}
    energies::Vector{BaselineEnergyReport}
    lb::LoadBalancer{Tuple{Int,MonteCarloSetup{N,T},SimulationSetup{RMinimumEnergy{N,T}}}}
end
function ShootingStarMinimizer{N}(nsteps::Int; length::Int=100, every::Int=1) where {N}
    n = nsteps ÷ every
    T = typeof_psystem(Val(N))
    positions = Vector{SimulationStep{N,T}}(undef, n)
    energies = Vector{BaselineEnergyReport}(undef, n)
    lb = LoadBalancer{Tuple{Int,MonteCarloSetup{N,T},SimulationSetup{RMinimumEnergy{N,T}}}}(nthreads()) do (ik, newmc, newsimu)
        let ik=ik, newmc=newmc, newsimu=newsimu, positions=positions, energies=energies
            run_montecarlo!(newmc, newsimu)
            positions[ik] = newsimu.record.minpos
            energies[ik] = newsimu.record.mine
        end
    end
    ShootingStarMinimizer{N,T}(every, length, positions, energies, lb)
end

function (star::ShootingStarMinimizer)(o::SimulationStep, e::BaselineEnergyReport, k::Int, mc::MonteCarloSetup, _)
    k == 0 && return
    ik, r = divrem(k, star.every)
    r == 0 || return
    newmc = MonteCarloSetup(mc; parallel=false)
    recordminimum = RMinimumEnergy(e, o)
    newsimu = SimulationSetup(300u"K", star.length; printevery=0, record=recordminimum)
    put!(star.lb, (ik, newmc, newsimu))
    nothing
end

Base.fetch(x::ShootingStarMinimizer) = wait(x.lb)


struct RainfallMinimizer{N,T} <: RecordFunction
    every::Int
    length::Int
    positions::Vector{SimulationStep{N,T}}
    energies::Vector{BaselineEnergyReport}
    tasks::Vector{Task}
end

function RainfallMinimizer{N}(nsteps::Int; length::Int=100, every::Int=1) where {N}
    n = nsteps ÷ every
    T = typeof_psystem(Val(N))
    positions = Vector{SimulationStep{N,T}}(undef, n)
    energies = Vector{BaselineEnergyReport}(undef, n)
    tasks = Vector{Task}(undef, n)
    RainfallMinimizer{N,T}(every, length, positions, energies, tasks)
end

function (rain::RainfallMinimizer)(o::SimulationStep, e::BaselineEnergyReport, k::Int, mc::MonteCarloSetup, _)
    k == 0 && return
    ik, r = divrem(k, rain.every)
    r == 0 || return
    newmc = MonteCarloSetup(mc; parallel=false)
    recordminimum = RMinimumEnergy(e, o)
    newsimu = SimulationSetup(300u"K", rain.length; printevery=0, record=recordminimum)
    rain.tasks[ik] = @spawn begin
        run_montecarlo!($newmc, $newsimu)
        $rain.positions[$ik] = $newsimu.record.minpos
        $rain.energies[$ik] = $newsimu.record.mine
    end
end

Base.fetch(x::RainfallMinimizer) = foreach(wait, x.tasks)
