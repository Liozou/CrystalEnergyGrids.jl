export SimulationSetup, run_montecarlo!

"""
    SimulationSetup{Trecord}

Definition of the simulation setup, which must provide the following information:
- `T` the temperature. Either a number (in K), a list of temperatures (in K, one per cycle)
  or a function. If a function, it must take two parameters: `j` the number of the current
  cycle and `n` the total number of cycles, and return a temperature (in K).
  Internally, the list of temperatures is stored as a field `temperatures` indexed by `j`,
  which varies between `1` and `n`.
- `pressure` the pressure (in Pa). Only used for GCMC, for which the Monte-Carlo moves must
  include a non-zero swap probability.
- `ncycles` the number of cycles of the simulation. Each cycle consists in `N` steps, where
  a step is a Monte-Carlo movement attempt on a random species, and `N` is the number of
  species.
- `ninit` the number of initialization cycles, where no output is written and no energy
  stored. The total number of cycles `n` used in the definition of `T` corresponds to
  `ncycles + ninit`. The function `record` defined below is still called during the
  initialization cycles.
- `outdir` the path to a directory in which the result is stored. If not provided (or set
  to an empty string), no record will be stored on file.
- `printevery`: one output of the positions is stored at the end of every `printevery`
  cycles. The first recorded position is that at the end of initialization. Using a value
  of 0 (the default) will record only that, as well as the last positions if `ncycles != 0`.
  It is also possible to output positions during the initialization (see `outtype` below):
  the only difference is that a value of `0` for `printevery` is then equivalent to `1`
- `outtype`: a `Vector{Symbol}` whose elements represent the different outputs that
  should be produced in `outdir`. The elements can be:
  + `:energies` to have the list of energies output as a serialized
    `Vector{BaselineEnergyReport}` in a "energies.serial" file.
  + `:pdb` to have the positions of the atoms output as a "trajectory.pdb" file.
  + `:stream` to have the positions of the atoms output as a [`StreamSimulationStep`](@ref)
    in a "steps.stream" file.
  + `:zst` to have the positions of the atoms output as a compressed
    [`StreamSimulationStep`](@ref) in a "steps.zst" file.
  + `:initial_pdb`, `:initial_stream`, `:initial_zst` to have the corresponding output even
    during initialization. The corresponding files have an "initial_" prefix.
  The default is `[:energies, :zst]`. See [`stream_to_pdb`] to convert a "steps.stream"
  or "steps.zst" file into the corresponding "trajectory.pdb" output.
- `record::Trecord` a function which can be used to record extra information at the end of
  each cycle. `record` can return `:stop` to end the computation at that point, otherwise
  its return value is ignored. Stopping may not be immediate and can take up to one
  additional cycle. `record` is called even during the initialization cycles.

If provided, `record` must have the signature
    `record(o::SimulationStep, e::BaselineEnergyReport, k::Int, mc::MonteCarloSetup, simu::SimulationSetup)`
where
- `o` contains the information on the positions of the species at this point. Read below
  for more information about what can be done with `o` depending on `needcomplete`
- `e` is the energy at this point.
- `k` is the cycle number (between `-ninit` and `ncycles`). Negative values of `k`
  correspond to the initialization cycles. `k == -ninit` corresponds to the starting
  configuration, before the initialization cycles. `k == 0` corresponds to the state before
  the start of the production cycles. `k == ncycles` is the last state of the simulation.
  !!! tip
      The cycle number `k` corresponds to the temperature at index `j = k + ninit`.
- `mc` is the input `MonteCarloSetup`. Note that the `positions` field may be incorrect,
  use that of `o`. It must not be modified by the `record` function.
- `simu` is this `SimulationSetup`. It must not be modified by the `record` function.

By default, the `record` function is given as an argument a value of `o` which is a copy
of `mc.step` which is thread-safe to read and modify (except for the `ff`, `charges`,
`isrigid` and `ffidx` fields which are never meant to change).
If the `record` function is guaranteed to only ever read or modify the `positions` and
`atoms` fields of `o`, then it can specify it by defining a method for
`needcomplete(::Trecord)` that returns `false` (the default is `true`).
If it reads any other fields, like `posidx`, no `needcomplete` method should be defined.

!!! warning
    If the value returned by `needcomplete` is `false` when it should be `true`, this may
    lead to race conditions and related subtle bugs, which may not be easy to identify.

Additionally, if a call to `record` can spawn asynchronous calls, a method for `Base.fetch`
should be implemented such that `Base.fetch(record)` blocks until all asynchronous calls
are done.

Upon creation of a `SimulationSetup`, a call to `initialize_record!(record, setup)` will be
realized. By default, `initialize_record!` does nothing, but you can specialize it for
for your appropriate record type if need be, i.e. by defining a method of signature:
`initialize_record!(record::T, setup::SimulationSetup{T}) where {T<:MyRecordType}`.

## Example

This is an example of a record function that stores the averages of the distance between
the particles and a given point, weighted by the current temperature:

```julia
struct AverageDistanceTemperatureRecord <: Function
    refpoint::SVector{3,typeof(1.0u"Å")}
    weightedaverages::Vector{typeof(1.0u"Å * K")}
end
AverageDistanceTemperatureRecord(refpoint) = AverageDistanceTemperatureRecord(refpoint, [])

function (x::AverageDistanceTemperatureRecord)(o, _, k, _, simu)
    k ≤ 0 && return # skip initialization cycles
    buffer, ortho, safemin = prepare_periodic_distance_computations(o.cell)
    weightedaverage = 0.0u"Å * K"
    for pos in o.positions
        buffer .= pos .- x.refpoint
        weightedaverage += periodic_distance!(buffer, o.cell, ortho, safemin)
    end
    x.weightedaverages[k] = weightedaverage * simu.temperatures[k] / length(o.atoms)
    nothing
end

# Initialize the internal `weightedaverages` vector to the right size
function CrystalEnergyGrids.initialize_record!(record::T, setup::SimulationSetup{T}) where {T<:AverageDistanceTemperatureRecord}
    resize!(record.weightedaverages, setup.ncycles)
end

# No need to give a complete copy of `o` since we only read the positions
CrystalEnergyGrids.needcomplete(::AverageDistanceTemperatureRecord) = false

# (optional since nothing is blocking here)
Base.fetch(::AverageDistanceTemperatureRecord) = nothing
```

To perform a simulation that runs for 1000 cycles with logarithmically increasing
temperatures, no trajectory stored but a record the average weighted distance to point
`(0.5, 0.9, 1.3)` (in Å) and energies stored every 2 cycles, use
`SimulationSetup((j,n)->log1p(j/n)*300u"K"/log(2), 1000, "", 2, AverageDistanceTemperatureRecord([0.5, 0.9, 1.3]u"Å"))`
"""
Base.@kwdef struct SimulationSetup{Trecord}
    temperatures::Vector{TK}
    pressure::typeof(1.0u"Pa")=-Inf*u"Pa"
    ncycles::Int
    ninit::Int=0
    outdir::String=""
    printevery::Int=0
    outtype::Vector{Symbol}=[:energies, :zst]
    record::Trecord=Returns(nothing)

    function SimulationSetup(T, P::Unitful.Pressure, ncycles::Int, ninit::Int=0, outdir::String="", printevery::Int=0, outtype::Vector{Symbol}=[:energies, :zst], record=Returns(nothing))
        @assert ncycles ≥ 0 && ninit ≥ 0
        n = ncycles + ninit
        temperatures = if T isa Number
            fill(T, n + (n == 0))
        elseif T isa Vector
            convert(Vector{TK}, T)
        else
            n ≤ 1 ? [T(1,2)] : T.(1:n, n)
        end
        known_outtype = Set((:energies, :pdb, :stream, :zst))
        union!(known_outtype, [Symbol(:initial_, x) for x in known_outtype])
        unknown_outtype = setdiff(outtype, known_outtype)
        isempty(unknown_outtype) || error(lazy"Unknown outtype: $(join(unknown_outtype, \", \", \" and \"))")
        ret = new{typeof(record)}(temperatures, P, ncycles, ninit, outdir, printevery, outtype, record)
        initialize_record!(record, ret)
        ret
    end
end
function SimulationSetup(T, ncycles::Int, ninit::Int=0, outdir::String="", printevery::Int=0, outtype::Vector{Symbol}=[:energies, :zst], record=Returns(nothing))
    SimulationSetup(T, -Inf*u"Pa", ncycles, ninit, outdir, printevery, outtype, record)
end
function SimulationSetup(T, ncycles::Integer; kwargs...)
    SimulationSetup(; temperatures=T, ncycles, kwargs...)
end

initialize_record!(::T, ::SimulationSetup{T}) where {T} = nothing
needcomplete(::Any) = true
needcomplete(::Returns) = false

mutable struct MoveKind
    trials::Int
    accepted::Int
end
MoveKind() = MoveKind(0, 0)
attempt!(x::MoveKind) = (x.trials += 1)
accept!(x::MoveKind) = (x.accepted += 1)
(x::MoveKind)() = x.accepted / x.trials

macro complete_mcmoves(def)
    append!(def.args[end].args, Expr(:const, Expr(:(::), f, :MoveKind)) for f in mcmovenames)

    arg_constructor = Expr(:call, :MoveStatistics, :dmax, :θmax, Expr(:call, :MoveKind))
    append!(arg_constructor.args, Expr(:call, :MoveKind) for _ in mcmovenames)
    constructor = Expr(:function,
        Expr(:call, :MoveStatistics, :dmax, :θmax),
        arg_constructor
    )

    arg_attempt = Expr(:if)
    expr = arg_attempt
    for (imove, move) in enumerate(mcmovenames)
        if imove != length(mcmovenames)
            push!(expr.args, Expr(:call, :(===), :move, QuoteNode(move)))
        end
        push!(expr.args, Expr(:call, :attempt!, Expr(:., :statistics, QuoteNode(move))))
        if imove == length(mcmovenames) - 1
            push!(expr.args, Expr(:block))
        elseif imove < length(mcmovenames)
            push!(expr.args, Expr(:elseif))
        end
        expr = last(expr.args)::Expr
    end
    attempt! = Expr(:function,
        Expr(:call, :attempt!, :statistics, :move),
        Expr(:block, arg_attempt, Expr(:call, :attempt!, Expr(:., :statistics, QuoteNode(:total))))
    )

    accept! = deepcopy(attempt!)
    Q = Expr[accept!]
    for ex in Q
        for (i, arg) in enumerate(ex.args)
            if arg === :attempt!
                ex.args[i] = :accept!
            elseif arg isa Expr
                push!(Q, arg)
            end
        end
    end

    esc(Expr(:block, def, constructor, attempt!, accept!))
end

@complete_mcmoves mutable struct MoveStatistics
    dmax::TÅ # maximum distance for a translation move
    θmax::typeof(1.0u"°") # maximum angle for a rotation move
    const total::MoveKind
    # const mcmove::Movekind ... are generated by the macro
end
# generated by the macro:
# function MoveStatistics(dmax, θmax)
#     MoveStatistics(dmax, θmax, MoveKind()...)
# end
# function attempt!(statistics, move)
#     if move === :translation
#         attempt!(statistics.translation)
#     elseif move === :rotation
#         ...
#     end
#     attempt!(statistics.total)
# end
# function accept!(statistics, move)
#     ...
# end


"""
    choose_newpos!(statistics::MoveStatistics, mc::MonteCarloSetup, pos::Vector{SVector{3,TÅ}}, randomdmax::TÅ, i::Int)

Choose an MC move for the molecule at positions given by `pos` and return a pair
`(newpos, movekind)` where `newpos` is the list of atom positions after the move, and
`movekind` is the kind of move that was chosen, as a `Symbol`.

`i` is the kind of the molecule

`randomdmax` is the maximum distance between two points in the unit cell (used in random
translations for instance).
"""
function choose_newpos!(statistics::MoveStatistics, mc::MonteCarloSetup,
                        pos::AbstractVector{SVector{3,TÅ}}, randomdmax::TÅ, i::Int, ffidxi,
                        movekind=mc.mcmoves[i](rand(mc.rng)))
    attempt!(statistics, movekind)
    if movekind === :swap
        movekind = ifelse(rand(mc.rng, Bool), :swap_deletion, :swap_insertion)
    end
    speciesblock = mc.speciesblocks[i]
    for _ in 1:100
        newpos = if movekind === :translation
            random_translation(mc.rng, pos, statistics.dmax)
        elseif movekind === :rotation
            random_rotation(mc.rng, pos, statistics.θmax, mc.bead[i])
        elseif movekind === :random_translation
            random_translation(mc.rng, pos, randomdmax)
        elseif movekind === :random_rotation
            random_rotation(mc.rng, pos, 180u"°", mc.bead[i])
        elseif movekind === :random_reinsertion || movekind === :swap_insertion
            random_rotation(mc.rng, random_translation(mc.rng, pos, randomdmax), 180u"°", mc.bead[i])
        elseif movekind === :swap_deletion
            return Vector(pos), :swap_deletion, false
        else
            error(lazy"Unknown move kind: $movekind")
        end
        blocked = inblockpocket(speciesblock, mc.atomblocks, ffidxi, newpos)
        if movekind === :swap_insertion
            # to be consistent with the excluded volume, only refuse insertion where the
            # bead is in the species blockpocket
            bead = mc.bead[i]
            bead += bead==0
            if !speciesblock[newpos[bead]]
                return newpos, movekind, blocked
            end
        elseif !blocked
            return newpos, movekind, false
        end
    end
    movekind === :swap_insertion || @warn "Trapped species did not manage to move out of a blocked situation. This could be caused by an impossible initial configuration."
    return Vector(pos), movekind, true # signal that the move was blocked
end


function try_swap!(mc::MonteCarloSetup, i::Int, statistics::MoveStatistics, randomdmax, temperature, energy, φPV_div_k::TK)
    ffidxi = mc.step.ffidx[i]
    pos, move, blocked = choose_newpos!(statistics, mc, mc.models[i], randomdmax, i, ffidxi, :swap)
    blocked && return energy # insertion failed due to blocking spheres

    isinsertion = move === :swap_insertion
    j = 0
    posidxi = mc.step.posidx[i]
    newpos = if isinsertion
        j = length(posidxi) + 1
        pos
    else
        isempty(posidxi) && return energy # deletion failed because there is no molecule to delete
        j = rand(1:length(posidxi))
        molpos = posidxi[j]
        mc.step.positions[molpos]
    end
    idx = (i, j)

    contrib = movement_energy(mc, idx, newpos)
    ediff = isdefined_ewald(mc) ? mc.ewald.ctx.energies[i] : 0.0u"K"
    if !isinsertion
        contrib = -contrib
        ediff = -ediff
    end

    swapinfo = SwapInformation(φPV_div_k, i, true, isinsertion)

    newenergy = first(handle_acceptation(mc, (i,j), MCEnergyReport(0.0u"K", 0.0u"K", 0.0u"K", ediff), contrib, temperature, newpos, move, statistics, nothing, energy, swapinfo))

    return newenergy
end

function risk_of_underflow(current::MCEnergyReport, diff::MCEnergyReport)
    vc = (current.framework.direct, current.framework.vdw, current.inter, current.reciprocal)
    vd = (diff.framework.direct, diff.framework.vdw, diff.inter, diff.reciprocal)
    for i in 1:4
        xc, xd = vc[i], vd[i]
        abs(xc) + abs(xd) > 1e-15u"K" && isapprox(xc, -xd; rtol=0.01) && return true
    end
    return false
end

function underflow_averted_warning(energy::MCEnergyReport, newenergy::MCEnergyReport, diff::MCEnergyReport)
    ve = (energy.framework.direct, energy.framework.vdw, energy.inter, energy.reciprocal)
    vn = (newenergy.framework.direct, newenergy.framework.vdw, newenergy.inter, newenergy.reciprocal)
    vd = (diff.framework.direct, diff.framework.vdw, diff.inter, diff.reciprocal)
    for i in 1:4
        xe, xn, xd = ve[i], vn[i], vd[i]
        abs(xe) + abs(xn) > 1e-15u"K" && !isapprox(xe + xd, xn; rtol=1e-8) && !isapprox(xn - xd, xe; rtol=1e-8) && (#=(@show i, xe, xn, xd);=# return true)
    end
    return false
end

function handle_acceptation(mc::MonteCarloSetup, idx, before, after, temperature, oldpos, move::Symbol, statistics::MoveStatistics, running_update, energy, swapinfo::SwapInformation)
    accepted = compute_accept_move(before, after, temperature, mc, swapinfo)
    if accepted
        parallel = !(running_update isa Nothing) && mc.step.parallel
        if parallel && !swapinfo.isswap # swap moves need to be updated before choosing the next idx
            # do not use newpos since it can be changed in the next iteration before the Task is run
            running_update = @spawn update_mc!($mc, $idx, $oldpos, $swapinfo)
        else
            update_mc!(mc, idx, oldpos, swapinfo)
        end
        diff = after - before
        accept!(statistics, ifelse(swapinfo.isswap, :swap, move))
        if risk_of_underflow(energy.er, diff)
            parallel && wait(running_update)
            newenergy = baseline_energy(mc) # reset computation to avoid underflows
            if underflow_averted_warning(energy.er, newenergy.er, diff)
                @error "Underflow mitigation stems from incorrect energy computation:" energy diff newenergy
            end
            newenergy
        elseif swapinfo.isswap # update the tail correction as well
            modify_species!(mc.tailcorrection, idx[1], ifelse(swapinfo.isinsertion, 1, -1))
            BaselineEnergyReport(energy.er + diff, mc.tailcorrection[])
        else
            energy + diff
        end
    else
        energy
    end, accepted, running_update
end

function print_report(simulation::SimulationSetup, time_begin, time_end, statistics::MoveStatistics)
    open(joinpath(simulation.outdir, "report.txt"), "a") do io
        if position(io) > 0
            println(io, "\n\n\n### NEW REPORT STARTING HERE ###\n")
        end
        println(io, "Simulation ran for ", time_end-time_begin, " s")
        print(io, "Number of cycles: ", simulation.ncycles)
        if simulation.ninit > 0
            println(io, " (+ ", simulation.ninit, " initialization cycles)")
        else
            println(io)
        end
        if allequal(simulation.temperatures)
            println(io, "Temperature: ", simulation.temperatures[1])
        else
            print(io, "Temperatures: ", simulation.temperatures[1], " (initial) – ")
            if simulation.ninit > 0
                print(io, simulation.temperatures[1+simulation.ninit], " (start of production) – ")
            end
            print(io, simulation.temperatures[end], " (final)")
        end
        println(io)

        println(io, "Accepted moves: ", statistics.total.accepted, '/', statistics.total.trials, " (attempted) = ", statistics.total(), " among which:")
        for (i, symb) in enumerate(mcmovenames)
            stat = getfield(statistics, 3+i)::MoveKind
            stat.trials == 0 && continue
            str = replace(String(symb), '_'=>' ')
            println(io, " - accepted ", str, stat.accepted > 1 ? "s: " : ": ", stat.accepted, '/', stat.trials, " (attempted) = ", stat())
        end
        println(io, "Final dmax: ", statistics.dmax)
        println(io, "Final θmax: ", statistics.θmax)
        println(io)
    end
end


"""
    run_montecarlo!(mc::MonteCarloSetup, simu::SimulationSetup)

Run a Monte-Carlo simulation. Return the list of energies of the system during the
simulation.

See [`MonteCarloSetup`](@ref) for the definition of the system and
[`SimulationSetup`](@ref) for the parameters of the simulation.
"""
function run_montecarlo!(mc::MonteCarloSetup, simu::SimulationSetup)
    time_begin = time()
    # energy initialization
    energy = baseline_energy(mc)
    if isinf(Float64(energy)) || isnan(Float64(energy))
        @error "Initial energy is not finite, this probably indicates a problem with the initial configuration."
    end
    energies = typeof(energy)[]
    initial_energies = typeof(energy)[]

    # record and outputs
    mkpath(simu.outdir)
    startstep = SimulationStep(mc.step, needcomplete(simu.record) ? :all : :output)
    output, output_task = output_handler(simu.outdir, simu.outtype, startstep, false)

    has_initial_output = any(startswith("initial")∘String, simu.outtype)
    initial_output, initial_output_task = output_handler(simu.outdir, simu.outtype, startstep, true)

    parallel = mc.step.parallel
    if parallel
        record_task = @spawn $simu.record(startstep, $energy, -simu.ninit, $mc, $simu)
    else
        simu.record(startstep, energy, -simu.ninit, mc, simu) === :stop && @goto end_cleanup
        record_task = Task(Returns(nothing)) # for type stability
    end

    # value initialisations
    old_idx = (0,0)
    accepted = false
    statistics = MoveStatistics(1.3u"Å", 30.0u"°")
    randomdmax = maximum(norm, eachcol(mc.step.mat))

    if simu.ninit == 0
        put!(output, startstep)
        push!(energies, energy)
        if simu.ncycles == 0 # single-point computation
            @goto end_cleanup
        end
    else
        put!(initial_output, startstep)
        push!(initial_energies, energy)
    end

    # GCMC handling
    has_swap = [i for (i, mcmoves) in enumerate(mc.mcmoves) if mcmoves[:swap] > 0]
    if !isempty(has_swap)
        allequal(simu.temperatures) || error("Cannot currently run GCMC with variation of temperature")
        simu.pressure == -Inf*u"Pa" && error("Cannot perform GCMC with unspecified pressure")
    end
    φPV_div_k = Vector{typeof(1.0u"K")}(undef, length(mc.mcmoves))
    for i_swap in has_swap
        T0 = first(simu.temperatures)
        P = simu.pressure
        φ = only(Clapeyron.fugacity_coefficient(mc.gcmcdata.model[i_swap], P, T0; phase=:stable, vol0=Clapeyron.volume(mc.gcmcdata.model0[i_swap], P, T0)))
        isnan(φ) && error("Specified gas not in gas form at the required temperature ($T0) and pressure ($P)!")
        PV_div_k = uconvert(u"K", P*mc.gcmcdata.volumes[i_swap]/u"k")
        φPV_div_k[i_swap] = φ*PV_div_k
    end

    # idle initialisations
    running_update = @spawn nothing
    local oldpos::Vector{SVector{3,TÅ}}
    local before::MCEnergyReport
    local after::MCEnergyReport

    # main loop
    for (counter_cycle, idx_cycle) in enumerate((-simu.ninit+1):simu.ncycles)
        temperature = simu.temperatures[counter_cycle]

        for i_swap in has_swap
            accepted && parallel && wait(running_update)
            for _ in 1:4
                energy = try_swap!(mc, i_swap, statistics, randomdmax, temperature, energy, φPV_div_k[i_swap])
            end
            old_idx = (0,0)
        end

        numsteps = max(20, length(mc.revflatidx))
        for idx_step in 1:numsteps
            isempty(mc.revflatidx) && break

            # choose the species on which to attempt a move
            idx = choose_random_species(mc)
            ffidxi = mc.step.ffidx[idx[1]]

            # currentposition is the position of that species
            # currentposition is either a Vector or a @view, that's OK
            currentposition = (accepted&(old_idx==idx)) ? oldpos : @view mc.step.positions[mc.step.posidx[idx[1]][idx[2]]]

            # newpos is the position after the trial move
            newpos, move, blocked = choose_newpos!(statistics, mc, currentposition, randomdmax, idx[1], ffidxi)

            # If, despite the multiple attempts, the move is blocked, skip early to the
            # next iteration.
            # No need to update oldpos, old_idx nor accepted thus.
            # The refused move is still taken into account in the statistics.
            blocked && continue

            isswap = move === :swap_insertion || move === :swap_deletion
            if move === :swap_insertion
                idx = (idx[1], length(mc.step.posidx[idx[1]]) + 1)
            end

            # if the previous move was accepted, wait for mc to be completely up to date
            # before computing
            accepted && parallel && wait(running_update)


            # Core computation: energy difference between after and before the move
            if old_idx == idx || isswap
                if accepted
                    before = after
                # else
                #     before = fetch(before_task) # simply keep its previous value
                end
                if !(old_idx == idx && ((move === :swap_deletion) & accepted)) # no need to compute it again in that specific case
                    after = movement_energy(mc, idx, newpos)
                end
                if move === :swap_deletion
                    after = -after
                    before = MCEnergyReport(0.0u"K", 0.0u"K", 0.0u"K", isdefined_ewald(mc) ? -mc.ewald.ctx.energies[idx[1]] : 0.0u"K")
                elseif move === :swap_insertion
                    before = MCEnergyReport(0.0u"K", 0.0u"K", 0.0u"K", isdefined_ewald(mc) ? mc.ewald.ctx.energies[idx[1]] : 0.0u"K")
                end
            else
                before, after = combined_movement_energy(mc, idx, newpos)
            end

            swapinfo = SwapInformation(φPV_div_k[idx[1]], idx[1], isswap, move === :swap_insertion)

            oldpos = newpos
            energy, accepted, running_update = handle_acceptation(mc, idx, before, after, temperature, oldpos, move, statistics, running_update, energy, swapinfo)
            old_idx = ifelse(move === :swap_deletion, 0, idx)
        end

        # end of cycle
        report_now = (idx_cycle ≥ 0 && (idx_cycle == 0 || (simu.printevery > 0 && idx_cycle%simu.printevery == 0))) ||
                     (idx_cycle < 0 && has_initial_output && (simu.printevery == 0 || (simu.printevery > 0 && idx_cycle%simu.printevery == 0)))
        if !(simu.record isa Returns) || report_now
            accepted && parallel && wait(running_update)
            if idx_cycle == 0 # start production with a precise energy
                energy = baseline_energy(mc)
            end
            o = SimulationStep(mc.step, :output)
            if report_now
                if idx_cycle < 0
                    put!(initial_output, o)
                    push!(initial_energies, energy)
                else
                    put!(output, o)
                    push!(energies, energy)
                end
            end
            if !(simu.record isa Returns)
                ocomplete = needcomplete(simu.record) ? SimulationStep(o, :complete_output) : o
                if parallel
                    # ret_record cannot be inferred, but that's fine
                    ret_record = fetch(record_task)
                    ret_record === :stop && break
                    record_task = @spawn $simu.record($ocomplete, $energy, $idx_cycle, $mc, $simu)
                else
                    simu.record(ocomplete, energy, idx_cycle, mc, simu)
                end
            end
            yield()
        end

        translation_ratio = statistics.translation()
        if !isnan(translation_ratio)
            statistics.dmax = clamp(statistics.dmax * (1 + (translation_ratio - 0.5)*sqrt(10/(99+counter_cycle))), 0.1u"Å", 3.0u"Å")
        end
        rotation_ratio = statistics.rotation()
        if !isnan(rotation_ratio)
            statistics.θmax = (counter_cycle-1)/counter_cycle*statistics.θmax + rotation_ratio/counter_cycle*120u"°"
        end
    end
    if parallel
        wait(record_task)
        accepted && wait(running_update)
    end

    @label end_cleanup
    if simu.printevery == 0 && simu.ncycles > 0
        put!(output, mc.step)
        push!(energies, energy)
    end
    ozero = SimulationStep(mc.step, :zero)
    put!(output, ozero)
    put!(initial_output, ozero)
    if !isempty(simu.outdir)
        for initial_ in (:initial_, Symbol(""))
            if Symbol(initial_, :energies) in simu.outtype
                serialize(joinpath(simu.outdir, string(initial_, "energies.serial")), initial_ == :initial_ ? initial_energies : energies)
            end
        end
    end
    if !(simu.ninit == 0 && simu.ncycles == 0) # not a single-point computation
        lastenergy = baseline_energy(mc)
        if !isapprox(Float64(energy), Float64(lastenergy), rtol=1e-9)
            @error "Energy deviation observed between actual ($lastenergy) and recorded ($energy), this means that the simulation results are wrong!" actual=lastenergy recorded=energy
        end
    end
    parallel && (wait(output_task[]); wait(initial_output_task[]))
    fetch(simu.record)
    time_end = time()
    print_report(simu, time_begin, time_end, statistics)
    energies
end
