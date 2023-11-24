export SimulationSetup

"""
    SimulationSetup{Trecord}

Definition of the simulation setup, which must provide the following information:
- `T` the temperature. Either a number (in K), a list of temperatures (in K, one per cycle)
  or a function. If a function, it must take two parameters: `j` the number of the current
  cycle and `n` the total number of cycles, and return a temperature (in K).
  Internally, the list of temperatures is stored as a field `temperatures` indexed by `j`,
  which varies between `1` and `n`.
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
- `outtype`: a `Vector{Symbol}` whose elements represent the different outputs that
  should be produced in `outdir`. The elements can be:
  + `:energies` to have the list of energies output as a serialized
    `Vector{BaselineEnergyReport}` in a "energies.serial" file.
  + `:pdb` to have the positions of the atoms output as a "trajectory.pdb" file.
  + `:stream` to have the positions of the atoms output as a [`StreamSimulationStep`](@ref)
    in a "steps.stream" file.
  + `:zst` to have the positions of the atoms output as a compressed
    [`StreamSimulationStep`](@ref) in a "steps.zst" file.
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
    ncycles::Int
    ninit::Int=0
    outdir::String=""
    printevery::Int=0
    outtype::Vector{Symbol}=[:energies, :zst]
    record::Trecord=Returns(nothing)

    function SimulationSetup(T, ncycles::Int, ninit::Int=0, outdir::String="", printevery::Int=0, outtype::Vector{Symbol}=[:energies, :zst], record=Returns(nothing))
        @assert ncycles ≥ 0 && ninit ≥ 0
        n = ncycles + ninit
        temperatures = if T isa Number
            fill(T, n + (n == 0))
        elseif T isa Vector
            convert(Vector{TK}, T)
        else
            n ≤ 1 ? [T(1,2)] : T.(1:n, n)
        end
        unknownouttype = setdiff(outtype, Set((:energies, :pdb, :stream, :zst)))
        isempty(unknownouttype) || error(lazy"Unknown outtype: $(join(unknownouttype, \", \", \" and \"))")
        ret = new{typeof(record)}(temperatures, ncycles, ninit, outdir, printevery, outtype, record)
        initialize_record!(record, ret)
        ret
    end
end
function SimulationSetup(T, ncycles; kwargs...)
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
                        pos::AbstractVector{SVector{3,TÅ}}, randomdmax::TÅ, i::Int, ffidxi)
    r = rand(mc.rng)
    movekind = mc.mcmoves[i](r)
    attempt!(statistics, movekind)
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
        elseif movekind === :random_reinsertion
            random_rotation(mc.rng, random_translation(mc.rng, pos, randomdmax), 180u"°", mc.bead[i])
        else
            error(lazy"Unknown move kind: $movekind")
        end
        if !inblockpocket(speciesblock, mc.atomblocks, ffidxi, newpos)
            return newpos, movekind, false
        end
    end
    @warn "Trapped species did not manage to move out of a blocked situation. This could be caused by an impossible initial configuration."
    return pos, movekind, true # signal that the move was blocked
end


function print_report(outdir, time_begin, time_end, statistics)
    open(joinpath(outdir, "report.txt"), "a") do io
        if position(io) > 0
            println(io, "\n\n\n### NEW REPORT STARTING HERE ###\n")
        end
        println(io, "Simulation ran for ", time_end-time_begin, " s")
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
    reports = typeof(energy)[]

    # record and outputs
    mkpath(simu.outdir)
    startstep = SimulationStep(mc.step, needcomplete(simu.record) ? :all : :output)
    output, output_task = output_handler(simu.outdir, simu.outtype, startstep)
    parallel = mc.step.parallel
    if parallel
        record_task = @spawn $simu.record(startstep, $energy, -simu.ninit, $mc, $simu)
    else
        simu.record(startstep, energy, -simu.ninit, mc, simu) === :stop && @goto end_cleanup
        record_task = Task(Returns(nothing)) # for type stability
    end

    if simu.ninit == 0
        put!(output, startstep)
        push!(reports, energy)
        if simu.ncycles == 0 # single-point computation
            @goto end_cleanup
        end
    end

    # idle initialisations
    running_update = @spawn nothing
    local oldpos::Vector{SVector{3,TÅ}}
    local before::MCEnergyReport
    local after::MCEnergyReport

    # value initialisations
    nummol = max(20, length(mc.indices))
    old_idx = (0,0)
    accepted = false
    statistics = MoveStatistics(1.3u"Å", 30.0u"°")
    randomdmax = maximum(norm, eachcol(mc.step.mat))

    # main loop
    for (counter_cycle, idx_cycle) in enumerate((-simu.ninit+1):simu.ncycles)
        temperature = simu.temperatures[counter_cycle]

        for idnummol in 1:nummol
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


            # if the previous move was accepted, wait for mc to be completely up to date
            # before computing
            accepted && parallel && wait(running_update)

            i, j = idx
            ij = mc.offsets[i]+j

            if old_idx == idx
                if accepted
                    before = after
                # else
                #     before = fetch(before_task) # simply keep its previous value
                end
                @spawnif parallel begin
                    singlereciprocal = @spawn single_contribution_ewald($mc.ewald, $ij, $newpos)
                    fer = @spawn framework_interactions($mc, $i, $newpos)
                    singlevdw = single_contribution_vdw(mc.step, idx, newpos)
                    after = MCEnergyReport(fetch(fer)::FrameworkEnergyReport, singlevdw, fetch(singlereciprocal)::TK)
                end
            else
                molpos = mc.step.posidx[i][j]
                positions = @view mc.step.positions[molpos]
                @spawnif parallel begin
                    singlevdw_before_task = @spawn single_contribution_vdw($mc.step, $(i,j), $positions)
                    singlereciprocal_after = @spawn single_contribution_ewald($mc.ewald, $ij, $newpos)
                    singlereciprocal_before = @spawn single_contribution_ewald($mc.ewald, $ij, nothing)
                    fer_before = @spawn framework_interactions($mc, $i, $positions)
                    fer_after = @spawn framework_interactions($mc, $i, $newpos)
                    singlevdw_before = fetch(singlevdw_before_task)::TK
                    singlevdw_after = single_contribution_vdw(mc.step, (i,j), newpos)
                    before = MCEnergyReport(fetch(fer_before)::FrameworkEnergyReport, singlevdw_before, fetch(singlereciprocal_before)::TK)
                    after = MCEnergyReport(fetch(fer_after)::FrameworkEnergyReport, fetch(singlevdw_after)::TK, fetch(singlereciprocal_after)::TK)
                end
            end

            old_idx = idx
            accepted = compute_accept_move(before, after, temperature, mc.rng)
            oldpos = newpos
            if accepted
                if parallel
                    # do not use newpos since it can be changed in the next iteration before the Task is run
                    running_update = @spawn update_mc!($mc, $idx, $oldpos)
                else
                    update_mc!(mc, idx, oldpos)
                end
                diff = after - before
                accept!(statistics, move)
                if abs(Float64(diff.framework)) > 1e50 # an atom left a blocked pocket
                    parallel && wait(running_update)
                    energy = baseline_energy(mc) # to avoid underflows
                else
                    energy += diff
                end
            end
        end

        # end of cycle
        report_now = idx_cycle ≥ 0 && (idx_cycle == 0 || (simu.printevery > 0 && idx_cycle%simu.printevery == 0))
        if !(simu.record isa Returns) || report_now
            accepted && parallel && wait(running_update)
            o = SimulationStep(mc.step, :output)
            if report_now
                put!(output, o)
                push!(reports, energy)
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
            statistics.dmax *= 1 + (translation_ratio - 0.5)*(1.0 + 3*translation_ratio)*100/(99+counter_cycle)
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
        push!(reports, energy)
    end
    put!(output, SimulationStep(mc.step, :zero))
    if !isempty(simu.outdir) && :energies in simu.outtype
        serialize(joinpath(simu.outdir, "energies.serial"), reports)
    end
    if !(simu.ninit == 0 && simu.ncycles == 0) # not a single-point computation
        lastenergy = baseline_energy(mc)
        if !isapprox(Float64(energy), Float64(lastenergy), rtol=1e-9)
            @error "Energy deviation observed between actual ($lastenergy) and recorded ($energy), this means that the simulation results are wrong!"
        end
    end
    parallel && wait(output_task[])
    fetch(simu.record)
    time_end = time()
    print_report(simu.outdir, time_begin, time_end, statistics)
    reports
end
