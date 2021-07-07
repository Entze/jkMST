using Logging
using Crayons.Box
using LightGraphs, SimpleWeightedGraphs
using JuMP
using GLPK, SCIP, CPLEX
using Random
using ProgressMeter
using PrettyTables

include("common.jl")
include("heuristic.jl")
include("util.jl")
include("mtz.jl")
include("scf.jl")
include("mcf.jl")
include("cec.jl")

@enum HeaderName Filename Graphsize Treesize Opt Sol

function colot(name::HeaderName;
                        f::Union{Int, Nothing}=nothing,
                        ks::Union{Int, Nothing}=nothing,
                        k::Union{Int, Nothing}=nothing,
                        ms::Union{Int, Nothing}=nothing,
                        m::Union{Int, Nothing}=nothing,
                        sols::Union{Int, Nothing}=nothing,
                        sol::Union{Int, Nothing}=nothing,
                        line::Union{Int, Nothing}=nothing)
    return columnoftable(name, f=f, ks=ks, k=k, ms=ms, m=m, sols=sols, sol=sol, line=line)
end

function columnoftable(name::HeaderName;
                        f::Union{Int, Nothing}=nothing,
                        ks::Union{Int, Nothing}=nothing,
                        k::Union{Int, Nothing}=nothing,
                        ms::Union{Int, Nothing}=nothing,
                        m::Union{Int, Nothing}=nothing,
                        sols::Union{Int, Nothing}=nothing,
                        sol::Union{Int, Nothing}=nothing,
                        line::Union{Int, Nothing}=nothing)
    if name == Filename
        return 1
    elseif name == Graphsize
        return 2
    elseif name == Treesize
        return 3
    elseif name == Opt
        return 4
    elseif name == Sol
        @assert !isnothing(sols) "sols is not defined."
        @assert !isnothing(sol) "sol is not defined."
        @assert !isnothing(ms) "ms is not defined."
        @assert !isnothing(m) "m is not defined."
        return 4 + Int(sols != 1) + (Int(sols != 1) + ms) * (sol - 1) + m
    else
        @assert false "Defensive assert: Illegal case."
    end
end

function rowot(name::HeaderName;
                    f::Union{Int, Nothing}=nothing,
                    ks::Union{Int, Nothing}=nothing,
                    k::Union{Int, Nothing}=nothing,
                    ms::Union{Int, Nothing}=nothing,
                    m::Union{Int, Nothing}=nothing,
                    sols::Union{Int, Nothing}=nothing,
                    sol::Union{Int, Nothing}=nothing,
                    line::Union{Int, Nothing}=nothing)
    return rowoftable(name, f=f, ks=ks, k=k, ms=ms, m=m, sols=sols, sol=sol, line=line)
end

function rowoftable(name::HeaderName;
                    f::Union{Int, Nothing}=nothing,
                    ks::Union{Int, Nothing}=nothing,
                    k::Union{Int, Nothing}=nothing,
                    ms::Union{Int, Nothing}=nothing,
                    m::Union{Int, Nothing}=nothing,
                    sols::Union{Int, Nothing}=nothing,
                    sol::Union{Int, Nothing}=nothing,
                    line::Union{Int, Nothing}=nothing)
    lines = 4
    if name == Filename || name == Graphsize
        @assert !isnothing(f) "f is not defined."
        @assert !isnothing(ks) "ks is not defined."
        return 1 + (f - 1) * ks * lines
    elseif name == Treesize || name == Opt
        @assert !isnothing(f) "f is not defined."
        @assert !isnothing(ks) "ks is not defined."
        @assert !isnothing(k) "k is not defined."
        return 1 + (f - 1) * ks * lines + (k - 1) * lines
    elseif name == Sol
        @assert !isnothing(f) "f is not defined."
        @assert !isnothing(ks) "ks is not defined."
        @assert !isnothing(k) "k is not defined."
        @assert !isnothing(line) "line is not defined."
        return 1 + (f - 1) * ks * lines + (k - 1) * lines + (line - 1)
    else
        @assert false "Defensive assert: Illegal case."
    end
end

function solver_from_string(x::AbstractString)
    i::Int = -1
    xl = lowercase(x)
    for solver in instances(Solver)
        if string(solver) == xl
            i = Int(solver)
            break
        end
    end
    if i == -1
        error("No solver with name $x")
    end
    return Solver(i)
end

function solvingmode_from_string(x::AbstractString)
    i::Int = -1
    for mode in instances(SolvingMode)
        if string(mode) == x
            i = Int(mode)
            break
        end
    end
    if i == -1
        error("No mode with name $x")
    end
    return SolvingMode(i)
end

function main(;files::Vector{String}=[],
                directories::Vector{String}=[],
                modestrings::Vector{String}=["mtz"],
                kstrings::Vector{String}=["1.0"],
                solverstrings::Vector{String}=["cplex"],
                timeout_sec::Float64=Inf64,
                generate_timeout_sec::Float64=Inf64,
                preprocessinstances::Bool=true,
                preprocesstimeout::Float64=Inf64,
                printsolutiongraphs::Bool=true,
                printintermediatetable::Bool=true,
                debugmodels::Bool=true,
                intermediatetableinterval::Float64=60.0,
                tablestyle::String="text"
                )
    fl::Int = length(files)
    dl::Int = length(directories)
    modes::Vector{SolvingMode} = map(solvingmode_from_string, modestrings)
    ml::Int = length(modes)
    ks::Vector{Union{Float64,Int}} = map(s -> begin
        v::Float64 = parse(Float64, s)
        if v > 1 && v == trunc(v)
        return Int(trunc(Int, v))
    end
        if !(0 <= v <= 1)
        error("Specified illegal k=$s. Should be Float ∈ [0,1] or Int if > 1.")
    end
        return v
    end, kstrings)
    kl::Int = length(ks)
    solvers::Vector{Solver} = map(solver_from_string, solverstrings)
    sl::Int = length(solvers)
    @debug "Got $(fl) file$(fl == 1 ? "" : "s"), $(dl) $(dl == 1 ? "directory" : "directories"), $(ml) mode$(ml == 1 ? "" : "s"), $(kl) size$(kl == 1 ? "" : "s"), $(sl) solver$(sl == 1 ? "" : "s"), $(timeout == Inf64 ? "no" : format_seconds_readable(timeout)) timeout."
    fl = enumeratefiles!(files, directories, fl=fl)
    if fl == 0
        return
    end
    @assert fl > 0
    tablestyle = lowercase(tablestyle)
    if tablestyle == "text"
        tablebackend = :text
    elseif tablestyle == "html"
        tablebackend = :html
    elseif tablestyle == "latex"
        tablebackend = :latex
    else
        error("Unknown tablestyle: $tablestyle")
    end
    instances::Int = fl * kl * ml * sl
    if timeout != Inf64
        @debug "$instances instance$(instances == 1 ? "" : "s") to solve, $(format_seconds_readable(instances * timeout)) runtime upperbound."
    end
    (header, data) = generatetable(files=files, fl=fl, solverstrings=solverstrings, sl=sl, modestrings=modestrings, ml=ml, ks=ks, kl=kl)
    kMSTs!(files=files,
            fl=fl,
            sizes=ks,
            kl=kl,
            solvers=solvers,
            sl=sl,
            modes=modes,
            ml=ml,
            timeout_sec=timeout_sec,
            generate_timeout_sec=generate_timeout_sec,
            header=header,
            data=data,
            preprocessinstances=preprocessinstances,
            preprocesstimeout=preprocesstimeout,
            printintermediatetable=printintermediatetable,
            intermediatetableinterval=intermediatetableinterval,
            printsolutiongraphs=printsolutiongraphs,
            tablebackend=tablebackend)
    return
end

function enumeratefiles!(files::Vector{String}=[], directories::Vector{String}=[];fl::Int=length(files))::Int
    enumfiles = ProgressUnknown(desc="Enumerating files:", dt=0.5, enabled=isdebug())
    ProgressMeter.update!(enumfiles, fl)
    while !isempty(directories)
        directory = pop!(directories)
        for (root, dirs, fs) in walkdir(directory)
            for dir in dirs
                push!(directories, dir)
            end
            for file in fs
                if endswith(file, ".dat")
                    push!(files, joinpath(root, file))
                    ProgressMeter.next!(enumfiles)
                end
    end
        end
    end
    ProgressMeter.finish!(enumfiles)
    fl = length(files)
    @debug "Found $fl file$(fl == 1 ? "" : "s")."
    return fl
end

function generatetable(;files::Vector{String}=[],
                        fl::Int=length(files),
                        solverstrings::Vector{String}=[],
                        sl::Int=length(solverstrings),
                        modestrings::Vector{String}=[],
                        ml::Int=length(modestrings),
                        ks::Vector{Union{Int, Float64}} = [],
                        kl::Int=length(ks))::Tuple{Vector{String},Matrix{Union{Int, Float64,String}}}
    header = ["Graph", "|V|", "k", "OPT"]
    if sl > 1
        for s in solverstrings
            push!(header, s)
            append!(header, modestrings)
        end
    else
        append!(header, modestrings)
    end
    lines = 4
    data::Matrix{Union{String, Int, Float64}} = fill("", fl * kl * lines, length(header))
    for (f, file) in enumerate(files)
        data[rowot(Filename, f=f, ks=kl), colot(Filename)] = file
        for (k, size) in enumerate(ks)
            data[rowot(Treesize, f=f, k=k, ks=kl), colot(Treesize)] = size
        end
    end
    return (header, data)
end

function kMSTs!(;
    files::Vector{String} = [],
    fl::Int=length(files),
    sizes::Vector{Union{Float64, Int}} = [],
    kl::Int=length(sizes),
    solvers::Vector{Solver} = [],
    sl::Int=length(solvers),
    modes::Vector{SolvingMode} = [],
    ml::Int=length(modes),
    timeout_sec::Float64 = Inf64,
    generate_timeout_sec::Float64 = Inf64,
    header::Vector{String}=[],
    data::Matrix{Union{String, Int, Float64}}=[],
    preprocessinstances::Bool=true,
    preprocesstimeout::Float64=Inf64,
    printintermediatetable::Bool = true,
    intermediatetableinterval::Float64 = 60.0,
    printsolutiongraphs::Bool=true,
    tablebackend::Symbol=:text
    )
    currenttime::Float64 = 0.0
    for (f, file) in enumerate(files)
        graph::SimpleWeightedGraph = read_file_as_simplegraph(file)
        n::Int = nv(graph) - 1
        data[rowot(Graphsize, f=f, ks=kl), colot(Graphsize, f=f, ks=kl)] = n
        for (k, size) in enumerate(sizes)
            treesize :: Int = 0
            if 0 <= size <= 1
                treesize = Int(floor(Int, (n * size)))
            else
                treesize = size
            end
            data[rowot(Treesize, f=f, k=k, ks=kl), colot(Treesize, f=f, k=k, ks=kl)] = treesize
            initialsolution::Union{Vector{Edge{Int}}, Nothing} = nothing
            lb::Int = 0
            ub::Int = typemax(Int)
            if preprocessinstances
                @debug "Preprocessing instance $file ($treesize nodes)."
                preprocesstime = @elapsed (initialsolution, lb, ub) = preprocess(graph, treesize, preprocesstimeout)
                @debug "Preprocessed instance $file ($treesize nodes) with lb: $(string(CYAN_FG(string(lb)))), ub: $(string(MAGENTA_FG(string(ub)))), rg: $(string(GREEN_FG(format_ratio_readable(1-(lb/ub))))), in $(format_seconds_readable(preprocesstime))."
            end
            opt::Union{Nothing,Int} = nothing
            for (s, solver) in enumerate(solvers)
                for (m, mode) in enumerate(modes)
                    line1::Union{String, Int, Float64} = ""
                    line2::Union{String, Int, Float64} = ""
                    line3::Union{String, Int, Float64} = ""
                    line4::Union{String, Int, Float64} = ""
                    @debug "Starting instance $file, with size $treesize, solver $solver and mode $mode."
                    kmstsolutionreport :: KMSTSolutionReport = kMST(graph, mode, treesize, solver, timeout_sec=timeout_sec, generate_timeout_sec=generate_timeout_sec, initialsolution=initialsolution, lowerbound=lb, upperbound=ub, printsolutiongraphs=printsolutiongraphs, debugmodels=debugmodels)
                    infostring = "Solved instance $file, with size $treesize, solver $solver and mode $mode."
                    if kmstsolutionreport.termination_status != MOI.OPTIMAL
                        infostring *=  " " * string(RED_FG("Did not find an optimal solution."))
                        if !isnothing(kmstsolutionreport.objective_value)
                            line1 = kmstsolutionreport.objective_value
                        end
                        if !isnothing(kmstsolutionreport.relative_gap)
                            line2 = format_ratio_readable(kmstsolutionreport.relative_gap)
                        end
                    else
                        line1 = format_seconds_readable(kmstsolutionreport.solve_time_sec)
                        @assert isnothing(opt) || opt == kmstsolutionreport.objective_value "Optimal value changed from $opt to $(kmstsolutionreport.objective_value) at $file with size $(treesize), solver $solver, mode $mode."
                        opt = kmstsolutionreport.objective_value
                    end
                    infostring *= " In $(format_seconds_readable(kmstsolutionreport.solve_time_sec)), with weight $(kmstsolutionreport.objective_value)"
                    if !isnothing(kmstsolutionreport.bnb_nodes)
                        line3 = kmstsolutionreport.bnb_nodes
                        infostring *= ", with $(kmstsolutionreport.bnb_nodes) BNB node$(kmstsolutionreport.bnb_nodes == 1 ? "" : "s")"
                    end
                    if !isnothing(kmstsolutionreport.violated_ineq)
                        line4 = kmstsolutionreport.violated_ineq
                        infostring *= ", with $(kmstsolutionreport.violated_ineq) violated inequalit$(kmstsolutionreport.bnb_nodes == 1 ? "y" : "ies")"
                    end
                    infostring *= "."
                    @info infostring
                    data[rowot(Sol, f=f, k=k, ks=kl, sol=s, sols=sl, m=m, ms=ml, line=1), colot(Sol, f=f, k=k, ks=kl, sol=s, sols=sl, m=m, ms=ml, line=1)] = line1
                    data[rowot(Sol, f=f, k=k, ks=kl, sol=s, sols=sl, m=m, ms=ml, line=2), colot(Sol, f=f, k=k, ks=kl, sol=s, sols=sl, m=m, ms=ml, line=2)] = line2
                    data[rowot(Sol, f=f, k=k, ks=kl, sol=s, sols=sl, m=m, ms=ml, line=3), colot(Sol, f=f, k=k, ks=kl, sol=s, sols=sl, m=m, ms=ml, line=3)] = line3
                    data[rowot(Sol, f=f, k=k, ks=kl, sol=s, sols=sl, m=m, ms=ml, line=4), colot(Sol, f=f, k=k, ks=kl, sol=s, sols=sl, m=m, ms=ml, line=4)] = line4
                    currenttime += kmstsolutionreport.solve_time_sec
                    if printintermediatetable && currenttime > intermediatetableinterval && (f != fl || k != kl || s != sl || m != ml)
                        pretty_table(data, header=header, backend=Val(tablebackend))
                        currenttime = 0.0
                    end
                    Base.GC.gc()
                end
            end
            o = "Unkn."
            if !isnothing(opt)
                o = opt
            end
            data[rowot(Opt, f=f, k=k, ks=kl, ms=ml, sols=sl), colot(Opt, f=f, k=k, ks=kl, ms=ml, sols=sl)] = o
        end
    end
    if tablebackend == :html
        pretty_table(data, header=header, backend=Val(tablebackend))
    else
        pretty_table(data, header=header, backend=Val(tablebackend), body_hlines=Int[4 * i for i in 1:(fl*kl)])
    end
end

struct KMSTSolutionReport
    termination_status :: MOI.TerminationStatusCode
    objective_value :: Union{Int, Nothing}
    relative_gap :: Union{Float64, Nothing}
    solve_time_sec :: Float64
    bnb_nodes :: Union{Int, Nothing}
    violated_ineq :: Union{Int, Nothing}
end

function kMST(graph :: SimpleWeightedGraph,
    mode::SolvingMode,
    k::Int,
    solver::Solver
    ;
    timeout_sec::Float64=Inf64,
    generate_timeout_sec::Float64=Inf64,
    initialsolution :: Union{Vector{Edge{Int}}, Nothing}=nothing,
    lowerbound::Int=0,
    upperbound::Int=typemax(Int),
    printsolutiongraphs::Bool=true,
    debugmodels::Bool=true) :: KMSTSolutionReport
    n::Int = nv(graph) - 1
    prelude_time::Float64 = 0
    prelude_time += @elapsed begin
        model = generate_model(graph, mode, k, solver, timeout_sec=timeout_sec, generate_timeout_sec=generate_timeout_sec, lowerbound=lowerbound, upperbound=upperbound,debugmodels=debugmodels)
    end
    if generate_timeout_sec != Inf64 && prelude_time > generate_timeout_sec
        return KMSTSolutionReport(MOI.TIME_LIMIT, upperbound, 1-(lowerbound/upperbound), prelude_time, 0, nothing)
    end
    nvars::Int = num_variables(model)
    if !isnothing(initialsolution) && solver != glpk
        prelude_time += @elapsed warmstart_model!(model, graph, mode, k, initialsolution)
    elseif isdebug()
        @debug "Model has $nvars variable$(nvars == 1 ? "" : "s")."
    end
    if isdebug()
        all_vars::Vector{VariableRef} = all_variables(model)
        bound::Int = count(v -> is_fixed(v) || (has_lower_bound(v) && has_upper_bound(v)), all_vars)
        if bound == nvars
            @debug "Model has no unbound variables."
        else
            free::Int = count(v -> !is_fixed(v) && !has_lower_bound(v) && !has_upper_bound(v), all_vars)
            onlyupperbound::Int = count(v -> !is_fixed(v) && !has_lower_bound(v) && has_upper_bound(v), all_vars)
            onlylowerbound::Int = count(v -> !is_fixed(v) && !has_upper_bound(v) && has_lower_bound(v), all_vars)
            @debug "Model has $bound ($(format_ratio_readable(bound,nvars))) bound variable$(bound == 1 ? "" : "s" ), $free ($(format_ratio_readable(free, nvars))) free variable$(free == 1 ? "" : "s"), $onlylowerbound ($(format_ratio_readable(onlylowerbound, nvars))) variable$(onlylowerbound == 1 ? "" : "s") with no upperbound, $onlyupperbound ($(format_ratio_readable(onlyupperbound, nvars))) variable$(onlyupperbound == 1 ? "" : "s") with no lowerbound."
        end

        constypes::Vector{Tuple{DataType, DataType}} = list_of_constraint_types(model)
        ncons::Int = 0
        debugstring::String = ""
        numcons::Dict{DataType, Int} = SortedDict()
        for (reftype, constype) in constypes
            c::Int = num_constraints(model, reftype, constype)
            if !(constype in keys(numcons))
                numcons[constype] = 0
            end
            numcons[constype] += c

            ncons += c
        end
        for (constype, c) in numcons
            if constype == MOI.LessThan{Float64}
                debugstring *= " $c ($(format_ratio_readable(c, ncons))) of type '≤'."
            elseif constype == MOI.GreaterThan{Float64}
                debugstring *= " $c ($(format_ratio_readable(c, ncons))) of type '≥'."
            elseif constype == MOI.EqualTo{Float64}
                debugstring *= " $c ($(format_ratio_readable(c, ncons))) of type '='."
            elseif constype == MOI.ZeroOne
                debugstring *= " $c ($(format_ratio_readable(c, ncons))) of type 'binary'."
            elseif constype == MOI.Integer
                debugstring *= " $c ($(format_ratio_readable(c, ncons))) of type 'integer'."
            elseif constype == MOI.Interval
                debugstring *= " $c ($(format_ratio_readable(c, ncons))) of type 'interval'."
            else
                debugstring *= " $c ($(format_ratio_readable(c, ncons))) of type '$constype'."
            end
        end
        @debug "Model has $ncons constraint$(ncons == 1 ? "" : "s")." * debugstring
    end

    @debug "Generated model in $(format_seconds_readable(prelude_time))."
    if generate_timeout_sec != Inf64 && prelude_time > generate_timeout_sec
        return KMSTSolutionReport(MOI.TIME_LIMIT, upperbound, 1-(lowerbound/upperbound), prelude_time, 0, nothing)
    end

    if mode == mcf && n >= 205
        return KMSTSolutionReport(MOI.MEMORY_LIMIT, upperbound, 1-(lowerbound/upperbound), prelude_time, 0, nothing)
    end

    kmstsolution = solve!(model, graph, k, mode)
    if printsolutiongraphs && !isnothing(kmstsolution.solution_graph)
        print_weighted_graph(kmstsolution.solution_graph)
    end

    objective_value :: Union{Int, Float64, Nothing} = kmstsolution.objective_value
    if objective_value == Inf64
        objective_value = nothing
    else
        objective_value = Int(round(Int, objective_value))
    end

    return KMSTSolutionReport(
        kmstsolution.termination_status,
        objective_value,
        kmstsolution.relative_gap,
        kmstsolution.solve_time_sec,
        kmstsolution.bnb_nodes,
        kmstsolution.violated_ineq)
end

function preprocess(graph::SimpleWeightedGraph, k::Int, timeout_sec::Float64=Inf64) :: Tuple{Vector{Edge{Int}}, Int, Int}
    nthreads :: Int = Base.Threads.nthreads()
    calc_time :: Float64 = 0

    lb::Int = 0
    ub::Int = typemax(Int)

    n::Int = nv(graph)

    lbfuture = Base.Threads.@spawn begin
    checked :: Set{Edge{Int}} = Set()
    tocheck_lock = Base.ReentrantLock()
    tocheck :: Vector{Edge{Int}} = []
    lowerbound_lock = Base.ReentrantLock()
    (lbtree, lb, s) = kruskal_heuristic(graph, k)
    append!(tocheck, setdiff(s, checked))
    for e in lbtree
        push!(checked, e)
    end

    while !isempty(tocheck)
        Threads.@threads for i = 1:max(1,div(nthreads-1,2))
            if !isempty(tocheck)
                unskippable = nothing
                lock(tocheck_lock)
                try
                    if !isempty(tocheck)
                        unskippable = pop!(tocheck)
                        while unskippable in checked && !isempty(tocheck)
                            unskippable = pop!(tocheck)
                        end
                        if unskippable in checked
                            unskippable = nothing
                        end
                    end
                finally
                    unlock(tocheck_lock)
                end
                if !isnothing(unskippable)
                    push!(checked, unskippable)
                    (lbtreecand, lbcand, s) = kruskal_heuristic(graph, k, unskippable=unskippable, bound=lb)
                    if !isnothing(lbtreecand)
                        for e in lbtreecand
                            push!(checked, e)
                        end
                    end
                    append!(tocheck, setdiff(s, checked))
                    if lbcand < lb
                        lock(lowerbound_lock)
                        try
                            if lbcand < lb
                                lbtree = lbtree
                                lb = lbcand
                            end
                        finally
                            unlock(lowerbound_lock)
                        end
                    end
                    if isempty(s)
                        lock(tocheck_lock)
                        try
                            empty!(tocheck)
                        finally
                            unlock(tocheck_lock)
                        end
                    end
                end
            end
        end
    end
    @assert !isnothing(lbtree) "No lower bound tree found."
    lb
    end


    ubtree::Vector{Edge{Int}} = []

    startnodes_lock = Base.ReentrantLock()
    upperbound_lock = Base.ReentrantLock()

    ubfuture = Base.Threads.@spawn begin
    startnodes::Vector{Int} = collect(Int,2:n)
    shuffle!(startnodes)
    checked :: Set{Int} = Set()

    while calc_time < timeout_sec && !isempty(startnodes)
        calc_time += @elapsed begin
            Threads.@threads for i = 1:max(1,div(nthreads-1,2)+((nthreads-1)%2))
                if !isempty(startnodes)
                    lock(startnodes_lock)
                    startnode = nothing
                    try
                        if !isempty(startnodes)
                            startnode = pop!(startnodes)
                            while startnode in checked && !isempty(startnodes)
                                startnode = pop!(startnodes)
                            end
                            if startnode in checked
                                startnode = nothing
                            end
                        end
                    finally
                        unlock(startnodes_lock)
                    end
                    if !isnothing(startnode)
                        heur::Union{Tuple{Vector{Edge{Int}}, Int}, Nothing} = prim_heuristic(graph, k, startnode=startnode, upperbound=ub)
                        if !isnothing(heur)
                            for e in heur[1]
                                for v in [src(e), dst(e)]
                                    push!(checked, v)
                                end
                            end
                            if heur[2] < ub
                                lock(upperbound_lock)
                                try
                                    if heur[2] < ub
                                        (ubtree, ub) = heur
                                    end
                                finally
                                    unlock(upperbound_lock)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    (ubtree,ub)
    end

    lb = fetch(lbfuture)
    (ubtree,ub) = fetch(ubfuture)

    @assert lb <= ub "Lowerbound has to be smaller than upperbound. !($lb < $ub)"

    return (ubtree, lb, ub)
end

function generate_model(graph::SimpleWeightedGraph,
        mode::SolvingMode,
        k::Int,
        solver::Solver;
        timeout_sec::Float64 = Inf64,
        generate_timeout_sec::Float64 = Inf64,
        lowerbound::Int = 0,
        upperbound::Int = typemax(Int),
        debugmodels::Bool=true)
    model = Model()
    if solver == glpk
        set_optimizer(model, GLPK.Optimizer)
    elseif solver == cplex
        set_optimizer(model, CPLEX.Optimizer)
    elseif solver == scip
        set_optimizer(model, SCIP.Optimizer)
    else
        @assert false "Defensive assert: Unhandled solver $solver."
    end
    if !isdebug()
        set_silent(model)
    else
        unset_silent(model)
        if solver == glpk
            set_optimizer_attribute(model, "msg_lev", GLPK.GLP_MSG_ON)
        end
    end
    if timeout_sec != Inf64
        if timeout < 0
            error("Timeout should be positive")
        end
        set_time_limit_sec(model, timeout)
    end
    if solver == cplex
        available_memory::Int = round(Int, Base.Sys.total_memory() // 2^20)
        work_memory::Int = max(1024, min(1024 * 8, div(available_memory, 3)))
        tree_limit::Int = max(1024 * 2, max(work_memory, div(9 * available_memory, 10)))
        @debug "There $(available_memory == 1 ? "is" : "are") $available_memory MiB available. Limiting work memory to $work_memory MiB, and setting tree limit to $tree_limit MiB."
        #set_optimizer_attribute(model, "CPX_PARAM_NODEFILEIND", 3)
        set_optimizer_attribute(model, "CPX_PARAM_TRELIM", tree_limit)
        set_optimizer_attribute(model, "CPX_PARAM_WORKMEM", work_memory)
    end
    generate_time::Float64 = 0
    @debug "Generating common variables."
    generate_time += @elapsed basic_kmst!(model, graph, k, lowerbound=lowerbound, upperbound=upperbound, isdebug=isdebug())
    if generate_time > generate_timeout_sec + 1.0
        return model
    end
    generate_time += @elapsed begin
    @debug "Generating $mode variables."
    if mode == mtz
        miller_tucker_zemlin(model, graph, k)
    elseif mode == scf
        single_commodity_flow!(model, graph, k)
    elseif mode == mcf
        multi_commodity_flow!(model, graph, k, generate_timeout_sec=generate_timeout_sec-generate_time, isdebug=isdebug())
    elseif mode == cec
        cycle_elimination_constraints!(model, graph, k)
    end
    end
    if debugmodels && generate_time < generate_timeout_sec + 1.0
        @debug model
    end
    return model
end

function warmstart_model!(model, graph::SimpleWeightedGraph, mode::SolvingMode, k::Int, solution::Vector{Edge{Int}})
    basic_kmst_warmstart!(model, graph, k, solution)
    if mode == mtz
        miller_tucker_zemlin_warmstart(model, graph, k, solution)
    end

    filled = check_kmst_warmstart!(model, graph, k)

    if mode == mtz
        filled += check_miller_tucker_zemlin_warmstart(model, graph, k)
    end

    @debug "Found $filled variable$(filled == 1 ? "" : "s") that could be filled."

    if isdebug()
        num_vars :: Int = JuMP.num_variables(model)
        vars :: Vector{VariableRef} = all_variables(model)
        start_values = start_value.(vars)
        with_startvalue :: Int = num_vars - count(isnothing, start_values)
        are_fixed::Int = count(is_fixed, vars)
        known::Int = count(v -> is_fixed(v) || !isnothing(start_value(v)), vars)
        @debug "Model has $num_vars variables, $known ($(format_ratio_readable(known, num_vars))) $(known == 1 ? "is" : "are") known, $with_startvalue ($(format_ratio_readable(with_startvalue, num_vars))) $(with_startvalue == 1 ? "has" : "have") a start value, $are_fixed ($(format_ratio_readable(are_fixed, num_vars))) $(are_fixed == 1 ? "is" : "are") fixed."
    end
end

function read_file_as_simplegraph(path::String)
    graph::Union{Nothing,SimpleWeightedGraph} = nothing
    vertexsize::Int = 0
    edgesize::Union{AbstractFloat,Int} = Inf
    edgeprog = nothing
    totaltime, totallines = open(path) do f
        linecounter::Int = 0
        timetaken = @elapsed for line in eachline(f)
            if linecounter <= 1
                if linecounter == 0
                    vertexsize = parse(Int, line)
                    graph = SimpleWeightedGraph(vertexsize)
                else
                    edgesize = parse(Int, line)
                    edgeProg = Progress(edgesize, dt=0.5,
                    barglyphs=BarGlyphs('[', '#', ['-','~','+','*','=','>','#'], ' ', ']'),
                    desc="Parsing edges:",
                    enabled=isdebug())
                end
            else
                elemsstrings = split(line, " ", limit=4)
                elems = [parse(Int, v) for v in elemsstrings]
                index = elems[1]
                node1 = elems[2]
                node2 = elems[3]
                weight = elems[4]
                if weight > 0
                    add_edge!(graph, node1 + 1, node2 + 1, weight)
                else
                    add_edge!(graph, node1 + 1, node2 + 1, min(1.0 / edgesize, 0.5 - (eps() * 2.0)))
                end
            end
            linecounter += 1
            if linecounter >= (edgesize + 2)
                break
            end
        end
        (timetaken, linecounter)
    end

    if totallines != 1
        linesplural = "s"
    else
        linesplural = ""
    end

    @debug ("Read $(string(YELLOW_FG(string(totallines)))) line$(linesplural) in $(string(BLUE_FG(format_seconds_readable(totaltime, 1)))) at $(string(CYAN_FG(format_seconds_readable(totaltime / totallines, 2)))) per line")
    es = ne(graph)
    if es != edgesize
        @warn "Input file claims $edgesize edge$(edgesize == 1 ? "" : "s") but got $es:"
        print_weighted_graph(graph, nothing, nv(graph), es)
    end
    return graph
end
