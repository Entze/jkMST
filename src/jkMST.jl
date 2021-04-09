
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
    if name == Filename || name == Graphsize
        @assert !isnothing(f) "f is not defined."
        @assert !isnothing(ks) "ks is not defined."
        return 1 + (f - 1) * ks * 2
    elseif name == Treesize || name == Opt
        @assert !isnothing(f) "f is not defined."
        @assert !isnothing(ks) "ks is not defined."
        @assert !isnothing(k) "k is not defined."
        return 1 + (f - 1) * ks * 2 + (k - 1) * 2
    elseif name == Sol
        @assert !isnothing(f) "f is not defined."
        @assert !isnothing(ks) "ks is not defined."
        @assert !isnothing(k) "k is not defined."
        @assert !isnothing(line) "line is not defined."
        return 1 + (f - 1) * ks * 2 + (k - 1) * 2 + (line - 1)
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
                preprocessinstances::Bool=true,
                preprocesstimeout::Float64=Inf64,
                printsolutiongraphs::Bool=true,
                printintermediatetable::Bool=true,
                intermediatetableinterval::Float64=60.0,
                tablestyle::String="text"
                )
    fl::Int = length(files)
    dl::Int = length(directories)
    modes::Vector{SolvingMode} = map(solvingmode_from_string, modestrings)
    ml::Int = length(modes)
    ks::Vector{Union{Float64,Int}} = map(s -> begin
        v::Float64 = parse(Float64, s)
        if v > 1 && v == trunc(Float64, v)
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
            header=header,
            data=data,
            preprocessinstances=preprocessinstances,
            preprocesstimeout=preprocesstimeout,
            printintermediatetable=printintermediatetable,
            intermediatetableinterval=intermediatetableinterval,
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
    data::Matrix{Union{String, Int, Float64}} = fill("", fl * kl * 2, length(header))
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
                @debug "Preprocessing instance $file."
                preprocesstime = @elapsed (initialsolution, lb, ub) = preprocess(graph, treesize, preprocesstimeout)
                @debug "Preprocessed instance $file with lb: $(string(CYAN_FG(string(lb)))), ub: $(string(MAGENTA_FG(string(ub)))), in $(format_seconds_readable(preprocesstime))."
            end
            opt::Union{Nothing,Int} = nothing
            for (s, solver) in enumerate(solvers)
                for (m, mode) in enumerate(modes)
                    line1::Union{String, Int, Float64} = ""
                    line2::Union{String, Int, Float64} = ""
                    @debug "Starting instance $file, with size $treesize, solver $solver and mode $mode."
                    kmstsolutionreport = kMST(graph, mode, treesize, solver, timeout_sec=timeout_sec, initialsolution=initialsolution, lowerbound=lb, upperbound=ub,printsolutiongraphs=printsolutiongraphs)
                    infostring = "Solved instance $file, with size $size, solver $solver and mode $mode."
                    if kmstsolutionreport.termination_status != MOI.OPTIMAL
                        infostring *=  string(RED_FG("Did not find an optimal solution."))
                        line1 = kmstsolutionreport.objective_value
                        line2 = format_ratio_readable(kmstsolutionreport.relative_gap)
                    else
                        line1 = format_seconds_readable(kmstsolutionreport.solve_time_sec)
                        @assert isnothing(opt) || opt == kmstsolutionreport.objective_value "Optimal value changed from $opt to $(kmstsolutionreport.objective_value)."
                        opt = kmstsolutionreport.objective_value
                    end
                    infostring *= " In $(format_seconds_readable(kmstsolutionreport.solve_time_sec)), with weight $(kmstsolutionreport.objective_value)."
                    #TODO: write into lines

                end
            end
            if isnothing(opt)
                data[rowot(Opt, f=f, k=k, ks=kl, ms=ml, sols=sl), colot(Opt, f=f, k=k, ks=kl, ms=ml, sols=sl)] = "Unkn."
            end
        end
    end
    pretty_table(data, header, backend=tablebackend)
end

struct KMSTSolutionReport
    termination_status :: MOI.TerminationStatusCode
    objective_value :: Union{Int, Nothing}
    relative_gap :: Float64
    solve_time_sec :: Float64
end

function kMST(graph :: SimpleWeightedGraph, mode::SolvingMode, k::Int, solver::Solver; timeout_sec::Float64=Inf64, initialsolution :: Union{Vector{Edge{Int}}, Nothing}=nothing, lowerbound::Int=0, upperbound::Int=typemax(Int), printsolutiongraphs::Bool=true) :: KMSTSolutionReport
    n::Int = nv(graph) - 1
    model = generate_model(graph, mode, k, solver,timeout_sec=timeout_sec, lowerbound=lowerbound, upperbound=upperbound)
    if !isnothing(initialsolution) && solver != glpk
        warmstart_model!(model, graph, mode, k, initialsolution)
    end

    kmstsolution = solve!(model, graph, k)
    if printsolutiongraphs
        print_weighted_graph(kmstsolution.solution_graph)
    end

    return KMSTSolutionReport(kmstsolution.termination_status, Int(round(Int, kmstsolution.objective_value)), kmstsolution.relative_gap, kmstsolution.solve_time_sec)
end

function presolve(graph::SimpleWeightedGraph, k::Int)
    @debug "Presolving instance."
    totaltime = @elapsed begin
        n = nv(graph)
        candidates::Vector{Int} = collect(Int, 2:n)
        shuffle!(candidates)
        bestweight = typemax(Int)
        best = nothing
        searchtime = 0
        optimal::Bool = true
        for c in candidates
        searchtime += @elapsed begin
            res = prim_heuristic(graph, k, startnode=c, upperbound=bestweight)
            if !isnothing(res)
            (a, w) = res
            if w < bestweight
                best = a
                bestweight = w
            end
        end
            end
        
        if searchtime > 9.98
            optimal = false
            break;
        end
    end
    end
    debugstring = "Presolved instance in $(format_seconds_readable(totaltime)) with weight $(bestweight)."

    if optimal
        debugstring *= " Solution is $(string(GREEN_FG("optimal")))."
    else
        debugstring *= " Solution is $(string(RED_FG("suboptimal")))."
    end

    @debug debugstring
    sort!(best, by=dst)
    sort!(best, by=src, alg=MergeSort)
    return (best, bestweight)
end

function preprocess(graph::SimpleWeightedGraph, k::Int, timeout_sec::Float64=Inf64) :: Tuple{Vector{Edge{Int}}, Int, Int}
    return ([], 0, typemax(Int))
end

function generate_model(graph::SimpleWeightedGraph,
        mode::SolvingMode,
        k::Int,
        solver::Solver;
        timeout_sec::Float64 = Inf64,
        lowerbound::Int = 0,
        upperbound::Int = typemax(Int))
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

    basic_kmst!(model, graph, k, lowerbound=lowerbound, upperbound=upperbound)
    if mode == mtz
        miller_tuckin_zemlin!(model, graph, k)
    end
    @debug model
    return model
end

function warmstart_model!(model, graph::SimpleWeightedGraph, mode::SolvingMode, k::Int, solution::Vector{Edge{Int}})
    basic_kmst_warmstart!(model, graph, k, solution)
    if mode == mtz
        miller_tuckin_zemlin_warmstart!(model, graph, k, solution)
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
                    edgeProg = Progress(edgesize, dt=0.5, desc="Parsing edges:", enabled=isdebug())
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


