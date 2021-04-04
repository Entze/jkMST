module jkMST

using Logging
using ArgParse
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

function isprinted(level, logger = global_logger())
    return Logging.min_enabled_level(logger) <= level 
end

function isdebug(logger = global_logger())
    return isprinted(Logging.Debug, logger)
end

@enum Solver glpk scip cplex

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

@enum SolvingMode mtz scf mcf cec dcc

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


function parse_cmdargs()
    s = ArgParseSettings(prog="jkMST", version="1.0", add_version = true)

    @add_arg_table! s begin
        "--directory", "-d"
            help = "A directory containing graph files. All subdirectories are considered. Only files with extension .dat are read."
            nargs = '*'
            arg_type = String
        "--file", "-f"
            help = "A graph represented as file. Any file extension is permitted."
            nargs = '*'
            arg_type = String
        "--mode", "-m"
            help = "The solving method. Every instance is solved with every argument. (Allowed values: $(string([string(v) for v in instances(SolvingMode)])))"
            default = ["mtz"]
            nargs = '+'
            arg_type = String
        "--size", "-k"
            help = "Size of the k-MST. If k ∈ [0,1] interpreted as fraction of |V|. If k > 1 interpreted as absolute value. Every instance is solved with every argument."
            required = true
            nargs = '+'
            arg_type = String
        "--verbose", "-v"
            help = "Increase verbosity."
            action = :count_invocations
        "--solver", "-s"
            help = "The solver. Every instance is solved with every argument. (Allowed values: $(string([string(v) for v in instances(Solver)]))"
            default = ["cplex"]
            nargs = '+'
            arg_type = String
        "--quiet", "-q"
            help = "Decrease verbosity."
            action = :count_invocations
    end

    return parse_args(s)
end

const loglevels = [Logging.Debug, Logging.Info, Logging.Warn, Logging.Error]

function main()
    setuptime = @elapsed begin
        parsed_args = parse_cmdargs()
        verbosity::Int = parsed_args["verbose"]
        quiet::Int = parsed_args["quiet"]
        level = loglevels[max(1, min(4, 2 + (quiet - verbosity)))]
        logger = Logging.ConsoleLogger(stdout, level)
        global_logger(logger)
        @debug "LogLevel: $(Logging.min_enabled_level(logger)), isdebug: $(isdebug())."
        files::Vector{String} = parsed_args["file"]
        fl::Int = length(files)
        directories::Vector{String} = parsed_args["directory"]
        dl::Int = length(directories)
        modestrings::Vector{String} = parsed_args["mode"]
        modes::Vector{SolvingMode} = map(solvingmode_from_string, modestrings)
        ml::Int = length(modes)
        sizestrings::Vector{String} = parsed_args["size"]
        sizes :: Vector{Union{Float64, Int}} = map(s -> begin
            v :: Union{Float64, Int} = parse(Float64, s)
            if v > 1 && v == trunc(v)
                v = Int(trunc(v))
                return v
            end
            if !(0 <= v <= 1)
                error("Specified illegal k=$s. Should be Float ∈ [0,1] or Int if > 1.")
            end
            return v
        end, sizestrings)
        kl::Int = length(sizes)
        solverstrings::Vector{String} = parsed_args["solver"]
        solvers::Vector{Solver} = map(solver_from_string, solverstrings)
        sl::Int = length(solvers)
        @debug "Got $(fl) file$(fl == 1 ? "" : "s"), $(dl) $(dl == 1 ? "directory" : "directories"), $(ml) mode$(ml == 1 ? "" : "s"), $(kl) size$(kl == 1 ? "" : "s"), $(sl) solver$(sl == 1 ? "" : "s")."
        enumdirs = ProgressUnknown(desc="Enumerating files:", dt=0.5, enabled=isdebug())
        ProgressMeter.update!(enumdirs,fl)
        while !isempty(directories)
            directory = pop!(directories)
            for (root, dirs, fs) in walkdir(directory)
                for dir in dirs
                    push!(directories, dir)
                end
                for file in fs
                    if endswith(file, ".dat")
                        push!(files, joinpath(root, file))
                        ProgressMeter.next!(enumdirs)
                    end
                end
            end
        end
        ProgressMeter.finish!(enumdirs)
        fl = length(files)
        if fl + dl <= 0
            error("At least one file or directory required.")
        end
        headers = ["Graph", "|V|", "k", "OPT"]
        if sl > 1
            for s in solverstrings
                push!(headers, s)
                append!(headers, modestrings)
            end
        else
            append!(headers, modestrings)
        end
        data :: Matrix{Union{String, Int}} = fill("", fl * kl, length(headers))
        for (i, file) in enumerate(files)
            data[1 + ((i-1) * kl), 1] = file
        end
    end
    alltime = @elapsed begin
        for (f,file) in enumerate(files)
            n :: Union{Nothing, Int} = nothing
            for (k,size) in enumerate(sizes)
                opt :: Union{Nothing, Int} = nothing
                for (s,solver) in enumerate(solvers)
                    for (m, mode) in enumerate(modes)
                        (solvingtime, objectivevalue, n_, k_) = kMST(file, mode, size, solver)
                        if isnothing(n)
                            n = n_
                        end
                        @assert n == n_ "Size of graph changed, should be $n is $n_. Graph: $file Solver: $solver Mode: $mode Size: $size."
                        if isnothing(opt)
                            opt = objectivevalue
                        end
                        @assert opt == objectivevalue "Objective value of graph changed, should be $opt is $objectivevalue. Graph: $file Solver: $solver Mode: $mode Size: $size."
                        data[((f-1) * kl) + k, 3] = k_
                        column = 4 + Int(sl != 1) + m + ((s-1) * (Int(sl != 1) + ml))
                        data[((f-1) * kl) + k, column] = format_seconds_readable(solvingtime)
                    end
                end
                data[((f-1) * kl) + k, 4] = opt
            end
            data[1 + ((f-1) * kl), 2] = n
        end
    end
    @debug "Finished in $(format_seconds_readable(alltime)) + $(format_seconds_readable(setuptime)) setup time."
    pretty_table(data, headers)
end


function kMST(path::String, mode::SolvingMode, k::Union{Int, Float64}, solver::Solver) :: Tuple{Float64, Int, Int, Int}
    @debug ("Solving with: $(string(YELLOW_FG(uppercase(string(solver))))).")
    @debug ("Reading file: $(string(CYAN_FG(path))).")
    @debug ("Solving in: $(string(MAGENTA_FG(string(mode)))) mode.")
    @debug ("Spanning tree-size: $(string(GREEN_FG(string(k)))).")
    graph = read_file_as_simplegraph(path)
    n::Int = nv(graph) - 1
    k_ = k
    if 0 <= k <= 1
        k_ = Int(floor(k * n))
    end
    inittime = @elapsed begin
        if solver != glpk
            (pre, ub) = presolve(graph, k_)
            model = generate_model(graph, mode, k_, solver, upperbound=ub)
            warmstart_model!(model, graph, mode, k_, pre)
        else
            model = generate_model(graph, mode, k_, solver)
        end
    end
    ov :: Int = -1
    if mode == mtz || mode == scf || mode == mcf
        (ov,solution_graph) = compact_formulation_solution!(model, graph, k_)
    end
    st = solve_time(model)
    @info "Found $k-MST of weight $(ov) for $path in $mode formulation with $solver in $(format_seconds_readable(st)) + $(format_seconds_readable(inittime)) presolving and initialization."
    print_weighted_graph(solution_graph, nothing)
    return (st, ov, n, k_)
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

function generate_model(graph::SimpleWeightedGraph,
        mode::SolvingMode,
        k::Int,
        solver::Solver;
        lowerbound=typemin(Int),
        upperbound=typemax(Int))
    @debug "Generating model."
    model = Model()
    modeltime = @elapsed begin
        if solver == glpk
            optimizer = GLPK.Optimizer
        elseif solver == cplex
            optimizer = CPLEX.Optimizer
        else
            optimizer = SCIP.Optimizer
        end
        set_optimizer(model, optimizer)
        if !isdebug()
            set_silent(model)
        end

        basic_kmst!(model, graph, k, lowerbound=lowerbound, upperbound=upperbound)
        if mode == mtz
            miller_tuckin_zemlin!(model, graph, k)
        end
    end
    @debug "Generated model in $(format_seconds_readable(modeltime)):"
    @debug model
    return model
end

function warmstart_model!(model, graph::SimpleWeightedGraph, mode::SolvingMode, k::Int, solution)
    @debug "Warmstart:"
    totaltime = @elapsed begin
        basic_kmst_warmstart!(model, graph, k, solution)
        if mode == mtz
            miller_tuckin_zemlin_warmstart!(model, graph, k, solution)
        end
    end
    @debug "Warmstart finished in $(format_seconds_readable(totaltime))."
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


end
