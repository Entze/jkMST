module jkMST

using Logging
using ArgParse
using Crayons.Box
using LightGraphs, SimpleWeightedGraphs
using JuMP
using GLPK, SCIP, CPLEX
using Random

include("common.jl")
include("heuristic.jl")
include("util.jl")
include("mtz.jl")

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
    return SolvingMode(i)
end


function parse_cmdargs()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--file", "-f"
            help = "A graph represented as file."
            required = true
        "--mode", "-m"
            help = "The solving method. (Allowed values: $(string([string(v) for v in instances(SolvingMode)])))"
            default = "mtz"
        "--size", "-k"
            help = "Size of the k-MST."
            required = true
            arg_type = Int
        "--verbose", "-v"
            help = "Increase verbosity."
            action = :count_invocations
        "--solver", "-s"
            help = "Choose solver. (Allowed values: $(string([string(v) for v in instances(Solver)]))"
            default = "scip"
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
        file::String = parsed_args["file"]
        modestring::String = parsed_args["mode"]
        mode::SolvingMode = solvingmode_from_string(modestring)
        size::Int = parsed_args["size"]
        solverstring::String = parsed_args["solver"]
        solver::Solver = solver_from_string(solverstring)
    end
    alltime = @elapsed kMST(file, mode, size, solver)
    @debug "Finished in $(format_seconds_readable(alltime)) + $(format_seconds_readable(setuptime)) setup time."
end


function kMST(path::String, mode::SolvingMode, k::Int, solver::Solver)
    @debug ("Solving with: $(string(YELLOW_FG(uppercase(string(solver))))).")
    @debug ("Reading file: $(string(CYAN_FG(path))).")
    @debug ("Solving in: $(string(MAGENTA_FG(string(mode)))) mode.")
    @debug ("Spanning tree-size: $(string(GREEN_FG(string(k)))).")
    graph = read_file_as_simplegraph(path)
    (pre, ub) = presolve(graph, k)
    model = generate_model(graph, mode, k, solver, upperbound=ub)
    warmstart_model!(model, graph, mode, k, pre)
    if mode == mtz || mode == scf || mode == mcf
        compact_formulation_solution!(model, graph, k)
    end
end

function presolve(graph::SimpleWeightedGraph, k::Int)
    @debug "Presolving instance."
    totaltime = @elapsed begin
        n = nv(graph)
        candidates :: Vector{Int} = collect(Int, 2:n)
        shuffle!(candidates)
        bestweight = typemax(Int)
        best = nothing
        searchtime = 0
        optimal::Bool = true
        for c in candidates
            searchtime += @elapsed begin
            res = prim_heuristic(graph, k, startnode=c, upperbound=bestweight)
            if !isnothing(res)
                (a,w) = res
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
    if mode == mtz
        mtz_start_values_from_heuristic!(model, graph, k, solution)
    end
end

function read_file_as_simplegraph(path::String)
    graph::Union{Nothing,SimpleWeightedGraph} = nothing
    vertexsize::Int = 0
    edgesize::Union{AbstractFloat,Int} = Inf
    totaltime, totallines = open(path) do f
        linecounter::Int = 0
        timetaken = @elapsed for line in eachline(f)
            if linecounter <= 1
                if linecounter == 0
                    vertexsize = parse(Int, line)
                    graph = SimpleWeightedGraph(vertexsize)
                else
                    edgesize = parse(Int, line)
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
        @warn "Input file claims $edgesize edges but got $es:" 
        print_weighted_graph(graph, nothing, nv(graph), es)
    end
    return graph
end


end
