module jkMST

using Logging
using ArgParse
using Crayons.Box
using LightGraphs, SimpleWeightedGraphs
using JuMP
using GLPK

include("solvingMode.jl")
include("Util.jl")
include("Common.jl")
include("Mtz.jl")

function ArgParse.parse_item(::Type{SolvingMode}, x::AbstractString)
    return read(SolvingMode, x)
end

function parse_cmdargs()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--file", "-f"
            help = "A graph represented as file."
            required = true
        "--mode", "-m"
            help = "The solving method. (Allowed values: " * string([string(v) for v in instances(SolvingMode)]) * ")"
            default = "mtz"
        "--size", "-k"
            help = "Size of the k-MST."
            required = true
            arg_type = Int
        "--verbose", "-v"
            help = "Increase verbosity."
            action = :count_invocations
        "--quiet", "-q"
            help = "Decrease verbosity."
            action = :count_invocations
    end

    return parse_args(s)
end

const loglevels = [Logging.Debug, Logging.Info, Logging.Warn, Logging.Error]

function main()
    parsed_args = parse_cmdargs()
    verbosity :: Int = parsed_args["verbose"]
    quiet :: Int = parsed_args["quiet"]
    level = loglevels[max(1,min(4,2 + (quiet - verbosity)))]
    logger = Logging.ConsoleLogger(stdout, level)
    global_logger(logger)
    file :: String = parsed_args["file"]
    modestring :: String = parsed_args["mode"]
    mode :: SolvingMode = solvingmode_from_string(modestring)
    size :: Int = parsed_args["size"]
    kMST(file, mode, size)
end


function kMST(path :: String, mode :: SolvingMode, k:: Int)
    @debug ("Reading file: " * string(CYAN_FG(path)))
    @debug ("Solving in: $(string(MAGENTA_FG(string(mode)))) mode")
    @debug ("Spanning tree-size: $(string(GREEN_FG(string(k))))")
    graph = read_file_as_simplegraph(path)
    model = Model(GLPK.Optimizer)
    basic_kmst!(model, graph, k)
    if mode == mtz
        miller_tuckin_zemlin!(model, graph, k)
    end
    @debug model
    if mode == mtz || mode == scf || mode == mcf
        compact_formulation_solution!(model, graph, k)
    end
end


function read_file_as_simplegraph(path :: String)
    graph :: Union{Nothing, SimpleWeightedGraph} = nothing
    vertexsize :: Int = 0
    edgesize :: Int = 0
    totaltime, totallines = open(path) do f
        linecounter :: Int = 0
        timetaken = @elapsed for line in eachline(f)
            if linecounter <= 1
                if linecounter == 0
                    vertexsize = parse(Int, line)
                    graph = SimpleWeightedGraph(vertexsize)
                else
                    edgesize = parse(Int, line)
                end
            else
                elemsstrings = split(line, " ", limit = 4)
                elems = [parse(Int, v) for v in elemsstrings]
                index = elems[1]
                node1 = elems[2]
                node2 = elems[3]
                weight = elems[4]
                if weight > 0
                    add_edge!(graph, node1 + 1, node2 + 1, weight)
                else
                    add_edge!(graph, node1 + 1, node2 + 1, min(1.0/edgesize, 1.0 - (eps() * 2.0)))
                end
            end
            linecounter += 1
        end
        (timetaken, linecounter)
    end

    if totallines != 1
        linesplural = "s"
    else
        linesplural = ""
    end

    @debug ("Read $(string(YELLOW_FG(string(totallines)))) line$(linesplural) in $(string(BLUE_FG(format_seconds(totaltime, 0)))) at $(string(CYAN_FG(format_seconds(totaltime/totallines, 0)))) per line")
    es = ne(graph)
    if es != edgesize
        @warn "Input file claims $edgesize edges but got $es:" 
        print_weighted_graph(graph, nothing, nv(graph), es)
    end
    return graph
end


end
