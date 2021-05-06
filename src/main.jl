
using ArgParse, Logging


function isprinted(level, logger=global_logger())
    return Logging.min_enabled_level(logger) <= level 
end

function isdebug(logger=global_logger())
    return isprinted(Logging.Debug, logger)
end

include("types.jl")

function parse_cmdargs()
    s = ArgParseSettings(prog="jkMST", version="1.0", add_version=true)

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
            default = ["1"]
            nargs = '+'
            arg_type = String
        "--solver", "-s"
            help = "The solver. Every instance is solved with every argument. (Allowed values: $(string([string(v) for v in instances(Solver)]))"
            default = ["cplex"]
            nargs = '+'
            arg_type = String
        "--timeout", "-t"
            help = "Maximum time in seconds spent on a single instance (3600s = 1h)."
            default = Inf64
            arg_type = Float64
        "--generate-timeout"
            help = "Maximum time in seconds spent on generating a model (600s = 10m, 3600s = 1h)."
            default = Inf64
            arg_type = Float64
        "--no-preprocess-instances"
            help = "Do not run Kruskal and Prim before solving file inorder to obtain better starting bounds."
            action = :store_true
        "--preprocess-timeout"
            help = "Maximum amount of time spent on preprocessing instances."
            default = Inf64
            arg_type = Float64
        "--no-print-solution-graphs", "--no-print-solutions"
            help = "Suppress the printing of the solution graphs after a solution is found."
            action = :store_true
        "--intermediate-table-interval"
            help = "Set the interval in seconds for printing the table in between instances."
            default = 60.0
            arg_type = Float64
        "--no-print-intermediate-table"
            help = "Do not print the table in between instances."
            action = :store_true
        "--no-print-models"
            help = "Do not print the models in debug mode."
            action = :store_true
        "--table-style"
            help = "Choose table style. (Allowed values: [\"text\", \"html\", \"latex\"])"
            default = "text"
            arg_type = String
        "--system-info"
            help = "Display what properties the system detects."
            action = :store_true
        "--verbose", "-v"
            help = "Increase verbosity."
            action = :count_invocations
        "--quiet", "-q"
            help = "Decrease verbosity."
            action = :count_invocations
    end

    return parse_args(s)
end

parsed_args = parse_cmdargs()

system_info = parsed_args["system-info"]


function test_system()
    rng = MersenneTwister()

    ops::Int = 10^7
    flops::Int = ops * 4
    operation_time::Float64 = 0.0
        for j in 1:ops
            s = rand(rng, Float64)
            m = rand(rng, Float64)
            operation_time += @elapsed r = m * s
            operation_time += @elapsed r = r - m
            operation_time += @elapsed r = s / m
            operation_time += @elapsed r = 5.0 + m
    end

    flopss = flops/operation_time
    prefix = ""
    if flopss > 1000
        prefix = "Kilo"
        flopss /= 1000
        if flopss > 1000
            prefix = "Mega"
            flopss /= 1000
            if flopss > 1000
                prefix = "Giga"
                flopss /= 1000
                if flopss > 1000
                    prefix = "Terra"
                    flopss /= 1000
                    if flopss > 1000
                        prefix = "Peta"
                        flopss /= 1000
                    end
                end
            end
        end
    end

    @info "CPUs: $(Base.Sys.CPU_THREADS) Threads: $(Base.Threads.nthreads()), estimated FLOPS: $(round(flopss; digits = 2)) $(prefix) FLOP/s, Total memory: $(round(Int,Base.Sys.total_memory() // 2^20)) MiB, Free memory: $(round(Int, Base.Sys.free_memory() // 2^20)) MiB"
    
end

if system_info
    using Random, Dates
    test_system()
    exit(0)
end

verbosity = parsed_args["verbose"]
quiet = parsed_args["quiet"]
level = loglevels[max(1, min(4, 2 + (quiet - verbosity)))]
logger = Logging.ConsoleLogger(stdout, level)
global_logger(logger)
@debug "LogLevel: $(Logging.min_enabled_level(logger)), isdebug: $(isdebug()), threads: $(Base.Threads.nthreads())."
files = parsed_args["file"]
directories = parsed_args["directory"]
modestrings = parsed_args["mode"]
kstrings = parsed_args["size"]
solverstrings = parsed_args["solver"]
timeout = parsed_args["timeout"]
generate_timeout = parsed_args["generate-timeout"]
preprocessinstances = !parsed_args["no-preprocess-instances"]
preprocesstimeout = parsed_args["preprocess-timeout"]
printsolutiongraphs = !parsed_args["no-print-solution-graphs"]
intermediatetableinterval = parsed_args["intermediate-table-interval"]
printintermediatetable = !parsed_args["no-print-intermediate-table"]
debugmodels = !parsed_args["no-print-models"]
tablestyle = parsed_args["table-style"]

if isempty(files) && isempty(directories)
    exit()
end

@debug "Compiling program."
includetime = @elapsed include("jkMST.jl")
@debug "Compiled program in $(format_seconds_readable(includetime))."

main(files=files,
    directories=directories,
    modestrings=modestrings,
    kstrings=kstrings,
    solverstrings=solverstrings,
    timeout_sec=timeout,
    generate_timeout_sec=generate_timeout,
    preprocessinstances=preprocessinstances,
    preprocesstimeout=preprocesstimeout,
    printsolutiongraphs=printsolutiongraphs,
    printintermediatetable=printintermediatetable,
    debugmodels=debugmodels,
    intermediatetableinterval=intermediatetableinterval,
    tablestyle=tablestyle
    )
