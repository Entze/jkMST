
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
            required = true
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
verbosity = parsed_args["verbose"]
quiet = parsed_args["quiet"]
level = loglevels[max(1, min(4, 2 + (quiet - verbosity)))]
logger = Logging.ConsoleLogger(stdout, level)
global_logger(logger)
@debug "LogLevel: $(Logging.min_enabled_level(logger)), isdebug: $(isdebug())."
files = parsed_args["file"]
directories = parsed_args["directory"]
modestrings = parsed_args["mode"]
kstrings = parsed_args["size"]
solverstrings = parsed_args["solver"]
timeout = parsed_args["timeout"]
preprocessinstances = false
preprocesstimeout = Inf64
printsolutiongraphs = true
intermediatetableinterval = 60.0
printintermediatetable = true
tablestyle = "text"

include("jkMST.jl")

main(files=files,
    directories=directories,
    modestrings=modestrings,
    kstrings=kstrings,
    solverstrings=solverstrings,
    timeout_sec=timeout,
    preprocessinstances=preprocessinstances,
    preprocesstimeout=preprocesstimeout,
    printsolutiongraphs=printsolutiongraphs,
    printintermediatetable=printintermediatetable,
    intermediatetableinterval=intermediatetableinterval,
    tablestyle=tablestyle
    )
