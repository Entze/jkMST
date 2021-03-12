module jkMST

include("solvingMode.jl")

using ArgParse
using Crayons.Box
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

    end

    return parse_args(s)
end

function main()
    parsed_args = parse_cmdargs()
    file = parsed_args["file"]
    modestring = parsed_args["mode"]
    mode = solvingmode_from_string(modestring)
    size = parsed_args["size"]
    println("Reading file: ", CYAN_FG(file))
    println("Solving in: ", MAGENTA_FG(string(mode)), " mode")
    println("Spanning tree-size: ",  GREEN_FG(string(size)))
end

end
