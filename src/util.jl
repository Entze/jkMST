const time_names_table = ["ns", "µs", "ms", "s", "m", "h", "d", "y"]
const time_step_table =     [1000, 1000, 1000, 60,  60,  24, 365]

using Formatting
using Logging
using LightGraphs, SimpleWeightedGraphs

function format_seconds_readable(seconds, digits=2)::String
    fspec = FormatSpec(".$(digits)f")
    if 0.1 <= seconds <= 120.0 || seconds == 0.0
        return "$(fmt(fspec, seconds))s"
    end

    if seconds > 60.0
        # Up
        minutes = seconds / 60.0
        if minutes > 120.0
            hours = minutes / 60.0
            return "$(fmt(fspec, hours))h"
        end
        return "$(fmt(fspec, minutes))m"
    else
        # Down
        milliseconds = seconds * 1000.0
        if milliseconds < 0.1
            microseconds = milliseconds * 1000.0
            if microseconds < 0.1
                nanoseconds = milliseconds * 1000.0
                return "$(fmt(fspec, nanoseconds))ns"
            end
            return "$(fmt(fspec, microseconds))µs"
        end
        return "$(fmt(fspec, milliseconds))ms"
    end
    return string(seconds)
end

function format_ratio_readable(amount::Union{Int, Float64}, max::Union{Int, Float64} = 1, digits=2)::String
    if amount == 0
        if max == 0
            return "100%"
        else
            return "0%"
        end
    end
    fspec = FormatSpec(".$(digits)f")
    perc::Union{Float64, Rational{Int}} = 0.0
    if typeof(amount) == Float64
        perc = amount/max
    else
        max::Int
        perc = amount//max
    end
    perc *= 100
    if 1//5 <= perc || perc == 0
        return "$(fmt(fspec, perc))%"
    end

    perc *= 100

    if 1//5 <= perc
        return "$(fmt(fspec, perc))‰"
    end

    perc *= 100

    if 1//5 <= perc
        return "$(fmt(fspec, perc))‱"
    end

    perc *= 100

    if 1//5 <= perc
        return "$(fmt(fspec, perc))pcm"
    end

    perc *= 100

    return "$(fmt(fspec, perc))ppm"

end

function print_weighted_graph(graph::SimpleWeightedGraph, level=nothing, numberofvertices=nv(graph), numberofedges=ne(graph))
    if level === nothing
        println(numberofvertices)
    else
        @logmsg level numberofvertices
    end
    if level === nothing
        println(numberofedges)
    else
        @logmsg level numberofedges
    end
    line::Int = 0
    for (i, e) in enumerate(edges(graph))
        msg::String = "$line $(src(e) - 1) $(dst(e) - 1) $(Int(round(weight(e))))"
        if level === nothing
            println(msg)
        else
            @logmsg level msg
        end
        line += 1
    end
end

function print_graph(graph::AbstractGraph, level=nothing, weight=nothing, numberofvertices=nv(graph), numberofedges=ne(graph))
    if level === nothing
        println(numberofvertices)
    else
        @logmsg level numberofvertices
    end
    if level === nothing
        println(numberofedges)
    else
        @logmsg level numberofedges
    end
    if weight !== nothing
        if level === nothing
            println(weight)
        else
            @logmsg level weight
        end
    end
    line::Int = 0
    for (i, e) in enumerate(edges(graph))
        msg::String = "$line $(src(e) - 1) $(dst(e) - 1)"
        if level === nothing
            println(msg)
        else
            @logmsg level msg
        end
        line += 1
    end
end
