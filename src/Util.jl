
const time_names_table = ["ns", "µs", "ms", "s", "m", "h", "d", "y"]
const time_step_table =     [1000, 1000, 1000, 60,  60,  24, 365]

function format_seconds(seconds, trail = -1)
    ns = seconds * 10.0^9.0
    times :: Vector{Int128} = []
    t = Int128(trunc(ns))
    for time_step :: Int128 in time_step_table
        (t::Int128,r::Int128) = divrem(t, time_step)
        append!(times, r)
        if t <= 0
            break
        end
    end
    if t > 0
        append!(times, t)
    end
    tlength = length(times)
    if trail < 0
        trail = tlength + 1
    end
    until = min(max(1, (tlength - trail)), tlength)
    times_with_unit = zip(times[until:end], time_names_table[until:tlength])
    times_with_units = [string(first(e)) * (e[2]) for e in times_with_unit]
    reverse!(times_with_units)
    return join(times_with_units, ", ")
end


function print_weighted_graph(graph :: SimpleWeightedGraph, level = nothing, numberofvertices = nv(graph), numberofedges = ne(graph), weights = weights(graph))
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
    line :: Int = 0
    for (i, e) in enumerate(edges(graph))
        msg :: String = "$line $(src(e) - 1) $(dst(e) - 1) $(Int(trunc(weight(e))))"
        if level === nothing
            println(msg)
        else
            @logmsg level msg
        end
        line += 1
    end
end
