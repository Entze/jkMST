
using DataStructures
using LightGraphs, SimpleWeightedGraphs

function prim_heuristic(graph :: SimpleWeightedGraph, k :: Int; startnode :: Union{Nothing, Int} = nothing, upperbound = typemax(Int))
    nvg = nv(graph)
    es = edges(graph)
    es = collect(SimpleWeightedEdge, es)
    sort!(es, by=weight)
    k -= 1
    if isnothing(startnode)
        candidateedges = []
        candidatenodes = Dict{Int, Int}()
        for e in es
            i = src(e)
            j = dst(e)
            if i > 1 && j > 1
                if isempty(candidateedges)
                    push!(candidateedges, e)
                    push!(candidatenodes, i => 0)
                    push!(candidatenodes, j => 0)
                else
                    w = weight(e)
                    wh = weight(candidateedges[1])
                    if w < wh
                        candidateedges = [e]
                        candidatenodes = Dict{Int,Int}(i => 0, j => 0)
                    elseif w == wh
                        push!(candidateedges, e)
                        push!(candidatenodes, i => 0)
                        push!(candidatenodes, j => 0)
                    end
                end
            end
        end
        @assert !isempty(candidatenodes) "No candidate nodes."
        for e in es
            i = src(e)
            j = dst(e)
            ks = keys(candidatenodes)
            if i in ks || j in ks
                if i in ks
                    candidatenodes[i] = get(candidatenodes, i, 0) + Int(round(weight(e)))
                end
                if j in ks
                    candidatenodes[j] = get(candidatenodes, j, 0) + Int(round(weight(e)))
                end
            end
        end
        (_, startnode) = findmin(candidatenodes)
    end
    @assert !isnothing(startnode) "No startnode"
    @assert startnode != 1 "Startnode is 1"
    @assert 1 <= startnode && startnode <= nvg "Startnode $startnode is out of bounds (1, $nvg)"
    time = @elapsed begin
        finished :: Set{Int} = Set([startnode])
        parents :: Vector{Int} = zeros(Int, nvg)

        treesize :: Int = 0
        treeweight :: Int = 0

        while treesize < k && treeweight < upperbound && !isempty(es)
            del = 0
            for e in es
                del += 1
                u = src(e)
                v = dst(e)
                ((u == 1 || v == 1) || (u in finished && v in finished) || !(u in finished || v in finished)) && continue
                w = Int(round(weight(e)))
                if u in finished
                    push!(finished, v)
                    parents[v] = u
                else
                    push!(finished, u)
                    parents[u] = v
                end
                treeweight += w
                treesize += 1
                break
            end
            deleteat!(es, del)
        end
        
    end

    if treeweight > upperbound
        @debug "Prim $startnode: >$upperbound after $treesize nodes in $(format_seconds_readable(time))."
        return nothing
    end

    kmst :: Vector{Edge{Int}} = [Edge{Int}(parents[v], v) for v in vertices(graph) if parents[v] != 0]

    @debug "Prim $startnode: $treeweight, in $(format_seconds_readable(time))."

    return (kmst,treeweight)

end


