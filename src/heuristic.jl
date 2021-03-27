
using DataStructures
using LightGraphs, SimpleWeightedGraphs

function prim_heuristic(graph :: SimpleWeightedGraph, k :: Int, startnode :: Union{Nothing, Int} = nothing)
    es = edges(graph)
    es = collect(SimpleWeightedEdge, es)
    sort!(es, by=weight)
    k -= 1
    if isnothing(startnode)
        candidateedges = []
        candidatenodes = []
        for e in es
            i = src(e)
            j = dst(e)
            if i > 1 && j > 1
                if isempty(candidateedges)
                    push!(candidateedges, e)
                else
                    w = weight(e)
                    wh = weight(candidateedges[1])
                    if w < wh
                        candidateedges = [e]
                        candidatenodes = Dict(i => 0, j => 0)
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
                    candidatenodes[i] += Int(round(weight(e)))
                end
                if j in ks
                    candidatenodes[j] += Int(round(weight(e)))
                end
            end
        end
        (startnode, neighbouredgeweights) = findmin(candidatenodes)
    end
    @assert !isnothing(startnode) "No startnode"
    nvg = nv(graph)
    finished = Set([startnode])
    parents = zeros(Int, nvg)
    distmx = Int.(round.(weights(graph)))

    treesize = 0
    treeweight = 0

    while treesize < k
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

    kmst = [SimpleWeightedEdge(parents[v], v, distmx[parents[v], v]) for v in vertices(graph) if parents[v] != 0]

    return (kmst,treeweight)

end


