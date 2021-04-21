
using DataStructures
using LightGraphs, SimpleWeightedGraphs

function prim_heuristic(graph::SimpleWeightedGraph,
    k::Int;
    startnode::Union{Nothing,Int}=nothing,
    upperbound=typemax(Int))::Union{Tuple{Vector{Edge{Int}}, Int}, Nothing}
    n = nv(graph)
    es = edges(graph)
    es = collect(SimpleWeightedEdge, es)
    sort!(es, by=weight)
    k -= 1
    if isnothing(startnode)
        candidateedges = []
        candidatenodes = Dict{Int,Int}()
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
    @assert 1 <= startnode && startnode <= n "Startnode $startnode is out of bounds (1, $n)"
    time = @elapsed begin
        finished::Vector{Bool} = zeros(Bool, n)
        finished[startnode] = true
        parents::Vector{Int} = zeros(Int, n)

        treesize::Int = 0
        treeweight::Int = 0
        lastnode::Int = -1

        while treesize < k && treeweight < upperbound && !isempty(es)
        del = 0
        for e in es
            del += 1
            u = src(e)
            v = dst(e)
            ((u == 1 || v == 1) || finished[u] == finished[v]) && continue
            w = Int(round(weight(e)))
            if finished[u]
                finished[v] = true
                parents[v] = u
                lastnode = u
            else
                finished[u] = true
                parents[u] = v
                lastnode = v
            end
            treeweight += w
            treesize += 1
            break
        end
        deleteat!(es, del)
    end
    end

    if treeweight > upperbound
        @debug "Prim $startnode: >$upperbound after $treesize nodes, in $(format_seconds_readable(time))."
        return nothing
    end

    kmst::Vector{Edge{Int}} = [Edge{Int}(parents[v], v) for v in vertices(graph) if parents[v] != 0]

    @debug "Prim $startnode: $treeweight, in $(format_seconds_readable(time))."

    return (kmst, treeweight)

end

function kruskal_heuristic(graph::SimpleWeightedGraph, k::Int; unskippable::Union{Nothing,Edge{Int}} = nothing, bound::Int = typemax(Int))::Tuple{Union{Vector{Edge{Int}}, Nothing}, Int, Set{Edge{Int}}}
    time = @elapsed begin
    skips::Int = 0
    skipped :: Set{Edge{Int}} = Set()
    n::Int= nv(graph)
    es::Vector{SimpleWeightedEdge} = collect(SimpleWeightedEdge, edges(graph))
    sort!(es, by=weight)

    intree::Vector{Bool} = zeros(Bool, n)
    parents::Vector{Int} = zeros(Int, n)

    s::Int = 0
    d::Int = 0
    w::Int = 0
    treeweight::Int = 0
    treesize::Int = 0

    if !isnothing(unskippable)
        s = src(unskippable)
        d = dst(unskippable)
        @assert 1 <= s <= n "s: $s is not in [1,n]."
        @assert 1 <= d <= n "d: $d is not in [1,n]."
        intree[s] = true
        intree[d] = true
        parents[d] = s
        treesize += 1
        for e in es
            (src(e) == s && dst(e) == d) || (src(e) == d && dst(e) == s) || continue
            treeweight += Int(round(Int, weight(e)))
        end
    end


    for e::SimpleWeightedEdge in Iterators.drop(es, n-1)
        treesize < (k-1) || break  

        s = src(e)
        d = dst(e)

        isnothing(unskippable) || !((src(e) == src(unskippable) && dst(e) == dst(unskippable)) || ((src(e) == dst(unskippable) && dst(e) == src(unskippable)))) || continue

        if intree[s] && intree[d]
            skips += 1
            push!(skipped, Edge{Int}(min(s,d), max(s,d)))
            continue
        end

        if intree[s] || !intree[d]
            parents[d] = s
        else
            parents[s] = d
        end

        treesize += 1

        intree[s] = true
        intree[d] = true
        treeweight += Int(round(Int,weight(e)))
        if treeweight > bound
            break
        end
    end
    end # time-measuring end

    if treeweight > bound
        @debug "Kruskal $(src(unskippable))-$(dst(unskippable)): >$bound after $treesize nodes, in $(format_seconds_readable(time)) with $(skips == 0 ? "no skips" : (skips == 1 ? "1 skip" : "$skips skips"))."
        return (nothing, treeweight, skipped)
    end

    @assert treesize == k-1 "Tree should have $(k-1) edges but has $treesize."

    @debug "Kruskal$(isnothing(unskippable) ? "" : " $(src(unskippable))-$(dst(unskippable))"): $treeweight, in $(format_seconds_readable(time)) with $(skips == 0 ? "no skips" : (skips == 1 ? "1 skip" : "$skips skips"))."

    return ([Edge{Int}(v, parents[v]) for v in vertices(graph) if parents[v] != 0], treeweight, skipped)

end