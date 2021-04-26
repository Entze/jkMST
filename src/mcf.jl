
using JuMP
using LightGraphs, SimpleWeightedGraphs


function multi_commodity_flow!(model, graph :: SimpleWeightedGraph, k :: Int)
    n::Int = nv(graph)
    es::Vector{SimpleWeightedEdge} = collect(edges(graph))
    sort!(es, by=dst)
    sort!(es, by=src, alg=MergeSort)

    @variables(model, begin
        f[m=2:n, i=1:n, j=1:n; has_edge(graph, i, j)], (lower_bound=0, upper_bound=1)
        z[2:n], (binary = true, lower_bound=0, upper_bound=1)
    end)

    alpha::Int = 1

    for i in 2:n
        if has_edge(graph, i, alpha)
            for m in 2:n
                fix(f[m, i, alpha], 0, force=true)
                fix(variable_by_name(model, "y[$i,$alpha]"), 0, force=true)
            end
        end
    end

    #@constraint(model, source_flow, sum(sum(f[m, alpha, j] for j in 2:n if has_edge(graph, alpha, j)) for m in 2:n) == k)
    @constraint(model, one_real_root, sum(variable_by_name(model, "y[$alpha,$j]") for j in 2:n if has_edge(graph, alpha, j)) == 1)
    for i in 1:n
        for j in 1:n
            if has_edge(graph, i, j)
                for m in 2:n
                    @constraint(model, f[m, i, j] <= variable_by_name(model, "y[$i,$j]"))
                end
            end
        end
    end
    for m in 2:n
        @constraint(model, sum(f[m, i, m] for i in 1:n if has_edge(graph, i, m))
                            - sum(f[m, m, j] for j in 2:n if has_edge(graph, m, j)) == z[m])
        @constraint(model, sum(f[m, alpha, j] for j in 2:n if has_edge(graph, alpha, j)) == z[m])
        for v in 2:n
            if m != v
                @constraint(model, sum(f[m, i, v] for i in 1:n if has_edge(graph, i, v))
                                    - sum(f[m, v, j] for j in 2:n if has_edge(graph, v, j)) == 0)
            end
        end
    end

    @constraint(model, treesize, sum(z[i] for i in 2:n) == k)

end
