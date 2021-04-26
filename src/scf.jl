
using JuMP
using LightGraphs, SimpleWeightedGraphs


function single_commodity_flow!(model, graph :: SimpleWeightedGraph, k :: Int)
    n::Int = nv(graph)
    m::Int = ne(graph)
    deg::Int = max(k, n-k)
    es::Vector{SimpleWeightedEdge} = collect(edges(graph))
    sort!(es, by=dst)
    sort!(es, by=src, alg=MergeSort)

    @variables(model, begin
        f[i=1:n, j=1:n; has_edge(graph, i, j)], (lower_bound=0, upper_bound=k)
        z[i=2:n], (binary=true, lower_bound = 0, upper_bound = 1)
    end)

    alpha::Int = 1

    for i in (alpha+1):n
        @constraint(model, 
        (sum(f[dst(e),i] for e in es if src(e) == i) + sum(f[src(e),i] for e in es if dst(e) == i))
        - (sum(f[i,dst(e)] for e in es if src(e) == i) + sum(f[i,src(e)] for e in es if dst(e) == i)) == z[i])

    end
    @constraint(model, source_flow, sum(f[alpha,j] for j = 2:n if has_edge(graph, alpha, j)) - sum(f[i,alpha] for i = 2:n if has_edge(graph, i, alpha)) == k)
    @constraint(model, one_real_root, sum(variable_by_name(model, "y[$alpha,$j]") for j in 2:n if has_edge(graph, alpha, j)) == 1)

    for i in 1:n
        for j in 1:n
            if has_edge(graph, i, j)
                y_ij = variable_by_name(model, "y[$i,$j]")
                if i != 1 && j != 1
                    @constraint(model,
                        f[i,j] <= (k-1) * y_ij
                    )
                elseif i == 1
                    @constraint(model,
                        f[i,j] <= k * y_ij
                    )
                end
            end
        end
    end

    for i in 2:n
        @constraint(model, sum(variable_by_name(model, "y[$j,$i]") for j in 2:n if has_edge(graph, j, i)) <= (length(neighbors(graph, i)) - 1) * z[i]) # z is one if at least one edge goes out
        # @constraint(model, sum(variable_by_name(model, "y[$j,$i]") for j in 2:n if has_edge(graph, j, i)) >= z[i]) # z is zero if no edges go in
    end

    @constraint(model, sum(z[i] for i = 2:n) == k)

end