
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
        d[i=2:n], (integer=true, lower_bound = 0, upper_bound = k)
        z[i=2:n], (binary=true, lower_bound = 0, upper_bound = 1)
    end)

    alpha::Int = 1

    for i in (alpha+1):n
        @constraint(model, 
        (sum(f[dst(e),i] for e in es if src(e) == i) + sum(f[src(e),i] for e in es if dst(e) == i))
        - (sum(f[i,dst(e)] for e in es if src(e) == i) + sum(f[i,src(e)] for e in es if dst(e) == i)) == 1)
    end
    @constraint(model, source_flow, sum(f[alpha,j] for j = 1:n if has_edge(graph, alpha, j)) - sum(f[i,alpha] for i = 1:n if has_edge(graph, i, alpha)) == n - 1)

    for i in 1:n
        for j in 1:n
            if has_edge(graph, i, j)
                y_ij = variable_by_name(model, "y[$i,$j]")
                @constraint(model,
                    f[i,j] <= (k-1) * y_ij
                )
            end
        end
    end

    for j in (alpha+1):n
        @constraint(model,
            sum(variable_by_name(model, "y[$i,$j]") for i in 2:n if has_edge(graph, i, j)) 
            + sum(variable_by_name(model, "y[$j,$k_]") for k_ in 2:n if has_edge(graph, j, k_)) == d[j])
        @constraints(model, begin
            z[j] * k >= d[j]
            z[j] <= d[j]
        end)
    end

    @constraint(model, sum(z[i] for i = 2:n) == k)

end