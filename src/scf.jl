
using JuMP
using LightGraphs, SimpleWeightedGraphs


function single_commodity_flow!(model, graph :: SimpleWeightedGraph, k :: Int)
    n::Int = nv(graph)
    es::Vector{SimpleWeightedEdge} = collect(edges(graph))
    sort!(es, by=dst)
    sort!(es, by=src, alg=MergeSort)

    @variables(model, begin
        f[i=1:n, j=1:n; has_edge(graph, i, j)], (lower_bound=0)
        z[i=2:n], (binary=true, lower_bound = 0, upper_bound = 1)
    end)

    alpha::Int = 1

    for i in 2:n
        if has_edge(graph, alpha, i)
            fix(variable_by_name(model, "y[$i,$alpha]"), 0, force=true) # Constraint (14)
            fix(f[i, alpha], 0, force=true) # Constraint (15)
            set_upper_bound(f[alpha,i], k)
        end
        for j in 2:n
            if has_edge(graph, i, j)
                set_upper_bound(f[i,j], k-1) # Constraint (16)
                set_upper_bound(f[j,i], k-1) # Constraint (16)
            end
        end
        @constraint(model, sum(variable_by_name(model, "y[$j,$i]") for j in 1:n if has_edge(graph, j, i)) == z[i]) # Constraint (17)
        #@constraint(model, sum(variable_by_name(model, "y[$i,$j]") for j in 2:n if has_edge(graph, j, i)) <= length(neighbors(graph, i)) - 1)
    end


    for i in 2:n
        @constraint(model,
        (sum(f[dst(e),i] for e in es if src(e) == i) + sum(f[src(e),i] for e in es if dst(e) == i)
        - (sum(f[i,dst(e)] for e in es if src(e) == i && dst(e) != alpha) + sum(f[i,src(e)] for e in es if dst(e) == i && src(e) != alpha)) == z[i])) # Constraint (18)
    end
    @constraint(model, source_flow, sum(f[alpha,j] for j = 2:n if has_edge(graph, alpha, j)) == k) # Constraint (19)
    @constraint(model, one_real_root, sum(variable_by_name(model, "y[$alpha,$j]") for j in 2:n if has_edge(graph, alpha, j)) == 1) # Constraint (20)

    for i in 1:n
        for j in 1:n
            if has_edge(graph, i, j)
                y_ij = variable_by_name(model, "y[$i,$j]")
                if i != 1 && j != 1
                    @constraint(model,
                        f[i,j] <= (k-1) * y_ij # Constraint (21)
                    )
                elseif i == 1
                    @constraint(model,
                        f[i,j] <= k * y_ij # Constraint (22)
                    )
                end
            end
        end
    end

    @constraint(model, sum(z[i] for i = 2:n) == k) # Constraint (23)

end