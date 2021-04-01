
using JuMP
using LightGraphs, SimpleWeightedGraphs

function miller_tuckin_zemlin!(model, graph :: SimpleWeightedGraph, k :: Int)
    n :: Int = nv(graph)
    m = max(k, n-k) + 1
    @variables(model, begin
        u[1:n], (integer=true, upper_bound = k+1)
        d[2:n], (integer=true, lower_bound = 0, upper_bound = m)
        y[2:n], (binary=true, lower_bound = 0, upper_bound = 1)
    end)
    omega :: Int = 1
    @constraint(model, u[omega] == 0)
    for i in (omega+1):n
        @constraints(model, begin
            u[i] >= 1
            d[i] <= k
        end)
        for j in 1:n
            if has_edge(graph, i, j)
                x_ij = variable_by_name(model, "x[$i,$j]")
                @constraint(model,
                    u[i] >= u[j] + x_ij - m * (1 - x_ij)
                )
            end
        end
    end

    for i in (omega+1):n
        @constraint(model, sum(variable_by_name(model, "x[$i,$j]") for j = 1:n if has_edge(graph, i, j)) == 1)
    end
    @constraint(model, sum(variable_by_name(model, "x[$omega,$j]") for j = (omega + 1):n) == 0)

    for j in (omega+1):n
        @constraint(model,
            sum(variable_by_name(model, "x[$i,$j]") for i in 2:n if has_edge(graph, i, j)) 
            + sum(variable_by_name(model, "x[$j,$k_]") for k_ in 2:n if has_edge(graph, j, k_)) == d[j])
        @constraints(model, begin
            y[j] * k >= d[j]
            y[j] <= d[j]
        end)
    end

    @constraint(model, sum(y[i] for i = 2:n) == k)

end

function mtz_start_values_from_heuristic!(model, graph :: SimpleWeightedGraph, k :: Int, solution)
    nvg = nv(graph)
    y = zeros(Int, nvg)
    connected :: Vector{Vector{Int}} = [[] for i=1:nvg]
    for e in solution
        i = src(e)
        j = dst(e)
        y[i] = 1
        y[j] = 1
        push!(connected[i], j)
        push!(connected[j], i)
    end
    s = nothing
    for c in 2:nvg
        if !isempty(connected[c])
            s = c
            break
        end
    end

    stack = [s]
    marked :: Vector{Bool} = zeros(Bool, nvg)
    dagsort :: Vector{Int} = []
    while !isempty(stack)
        node = pop!(stack)
        neighbours = connected[node]
        for n in neighbours
            if !marked[n]
                push!(stack, n)
                set_start_value(variable_by_name(model, "x[$n,$node]"), 1)
                set_start_value(variable_by_name(model, "x[$node,$n]"), 0)
            end
        end
        marked[node] = true
        pushfirst!(dagsort, node)
    end

    for i in 2:nvg
        set_start_value(variable_by_name(model, "y[$i]"), y[i])
        if y[i] == 0
            set_start_value(variable_by_name(model, "x[$i,1]"), 1)
            for j in 2:nvg
                if has_edge(graph, i, j)
                    set_start_value(variable_by_name(model, "x[$i,$j]"), 0)
                    set_start_value(variable_by_name(model, "x[$j,$i]"), 0)
                end
            end
        end
    end
end
