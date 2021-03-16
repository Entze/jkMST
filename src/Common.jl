
using Logging
using JuMP
using LightGraphs, SimpleWeightedGraphs


function basic_kmst!(model, graph :: SimpleWeightedGraph, k :: Int)
    n :: Int = nv(graph)
    @variables(model, begin
        x[i=1:n, (i+1):n], (binary=true)
        u[2:n], (integer=true, lower_bound=1, upper_bound=n-1)
    end)

    for i in 1:n
        for j in (i+1):n
            if !has_edge(graph, i, j)
                @constraint(model, x[i,j] == 0)
            end
        end
    end

    @constraint(model, sum(x[i, j] for i = 2:n for j=(i+1):n) == k-1)
    rawdistances = weights(graph)
    distances = trunc.(Int, rawdistances)
    @objective(model, Min, sum(sum(distances[i,j] * x[i,j] for j = (i+1):n) for i = 1:n))

end

function compact_formulation_solution!(model, graph :: SimpleWeightedGraph, k :: Int)
    optimize!(model)
    ts = termination_status(model)
    if ts != MOI.OPTIMAL
        @warn "Optimizer did not find an optimal solution."
        return nothing
    end
    n :: Int = nv(graph)
    solution_graph :: SimpleWeightedGraph= get_solution_graph(model, graph, n)
    print_solution(model, solution_graph, Logging.Info)
end

function print_solution(model, solution_graph :: SimpleWeightedGraph, level=Logging.Debug)
    ov = objective_value(model)
    @logmsg level "Optimal solution with weight $(ov) found."
    print_weighted_graph(solution_graph,level)
end

function get_solution_graph(model, graph :: SimpleWeightedGraph, n :: Int)
    
    solution_graph = SimpleWeightedGraph(n)
    variables = all_variables(model)
    x = variables[1:div(n * (n-1), 2)]

    distances = weights(graph)
    for i in 1:n
        for j in (i+1):n
            if value(x[(j - 1) + ((i-1) * (n - (i-1)))]) == 1
                add_edge!(solution_graph, i, j, distances[i,j])
            end
        end
    end
    return graph
end