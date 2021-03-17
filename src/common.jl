
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
    @debug "Starting to optimize."
    optimizationtime = @elapsed optimize!(model)
    @debug "Finished optimizing in $(format_seconds_readable(optimizationtime, 2))."
    ts = termination_status(model)
    if ts != MOI.OPTIMAL
        @warn "Optimizer did not find an optimal solution."
        return nothing
    end
    n :: Int = nv(graph)
    solution_graph :: SimpleGraph= get_solution_graph(model, n, k)
    components = connected_components(solution_graph)
    if ne(solution_graph) != k-1
        @error "Solution has wrong number of edges."
        @debug "Should have $(k - 1) edges but has $(ne(solution_graph))."
    end
    if n - k + 1 != length(components)
        @error "Solution is not connected." maxlog=1
        @debug "Should have $(n - k + 1) components but has $(length(components))."
    end
    for component in components
        cl = length(component)
        if cl != 1 && cl != k
            @error "Solution is not connected."
            @debug "Should only have components with size 1 or $k but found $cl: $(component)"
        end
    end
    ov = objective_value(model)
    @info "Found $k-MST with weight $(Int(trunc(ov)))."
    print_graph(solution_graph)
end

function get_solution_graph(model, n :: Int, k :: Int)
    solution_graph = SimpleGraph(n)
    variables = all_variables(model)
    x = variables[1:div(n * (n-1), 2)]

    for i in 2:n
        for j in (i+1):n
            index = (j - i) + sum(n-(a - 1) for a = 2:i)
            xv = Int(trunc(value(x[index])))
            if xv == 1
                add_edge!(solution_graph, i, j)
            end
        end
    end
    return solution_graph
end