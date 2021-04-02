
using Logging
using JuMP
using LightGraphs, SimpleWeightedGraphs
using Crayons.Box

function value_rounded(v)
    return Int(round(value(v)))
end

function basic_kmst!(model,
        graph :: SimpleWeightedGraph,
        k :: Int;
        lowerbound :: Int = typemin(Int),
        upperbound :: Int = typemax(Int))
    n :: Int = nv(graph)
    es = collect(SimpleWeightedEdge, edges(graph))
    sort!(es, by=weight)

    lb = max(sum(Int(round(weight(e))) for e in Iterators.take(Iterators.drop(es, n-1), k-1)), lowerbound)
    ub = min(sum(Int(round(weight(e))) for e in Iterators.drop(es, n-(k-1))), upperbound)
    @debug "Objective has $(string(CYAN_FG(string(lb)))) as natural lowerbound."
    @debug "Objective has $(string(MAGENTA_FG(string(ub)))) as natural upperbound."
    @variables(model,begin
        x[i=1:n, j=1:n; has_edge(graph, i, j)], (binary=true, lower_bound=0, upper_bound=1) # Set bounds, that solver doesn't have to set it implicitly
        o, (integer=true, lower_bound=lb, upper_bound = ub)
    end)

    rawdistances = weights(graph)
    distances = round.(Int, rawdistances)
    @constraints(model, begin
        sum(x[i,j] for i = 2:n for j=2:n if has_edge(graph, i ,j)) == k-1 # there are k-1 edges in a tree of size k
        sum(sum(distances[i,j] * x[i,j] for j = 1:n if has_edge(graph,i,j)) for i = 1:n) == o
    end)
    @objective(model, Min, o)

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
    solution_graph :: SimpleWeightedGraph= get_solution_graph(model, graph)
    components = connected_components(solution_graph)
    if ne(solution_graph) != k-1
        @error "Solution has wrong number of edges."
        @debug "Should have $(k - 1) edges but has $(ne(solution_graph))."
    end
    if n - k + 1 != length(components)
        @error "Solution has wrong number of components."
        @debug "Should have $(n - k + 1) components but has $(length(components))."
    end
    for component in components
        cl = length(component)
        if cl != 1 && cl != k
            @error "Solution is not connected." maxlog=1
            @debug "Should only have components with size 1 or $k but found with size $cl: $(component)."
        end
    end
    ov = objective_value(model)
    ov_ = Int(trunc(ov))
    @info "Found $k-MST with weight $(ov_)."
    print_weighted_graph(solution_graph, nothing)
    w = sum(weight(e) for e in edges(solution_graph))
    w_ = Int(trunc(w))
    if w_ != ov_
        @warn "MILP solver claims objective value is $(ov_), however got $(w_)."
    end
end

function get_solution_graph(model, graph :: SimpleWeightedGraph)
    n :: Int = nv(graph)
    solution_graph = SimpleWeightedGraph(n)
    es = edges(graph)
    for e in es
        i = src(e)
        j = dst(e)
        w = weight(e)
        if i > 1 && j > 1
            x_ij = variable_by_name(model, "x[$i,$j]")
            xv = Int(round(value(x_ij)))
            if xv == 1
                add_edge!(solution_graph, i, j, w)
            else
                x_ji = variable_by_name(model, "x[$j,$i]")
                xv = Int(round(value(x_ji)))
                if xv == 1
                    add_edge!(solution_graph, i, j, w)
                end
            end
        end
    end
    
    return solution_graph
end