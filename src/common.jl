
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
    m :: Int = ne(graph)
    es = collect(SimpleWeightedEdge, edges(graph))
    sort!(es, by=weight)

    lb = max(sum(Int(round(weight(e))) for e in Iterators.take(Iterators.drop(es, n-1), k-1)), lowerbound)
    ub = min(sum(Int(round(weight(e))) for e in Iterators.drop(es, n-(k-1))), upperbound)
    @debug "Objective has $(string(CYAN_FG(string(lb)))) as natural lowerbound."
    @debug "Objective has $(string(MAGENTA_FG(string(ub)))) as natural upperbound."
    sort!(es, by=dst, alg=MergeSort)
    sort!(es, by=src, alg=MergeSort)
    @variables(model,begin
        y[i=1:n, j=1:n; has_edge(graph, i, j)], (binary=true, lower_bound=0, upper_bound=1) # Set bounds, that solver doesn't have to set it implicitly
        x[e=1:m], (binary=true, lower_bound=0, upper_bound=1)
        o, (integer=true, lower_bound=lb, upper_bound=ub)
    end)

    for i in 1:m
        e = es[i]
        s = src(e)
        d = dst(e)
        if s != 1 && d != 1
            @constraint(model, y[s, d] + y[d, s] == x[i])
        end
    end

    @constraints(model, begin
        sum(y[i,j] for i = 2:n for j=2:n if has_edge(graph, i ,j)) == k-1 # there are k-1 edges in a tree of size k
        sum(x[i] for i in 1:m if src(es[i]) != 1 && dst(es[i]) != 1) == k-1
        sum(x[e] * Int(round(Int, weight(es[e]))) for e in 1:m) == o
    end)
    @objective(model, Min, o)
end

function basic_kmst_warmstart!(model, graph :: SimpleWeightedGraph, k :: Int, solution :: Vector{Edge{Int}})
    n :: Int = nv(graph)
    m :: Int = ne(graph)
    es :: Vector{SimpleWeightedEdge} = collect(edges(graph))
    sort!(es, by=dst)
    sort!(es, by=src, alg=MergeSort)
    x :: Vector{Int} = zeros(Int, m)
    totalweight :: Int = 0
    for s in solution
        for r in 1:m
            e = es[r]
            i = src(e)
            j = dst(e)
            w = weight(e)
            if (src(s) == i && dst(s) == j) || (dst(s) == i && src(s) == j)
                x[r] = 1
                totalweight += Int(round(Int ,w))
                break
            end
        end
    end
    for r in 1:m
        e = es[r]
        s = src(e)
        d = dst(e)
        if src(e) != 1 && dst(e) != 1
            set_start_value(variable_by_name(model, "x[$r]"), x[r])
            y_ij = start_value(variable_by_name(model, "y[$s,$d]"))
            y_ji = start_value(variable_by_name(model, "y[$d,$s]"))
            if !isnothing(y_ij) && !isnothing(y_ji) && y_ij + y_ji != x[r]
                @warn "Start values for constraint y[$s,$d] + y[$d,$s] = x[$r] do not hold."
            end
        end
    end
    set_start_value(variable_by_name(model, "o"), totalweight)
end

struct KMSTSolution
    termination_status :: MOI.TerminationStatusCode
    solution_graph :: Union{SimpleWeightedGraph, Nothing}
    objective_value :: Float64
    relative_gap :: Float64
    solve_time_sec :: Float64
end

function solve!(model, graph :: SimpleWeightedGraph, k :: Int) :: KMSTSolution
    ts :: MOI.TerminationStatusCode = MOI.INTERRUPTED
    try
        @debug "Starting to optimize."
        optimizationtime = @elapsed optimize!(model)
        @debug "Finished optimizing in $(format_seconds_readable(optimizationtime, 2))."
        ts = termination_status(model)
    catch e
        @warn "Exception while optimizing. $e"
    end
    if ts != MOI.OPTIMAL
        warning = "Optimizer did not find an optimal solution."
        if ts == MOI.TIME_LIMIT
            warning *= " Time limit of $(format_seconds_readable(time_limit_sec(model))) reached."
        elseif ts == MOI.INFEASIBLE
            warning *= " Problem is infeasible."
            compute_conflict!(model)
            cs = MOI.get(model, MOI.ConflictStatus())
            if MOI.CONFLICT_FOUND

            end
        else
            warning *= " $ts"
        end
        @warn warning
    end

    st :: Float64 = solve_time(model)
    if !has_values(model)
        return KMSTSolution(ts, nothing, Inf64, Inf64, st)
    end
    rg :: Float64 = relative_gap(model)
    n::Int = nv(graph)
    solution_graph::SimpleWeightedGraph = get_solution_graph(model, graph)
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
        cl::Int = length(component)
        if cl != 1 && cl != k
            @error "Solution is not connected." maxlog=1
            @debug "Should only have components with size 1 or $k but found with size $cl: $(component)."
        end
    end
    ov::Float64 = objective_value(model)
    ovrounded::Int = Int(round(ov))
    w::Int = sum(Int(round(Int,weight(e))) for e in edges(solution_graph))
    if w != ovrounded
        @warn "MILP solver claims objective value is $ovrounded ($ov), however got $(w)."
    end
    return KMSTSolution(ts, solution_graph, ov, rg, st)
end

function get_solution_graph(model, graph :: SimpleWeightedGraph)
    n :: Int = nv(graph)
    solution_graph = SimpleWeightedGraph(n)
    es :: Vector{SimpleWeightedEdge} = collect(edges(graph))
    sort!(es, by=weight)
    sort!(es, by=dst, alg=MergeSort)
    sort!(es, by=src, alg=MergeSort)
    for e in es
        i = src(e)
        j = dst(e)
        w = weight(e)
        if i > 1 && j > 1
            x_ij = variable_by_name(model, "y[$i,$j]")
            xv = Int(round(value(x_ij)))
            if xv == 1
                add_edge!(solution_graph, i, j, w)
            else
                x_ji = variable_by_name(model, "y[$j,$i]")
                xv = Int(round(value(x_ji)))
                if xv == 1
                    add_edge!(solution_graph, i, j, w)
                end
            end
        end
    end
    
    return solution_graph
end