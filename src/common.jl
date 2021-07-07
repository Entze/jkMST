using Logging
using JuMP, CPLEX
using LightGraphs, SimpleWeightedGraphs
using PrettyTables, Formatting, ProgressMeter, Crayons.Box

function value_rounded(v)
    return round(Int, value(v))
end

function basic_kmst!(model,
        graph :: SimpleWeightedGraph,
        k :: Int;
        lowerbound :: Int = typemin(Int),
        upperbound :: Int = typemax(Int),
        isdebug=false)
    n :: Int = nv(graph)
    m :: Int = ne(graph)
    progress::Progress = Progress(1 + m + 2 * m + m + 3 + 1,
                                    dt=0.5,
                                    barglyphs=BarGlyphs('[', '#', ['-','~','+','*','=','>','#'], ' ', ']'),
                                    enabled=isdebug)
    es = collect(SimpleWeightedEdge, edges(graph))
    sort!(es, by=weight)

    lb::Int = sum(round(Int, weight(e)) for e in Iterators.take(Iterators.drop(es, n-1), k-1))
    if lb > lowerbound
        @debug "Objective has $(string(CYAN_FG(string(lb)))) as natural lowerbound."
    end
    lb = max(lb, lowerbound)
    ub::Int = sum(round(Int, weight(e)) for e in Iterators.drop(es, n-(k-1)))
    if ub < upperbound
        @debug "Objective has $(string(MAGENTA_FG(string(ub)))) as natural upperbound."
    end
    ub = min(ub, upperbound)
    sort!(es, by=dst, alg=MergeSort)
    sort!(es, by=src, alg=MergeSort)
    @variable(model, o, Int, lower_bound=lb, upper_bound=ub)
    next!(progress)
    @variable(model, x[e=1:m], Bin, lower_bound=0, upper_bound=1)
    ProgressMeter.update!(progress, m)
    y::Dict{Tuple{Int,Int}, VariableRef} = Dict{Tuple{Int,Int}, VariableRef}()
    for e in es
        i::Int = src(e)
        j::Int = dst(e)
        y[i,j] = @variable(model, base_name="y[$i,$j]", binary=true, lower_bound=0, upper_bound=1) # Set bounds, that solver doesn't have to set it implicitly
        next!(progress)
        y[j,i] = @variable(model, base_name="y[$j,$i]", binary=true, lower_bound=0, upper_bound=1)
        next!(progress)
    end

    for i in 1:m
        e = es[i]
        s = src(e)
        d = dst(e)
        @constraint(model, y[s, d] + y[d, s] == x[i]) # Constraint 2
        next!(progress)
    end

    @constraints(model, begin
        sum(x[i] for i in 1:m if src(es[i]) != 1 && dst(es[i]) != 1) == k-1 # Constraint 3
        sum(x[e] * Int(round(Int, weight(es[e]))) for e in 1:m) == o # Constraint 4
    end)
    ProgressMeter.update!(progress, 3)
    @objective(model, Min, o)
    next!(progress)
    finish!(progress)
    #MOI.set(model, MOI.HeuristicCallback(), heuristic!(model, graph, k))
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
        end
    end
    set_start_value(variable_by_name(model, "o"), totalweight)
end

function check_kmst_warmstart!(model, graph :: SimpleWeightedGraph, k :: Int) :: Int
    n :: Int = nv(graph)
    m :: Int = ne(graph)
    es :: Vector{SimpleWeightedEdge} = collect(edges(graph))
    sort!(es, by=dst)
    sort!(es, by=src, alg=MergeSort)
    x :: Vector{Union{Nothing,Int}} = fill(nothing, m)
    y :: Matrix{Union{Nothing,Int}} = fill(nothing, n, n)
    filled :: Int = 0
    totalweight :: Int = 0
    for r in 1:m
        e = es[r]
        i = src(e)
        j = dst(e)
        x[r] = start_value(variable_by_name(model, "x[$r]"))
        if x[r] == 1
            totalweight += Int(round(Int, weight(e)))
        end
        y[i, j] = start_value(variable_by_name(model, "y[$i,$j]"))
        y[j, i] = start_value(variable_by_name(model, "y[$j,$i]"))
        if isnothing(y[i,j]) && !isnothing(y[j,i]) && !isnothing(x[r])
            y[i,j] = x[r] - y[j,i]
            set_start_value(variable_by_name(model, "y[$i,$j]"), y[i,j])
            filled += 1
        elseif !isnothing(y[i,j]) && isnothing(y[j,i]) && !isnothing(x[r])
            y[j,i] = x[r] - y[i,j]
            set_start_value(variable_by_name(model, "y[$j,$i]"), y[j,i])
            filled += 1
        elseif !isnothing(y[i,j]) && !isnothing(y[j,i]) && isnothing(x[r])
            if 0 <= y[i,j] + y[j,i] <= 1
                x[r] = y[i,j] + y[j,i]
                set_start_value(variable_by_name(model, "x[$r]"), x[r])
                filled += 1
            else
                @warn "0 > y[$i,$j] + y[$j,$i] > 1 | ($(y[i,j]) + $(y[j,i]))"
            end
        elseif !isnothing(y[i,j]) && !isnothing(y[j,i]) && !isnothing(x[r]) && x[r] != y[i,j] + y[j,i]
            @warn "x[$r] != y[$i,$j] + y[$j,$i] | ($(x[r]) != $(y[i,j]) + $(y[j,i])"
        end
    end

    return filled

end

function heuristic!(model,
                    graph :: SimpleWeightedGraph,
                    k :: Int)
    y::Dict{Tuple{Int, Int}, VariableRef} = Dict{Tuple{Int, Int}, VariableRef}()
    y_val::Dict{Tuple{Int, Int}, Float64} = Dict{Tuple{Int, Int}, Float64}()
    x::Dict{Int, VariableRef} = Dict{Int, VariableRef}()
    x_val::Dict{Int, Float64} = Dict{Int, Float64}()
    o::VariableRef = variable_by_name(model, "o")
    o_val::Float64 = 0
    es = collect(edges(graph))
    sort!(es, by=dst)
    sort!(es, by=src, alg=MergeSort)
    r::Int = 0
    for e in es
        r += 1
        x[r] = variable_by_name(model, "x[$r]")
        i::Int = src(e)
        j::Int = dst(e)
        y[i,j] = variable_by_name(model, "y[$i,$j]")
        y[j,i] = variable_by_name(model, "y[$j,$i]")
    end
    return function(cb_data)
        variables::Vector{VariableRef} = Vector{VariableRef}()
        values::Vector{Float64} = Vector{Float64}()
        status = callback_node_status(cb_data, model)
        o_val = callback_value(cb_data, o)
        for ((i,j), v) in y
            y_val[i,j] = callback_value(cb_data, v)
        end
        for (k, v) in x
            x_val[k] = callback_value(cb_data, v)
        end
        if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
            @debug "Found fractional solution. Attempting to round."
            arcs::Vector{Tuple{Int,Int,Int}} = Vector{Tuple{Int,Int,Int}}()
            r = 0
            for e in es
                r += 1
                i::Int = src(e)
                j::Int = dst(e)
                push!(arcs, (r,i,j))
                push!(arcs, (r,j,i))
            end

            sort!(arcs, by=(a -> a[3]))
            sort!(arcs, by=(a -> a[2]), alg=MergeSort)
            sort!(arcs, by=(a -> weight(es[a[1]])), alg=MergeSort)
            sort!(arcs, by=(a -> y_val[a[2],a[3]]), alg=MergeSort, rev=true)
            sort!(arcs, by=(a -> a[2] == 1 || a[3] == 1), alg=MergeSort)

            visited::Vector{Bool} = zeros(Bool, nv(graph))
            parents::Vector{Int} = zeros(Int, nv(graph))
            lastnode::Int = 1
            treesize::Int = 1
            visited[1] = true
            parents[1] = -1

            includededges::Vector{Bool} = zeros(Bool, length(es))

            while treesize < k+1
                for (r,i,j) in arcs
                    (visited[i] && !visited[j]) || continue
                    visited[j] = true
                    includededges[r] = true
                    parents[i] = lastnode
                    parents[j] = i
                    lastnode = j
                    treesize += 1
                    @debug "Including $r ($i -> $j)"
                    break
                end
                filter!(v -> !visited[v[3]], arcs)
            end
            @assert count(visited) == k+1 "Tree should have size $k, but has size $(count(visited))."

            push!(variables, o)
            push!(values, 0.0)

            r = 0
            for e in es
                r += 1
                if includededges[r]
                    values[1] += round(weight(e))
                    push!(variables, x[r])
                    push!(values, 1.0)
                end
            end

            # for ((i,j),v) in y
            #    push!(variables, v)
            #    push!(values, Float64(visited[i] && visited[j] && (parents[j] == i || parents[i] == j)))
            # end

            status = MOI.submit(model, MOI.HeuristicSolution(cb_data), variables, values)

            debugstring = "Found solution with weight $(round(Int, values[1])), "
            if status == MOI.HEURISTIC_SOLUTION_ACCEPTED
                debugstring *= string(GREEN_FG("accepted"))
            elseif status == MOI.HEURISTIC_SOLUTION_REJECTED
                debugstring *= string(RED_FG("rejected"))
            else
                @assert status == MOI.HEURISTIC_SOLUTION_UNKNOWN
                debugstring *= string(YELLOW_FG("unknown status"))
            end
            @debug debugstring * "."

        elseif status == MOI.CALLBACK_NODE_STATUS_INTEGER
            @debug "Found integer solution. Attempting to permute."
        else
            @assert status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
        end
    end
end


struct KMSTSolution
    termination_status :: MOI.TerminationStatusCode
    solution_graph :: Union{SimpleWeightedGraph, Nothing}
    objective_value :: Float64
    relative_gap :: Float64
    solve_time_sec :: Float64
    bnb_nodes :: Union{Int, Nothing}
    violated_ineq :: Union{Int, Nothing}
end

function print_variable_table(model, graph :: SimpleWeightedGraph)
    if !has_values(model)
        @debug "Cannot print values, model has no values."
        return
    end
    @assert has_values(model) "Model should have values."
    n::Int= nv(graph)
    m::Int= ne(graph)
    es::Vector{SimpleWeightedEdge}= collect(edges(graph))
    sort!(es, by=dst)
    sort!(es, by=src, alg=MergeSort)

    all_vars::Vector{VariableRef} = all_variables(model)

    different_variables::Int = 0

    header :: Vector{String} = ["", "1-Vars"]
    names :: Vector{String} = []
    multidimvarslines :: Dict{String, Int} = Dict()
    multidimvarsdims :: Dict{String, Int} = Dict()
    multilines::Int = 0
    onevarsline::Int = 1
    onevarslines::Int = 0

    for var in all_vars
        name::String = JuMP.name(var)

        delim::Union{Int, Nothing} = findfirst('[', name)

        if isnothing(delim)
            onevarslines += 1
        elseif delim > 1
            basename::String = name[1:(delim-1)]
            if basename in keys(multidimvarsdims)
                @assert basename in keys(multidimvarslines) "If basename $basename is in multidimvarsdims it should be in multidimvarslines."
                @assert basename in names "If basename $basename is in multidimvarsdims it should be in names."
                multidimvarslines[basename] += 1
            else
                @assert !(basename in keys(multidimvarslines)) "If basename $basename is not in multidimvarsdims it should not be in multidimvarslines."
                @assert !(basename in names) "If basename $basename is in multidimvarsdims it should not be in names."
                multidimvarslines[basename] = 1
                multidimvarsdims[basename] = count(v -> v == ',', name) + 1
                push!(header, "")
                push!(header, basename)
                push!(names, basename)
            end
            @assert basename in keys(multidimvarslines) "Basename $basename should be present in multidimvarslines"
            @assert basename in keys(multidimvarsdims) "Basename $basename should be present in multidimvarsdims"
        end
    end
    if !isempty(multidimvarslines)
        multilines = maximum(values(multidimvarslines))
    end

    width::Int = length(header)
    height::Int = max(onevarslines, multilines)

    data::Matrix{Union{String, Int, Float64}} = fill("", height, width)

    onevarslines = 1
    for k in keys(multidimvarslines)
        multidimvarslines[k] = 1
    end

    for var in all_vars
        name :: String = JuMP.name(var)
        delim::Union{Int, Nothing} = findfirst('[', name)
        if isnothing(delim)
            @assert onevarslines <= height "1-vars exceeds height."
            data[onevarslines, 1] = name
            data[onevarslines, 2] = value(var)
            onevarslines += 1
        elseif delim > 1
            if round(value(var)) != 0
                basename::String = name[1:(delim-1)]
                @assert basename in keys(multidimvarslines) "Basename $basename never seen before. Not found in multidimvarslines."
                dims::String = name[delim:end]
                @assert basename in names "Basename $basename never seen before. Not found in names."
                index::Int = filter(v -> v[2] == basename, collect(enumerate(names)))[1][1]
                @assert 1 <= index * 2 + 1 <= width "$basename index ($index) exceeds width."
                @assert 1 <= multidimvarslines[basename] <= height "$name exceeds height."
                data[multidimvarslines[basename], index * 2 + 1] = dims
                if basename == "x"
                    e = es[parse(Int, dims[2:(end-1)])]
                    data[multidimvarslines[basename], index * 2 + 1] *= " = [$(src(e)),$(dst(e))]"
                end
                data[multidimvarslines[basename], index * 2 + 2] = value(var)
                multidimvarslines[basename] += 1
            end
        end
    end

    if !isempty(multidimvarslines)
        multilines = maximum(values(multidimvarslines))
    end
    height = max(onevarslines, multilines)

    data_view = view(data, 1:(height-1), :)

    if height > 45
        pretty_table(data_view, header=header, backend=Val(:html))
    else
        pretty_table(data_view, header=header)
    end

end

function solve!(model,
                graph :: SimpleWeightedGraph,
                k :: Int,
                mode :: SolvingMode) :: KMSTSolution
    ts :: MOI.TerminationStatusCode = MOI.INTERRUPTED
    @debug "Starting to optimize."
    optimizationtime = @elapsed optimize!(model)
    @debug "Finished optimizing in $(format_seconds_readable(optimizationtime, 2))."
    ts = termination_status(model)
    st :: Float64 = solve_time(model)
    if ts != MOI.OPTIMAL
        warning = "Optimizer did not find an optimal solution."
        if ts == MOI.TIME_LIMIT
            tl::Float64 = time_limit_sec(model)
            warning *= " Time limit of $(format_seconds_readable(tl)) reached"
            if tl < st
                warning *= " after $(format_seconds_readable(st)), overdrawn by $(format_seconds_readable(st - tl))"
            end
            warning *= "."
        elseif ts == MOI.MEMORY_LIMIT
            warning *= " Memory limit reached."
        else
           if ts == MOI.INFEASIBLE && mode != cec && mode != dcc
              if lowercase(solver_name(model)) != "cplex"
                  @debug "Setting optimizer from $(solver_name(model)) to CPLEX."
                  set_optimizer(model, CPLEX.Optimizer)
              end
              compute_conflict!(model)

              constraint_types::Vector{Tuple{DataType, DataType}} = list_of_constraint_types(model)

              debugstring::String = "Following constraints are part of the conflict:"

              for (ref, typ) in constraint_types
                  !(ref == VariableRef && typ == MOI.Integer) || continue
                  !(ref == VariableRef && typ == MOI.ZeroOne) || continue
                  !(ref == GenericAffExpr{Float64, VariableRef} && typ == MOI.Interval{Float64}) || continue
                  all_cons::Vector{ConstraintRef} = all_constraints(model, ref, typ)
                  @debug "$ref: $typ"
                  for con in all_cons
                      if MOI.get(model.moi_backend, MOI.ConstraintConflictStatus(), con.index) != MOI.NOT_IN_CONFLICT
                          debugstring *= "\n$con"
                      end
                  end
              end
              warning *= " Problem is infeasible."
              @debug debugstring
          else
              warning *= " $ts"
          end
          print_variable_table(model, graph)
        end
        @warn warning
    end

    bnbn::Union{Nothing, Int} = nothing
    if lowercase(solver_name(model)) == "cplex"
        bnbn = node_count(model)
    end
    lc::Union{Nothing, Int} = nothing
    if mode == cec
        lc = lazy_cec_constraints
    end
    if mode == dcc
        lc = lazy_dcc_constraints
    end
    if !has_values(model)
        return KMSTSolution(ts, nothing, Inf64, Inf64, st, bnbn, lc)
    end

    rg::Float64 = relative_gap(model)
    n::Int = nv(graph)
    solution_graph::SimpleWeightedGraph = get_solution_graph(model, graph)
    components = connected_components(solution_graph)
    bug::Bool = false
    if ne(solution_graph) != k-1
        @error "Solution has wrong number of edges."
        @debug "Should have $(k - 1) edges but has $(ne(solution_graph))."
        bug = true
    end
    if n - k + 1 != length(components)
        @error "Solution has wrong number of components."
        @debug "Should have $(n - k + 1) components but has $(length(components))."
        bug = true
    end
    for component in components
        cl::Int = length(component)
        if cl != 1 && cl != k
            @error "Solution is not connected." maxlog=1
            @debug "Should only have components with size 1 or $k but found with size $cl: $(component)."
            bug = true
        end
    end
    ov::Float64 = objective_value(model)
    ovrounded::Int = Int(round(ov))
    w::Int = sum(Int(round(Int,weight(e))) for e in edges(solution_graph))
    if w != ovrounded
        @warn "MILP solver claims objective value is $ovrounded ($ov), however got $(w)."
        bug = true
    end
    if bug
        print_variable_table(model, graph)
    end
    return KMSTSolution(ts, solution_graph, ov, rg, st, bnbn, lc)
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
