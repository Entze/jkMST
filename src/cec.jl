using JuMP
using LightGraphs, SimpleWeightedGraphs
using Random

lazy_cec_constraints = 0

function cycle_elimination_constraints!(
    model,
    graph::SimpleWeightedGraph,
    k::Int)

    global lazy_cec_constraints = 0
    @constraint(model, sum(variable_by_name(model, "y[1,$j]") for j in 2:nv(graph)) == 1)
    @constraint(model, sum(variable_by_name(model, "y[$j,1]") for j in 2:nv(graph)) >= 1)
    @constraint(model, sum(variable_by_name(model, "y[$j,1]") for j in 2:nv(graph)) <= k-1)
    for i in 2:nv(graph)
        @constraint(model, sum(variable_by_name(model, "y[$i,$j]")
                               for j in 1:nv(graph) if has_edge(graph, i, j))
                    >= sum(variable_by_name(model, "y[$j,$i]")
                           for j in 1:nv(graph) if has_edge(graph, j, i)))

        @constraint(model, sum(variable_by_name(model, "y[$i,$j]")
                               for j in 1:nv(graph) if has_edge(graph, i, j))
                    <= sum(variable_by_name(model, "y[$j,$i]")
                           for j in 1:nv(graph) if has_edge(graph, j, i)) * (k-1))
        @constraint(model, sum(variable_by_name(model, "y[$i,$j]") for j in 1:nv(graph) if has_edge(graph, i,j)) <= k-1)
        @constraint(model, sum(variable_by_name(model, "y[$j,$i]") for j in 1:nv(graph) if has_edge(graph, i,j)) <= 1)
    end

    MOI.set(model,
            MOI.LazyConstraintCallback(),
            cec_lazy_clause_generator(model, graph, k))
end

function cec_lazy_clause_generator(model, graph::SimpleWeightedGraph, k::Int)
    es = collect(SimpleWeightedEdge, edges(graph))
    sort!(es, by=dst)
    sort!(es, by=src, alg=MergeSort)
    x::Dict{Int,VariableRef} = Dict{Int,VariableRef}()
    y::Dict{Tuple{Int,Int},VariableRef} = Dict{Tuple{Int,Int},VariableRef}()
    r::Int = 0
    for e in es
        r += 1
        i::Int = src(e)
        j::Int = dst(e)
        x[r] = variable_by_name(model, "x[$r]")
        y[i,j] = variable_by_name(model, "y[$i,$j]")
        y[j,i] = variable_by_name(model, "y[$j,$i]")
    end
    return function (cb_data)
        status = callback_node_status(cb_data, model)
        if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
            @debug "Found fractional solution"
            return
        elseif status == MOI.CALLBACK_NODE_STATUS_INTEGER
            @debug "Found integer solution"
            # x_val::Dict{Int, } = Dict{Int, Float64}()
            y_val::Dict{Tuple{Int,Int},Int} = Dict{Tuple{Int,Int},Int}()
            r = 0
            for e in es
                r += 1
                i::Int = src(e)
                j::Int = dst(e)
                y_val[i,j] = round(callback_value(cb_data, y[i,j]))
                y_val[j,i] = round(callback_value(cb_data, y[j,i]))
            end
            cycle = find_cycle(graph, y_val,k)
            if !isnothing(cycle)
                cyclelength::Int = length(cycle)
                @debug "Detected cycle with length $cyclelength: $cycle."
                constraint = @build_constraint(sum(y[i,j]
                                                    for i in cycle
                                                        for j in cycle
                                                            if has_edge(graph, i, j))
                                               <= cyclelength - 1)

                # for debugging:
                # arcs::Vector{String} = Vector{String}()
                # for i in cycle
                #    for j in cycle
                #        if has_edge(graph, i, j)
                #            push!(arcs, "y[$i,$j]")
                #        end
                #    end
                # end
                # @debug join(arcs, " + ") * " <= $(cyclelength - 1)"

                global lazy_cec_constraints += 1
                MOI.submit(model, MOI.LazyConstraint(cb_data), constraint)
            end
        else
            @assert status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
            @debug "Found either fractional or integer solution"
        end
    end
end

function find_cycle(
    graph::SimpleWeightedGraph,
    y_val::Dict{Tuple{Int,Int},Int64},
    k::Int)::Union{Nothing,Set{Int}}
    visited::Set{Int} = Set{Int}()
    queue::Set{Int} = Set{Int}([1])
    while !isempty(queue)
        @assert !isempty(queue) "queue should not be empty."
        c::Int = pop!(queue)
        ns = neighbors(graph, c)
        for n::Int in ns
            if y_val[c, n] == 1 && n != 1
                push!(queue, n)
            end
        end
        push!(visited, c)
    end
    if length(visited) == k + 1
        return # no cycles
    end
    @assert length(visited) <= k "Solution should be a forest, but only has one tree."

    to_check::Set{Int} = Set{Int}(2:nv(graph))

    filter!(v -> !(v in visited), to_check)
    @assert length(to_check) == (nv(graph) - 1) - (length(visited) - 1) "Should have pruned $(length(visited) - 1) of $(nv(graph) - 1) items. However only $(length(to_check)) left."

    while !isempty(to_check)
        @assert isempty(queue) "queue should be empty."
        push!(queue, pop!(to_check))
        visited = Set{Int}()
        @assert isempty(visited) "visited should be empty."
        while !isempty(queue)
            @assert !isempty(queue) "queue should not be empty."
            c::Int = pop!(queue)
            # @debug "Visiting node $c"
            push!(visited, c)
            ns = neighbors(graph, c)
            for n::Int in ns
                if y_val[c, n] == 1 && n != 1
                    if n in visited
                        return visited
                    end
                    push!(queue, n)
                end
            end
        end
        filter!(v -> !(v in visited), to_check)
    end

    @warn "No cycle detected. Something is wrong with the implementation."

end
