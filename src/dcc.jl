using JuMP
using LightGraphs, SimpleWeightedGraphs
using Random

lazy_dcc_constraints = 0

function directed_cutset_constraints!(
    model,
    graph::SimpleWeightedGraph,
    k::Int)
    global lazy_dcc_constraints = 0
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
        # @constraint(model, sum(variable_by_name(model, "y[$i,$j]") for j in 1:nv(graph) if has_edge(graph, i,j)) <= k-1)
        @constraint(model, sum(variable_by_name(model, "y[$j,$i]") for j in 1:nv(graph) if has_edge(graph, i,j)) <= 1)
    end
    MOI.set(model,
            MOI.LazyConstraintCallback(),
            dcc_lazy_clause_generator(model, graph, k))
end


function dcc_lazy_clause_generator(model, graph::SimpleWeightedGraph, k::Int)
    es = collect(SimpleWeightedEdge, edges(graph))
    sort!(es, by=dst)
    sort!(es, by=src, alg=MergeSort)
    y::Dict{Tuple{Int,Int},VariableRef} = Dict{Tuple{Int,Int},VariableRef}()
    for e in es
        i::Int = src(e)
        j::Int = dst(e)
        y[i,j] = variable_by_name(model, "y[$i,$j]")
        y[j,i] = variable_by_name(model, "y[$j,$i]")
    end
    return function(cb_data)
        status = callback_node_status(cb_data, model)
        if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
            @debug "Found fractional solution"
            return
        elseif status == MOI.CALLBACK_NODE_STATUS_INTEGER
            @debug "Found integer solution"
            y_val::Dict{Tuple{Int,Int},Int} = Dict{Tuple{Int,Int},Int}()
            for e in es
                i::Int = src(e)
                j::Int = dst(e)
                y_val[i,j] = round(callback_value(cb_data, y[i,j]))
                y_val[j,i] = round(callback_value(cb_data, y[j,i]))
            end
            cutset = find_cutset(graph, y_val, k)
            if isnothing(cutset)
                return
            end
            @debug "Found cutset of size $(length(cutset)): $cutset."
            complement::Set{Int} = Set{Int}(2:nv(graph))
            filter!(v -> !(v in cutset), complement)
            @assert length(complement) + length(cutset) == nv(graph) "Not all nodes in cutset and complement"
            filter!(v -> v != 1, cutset)
            constraint = @build_constraint(sum(y[i,j] for i in cutset for j in complement if has_edge(graph, i, j)) >= sum(y[1, i] for i in cutset)) # Constraint 39
            global lazy_dcc_constraints += 1
            MOI.submit(model, MOI.LazyConstraint(cb_data), constraint)
        end
    end
end

function find_cutset(
    graph::SimpleWeightedGraph,
    y_val::Dict{Tuple{Int,Int},Int64},
    k::Int)::Union{Nothing,Set{Int}}
    visited::Set{Int} = Set{Int}()
    queue::Set{Int} = Set{Int}([1])
    while !isempty(queue)
        @assert !isempty(queue) "queue should not be empty."
        c::Int = pop!(queue)
        push!(visited, c)
        ns = neighbors(graph, c)
        for n::Int in ns
            if y_val[c, n] == 1 && n != 1
                push!(queue, n)
            end
        end
    end
    if length(visited) == k+1
        return # no cutset needed
    end
    @assert length(visited) <= k "Solution should be a forest, but only has one tree."
    return visited

end
