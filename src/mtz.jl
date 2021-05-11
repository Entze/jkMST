
using JuMP
using LightGraphs, SimpleWeightedGraphs

function miller_tucker_zemlin(model, graph :: SimpleWeightedGraph, k :: Int)
    n::Int = nv(graph)
    @variables(model, begin
        u[1:n], (integer=true, upper_bound = k+1)
        d[2:n], (integer=true, lower_bound = 0, upper_bound = k)
        z[2:n], (binary=true, lower_bound = 0, upper_bound = 1)
    end)
    omega :: Int = 1
    fix(u[omega], 0, force=true) # Constraint (6)
    for i in (omega+1):n
        set_lower_bound(u[i], 1)
        for j in 1:n
            if has_edge(graph, i, j)
                y_ij = variable_by_name(model, "y[$i,$j]") 
                @constraint(model,
                    u[i] >= u[j] + y_ij - k * (1 - y_ij) # Constraint (7)
                )
            end
        end
    end

    for i in (omega+1):n
        @constraint(model, sum(variable_by_name(model, "y[$i,$j]") for j = 1:n if has_edge(graph, i, j)) == 1) # Constraint (8)
    end
    @costraint(model, sum(variable_by_name(model, "y[$omega,$j]") for j = (omega + 1):n) == 0) # Constraint (9)

    for j in (omega+1):n
        @constraint(model,
            sum(variable_by_name(model, "y[$i,$j]") for i in 2:n if has_edge(graph, i, j)) 
            + sum(variable_by_name(model, "y[$j,$k_]") for k_ in 2:n if has_edge(graph, j, k_)) == d[j]) # Constraint (10)
        @constraints(model, begin
            z[j] * k >= d[j] # Constraint (11)
            z[j] <= d[j] # Constraint (12)
        end)
    end

    @constraint(model, sum(z[i] for i = 2:n) == k) # Constraint (13)

end

function miller_tucker_zemlin_warmstart(model, graph :: SimpleWeightedGraph, k :: Int, solution)
    n :: Int = nv(graph)
    z :: Vector{Int} = zeros(Int, n)
    d :: Vector{Int} = zeros(Int, n)
    connected :: Vector{Vector{Int}} = [[] for i=1:n]
    for e in solution
        i = src(e)
        j = dst(e)
        z[i] = 1
        z[j] = 1
        push!(connected[i], j)
        push!(connected[j], i)
        d[i] += 1
        d[j] += 1
    end
    s = nothing
    for c in 2:n
        if !isempty(connected[c])
            s = c
            break
        end
    end

    stack = [s]
    marked :: Vector{Bool} = zeros(Bool, n)
    dagsort :: Vector{Int} = []
    while !isempty(stack)
        node = pop!(stack)
        neighbours = connected[node]
        for n in neighbours
            if !marked[n]
                push!(stack, n)
                set_start_value(variable_by_name(model, "y[$n,$node]"), 1)
                set_start_value(variable_by_name(model, "y[$node,$n]"), 0)
            end
        end
        marked[node] = true
        pushfirst!(dagsort, node)
    end

    for i in 2:n
        set_start_value(variable_by_name(model, "z[$i]"), z[i])
        set_start_value(variable_by_name(model, "d[$i]"), d[i])
        if z[i] == 0
            set_start_value(variable_by_name(model, "y[$i,1]"), 1)
            for j in 2:n
                if has_edge(graph, i, j)
                    set_start_value(variable_by_name(model, "y[$i,$j]"), 0)
                    set_start_value(variable_by_name(model, "y[$j,$i]"), 0)
                end
            end
        end
    end
end

function check_miller_tucker_zemlin_warmstart(model, graph :: SimpleWeightedGraph, k :: Int) :: Int
    n :: Int = nv(graph)
    m :: Int = ne(graph)
    es :: Vector{SimpleWeightedEdge} = collect(edges(graph))
    sort!(es, by=dst)
    sort!(es, by=src, alg=MergeSort)
    x :: Vector{Union{Nothing, Int}} = fill(nothing, m)
    y :: Matrix{Union{Nothing, Int}} = fill(nothing, n, n)
    z :: Vector{Union{Nothing, Int}} = fill(nothing, n)
    u :: Vector{Union{Nothing, Int}} = fill(nothing, n)
    d :: Vector{Union{Nothing, Int}} = fill(nothing, n)
    filled :: Int = 0
    totalweight :: Int = 0
    for r in 1:m
        e::SimpleWeightedEdge = es[r]
        i::Int = src(e)
        j::Int = dst(e)
        x[r] = start_value(variable_by_name(model, "x[$r]"))
        if x[r] == 1
            totalweight += Int(round(Int,weight(e)))
        end
        for v in [i,j]
            if v != 1
                if isnothing(z[v])
                    z[v] = start_value(variable_by_name(model, "z[$v]"))
                end
                if isnothing(u[v])
                    u[v] = start_value(variable_by_name(model, "u[$v]"))
                end
                if isnothing(d[v])
                    d[v] = start_value(variable_by_name(model, "d[$v]"))
                end

            end
        end
        y[i, j] = start_value(variable_by_name(model, "y[$i,$j]"))
        y[j, i] = start_value(variable_by_name(model, "y[$j,$i]"))

        for v in [i,j]
            if v != 1
                if !isnothing(d[v]) && !isnothing(z[v])
                    if (d[v] >= 1 && z[v] == 0) || (d[v] == 0 && z[v] == 1)
                        @warn "z[$v] != sgn(d[$v]) | ($(z[v]) != sgn($(d[v])))"
                    end
                elseif !isnothing(d[v]) && isnothing(z[v])
                    if d[v] == 0
                        z[v] = 0
                    else
                        z[v] = 1
                    end
                    set_start_value(variable_by_name(model, "z[$v]"), z[v])
                    filled += 1
                elseif isnothing(d[v]) && !isnothing(z[v])
                    if z[v] == 0
                        d[v] = 0
                        set_start_value(variable_by_name(model, "d[$v]"), d[v])
                    end
                end
            end
        end
    end

    for i in 2:n
        if !isnothing(u[i])
            for j in 2:n
                if !isnothing(u[j]) && !isnothing(y[i,j]) && has_edge(graph, i, j)
                    if u[i] < u[j] + y[i,j] - k * (1 - y[i,j])
                        @warn "u[$i] < u[$j] + y[$i,$j] - k * (1 - y[$i,$j]) | ($(u[i]) < $(u[j]) + $(y[i,j]) - $deg * (1 - $(y[i,j])))"
                    end
                end
            end
        end
    end

    return filled
end