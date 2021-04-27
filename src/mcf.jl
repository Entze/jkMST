
using JuMP
using LightGraphs, SimpleWeightedGraphs


function multi_commodity_flow!(model, graph :: SimpleWeightedGraph, k :: Int; generate_timeout_sec::Float64 = Inf64, isdebug=false)
    n::Int = nv(graph)
    nes::Int = ne(graph)
    es::Vector{SimpleWeightedEdge} = collect(edges(graph))
    sort!(es, by=dst)
    sort!(es, by=src, alg=MergeSort)

    progress::Progress = Progress((n-1) * nes + n - 1 + 2 * (n-1) + 1 + (n-1) * nes + (n-1) * (2 + (n-1)) + 1 ,
                                    dt=0.5,
                                    barglyphs=BarGlyphs('[', '#', ['-','~','+','*','=','>','#'], ' ', ']'),
                                    enabled=isdebug)
    generate_time::Float64 = 0

    generate_time += @elapsed begin
        @variable(model, z[2:n], Bin, lower_bound=0, upper_bound=1)
    end
    ProgressMeter.update!(progress, n-1)
    f::Dict{Tuple{Int, Int, Int},VariableRef} = Dict()
    for m in 2:n
        for e in es
            i::Int = src(e)
            j::Int = dst(e)
            generate_time += @elapsed begin
                f[m, i, j] = @variable(model, base_name="f[$m,$i,$j]", lower_bound=0, upper_bound=1)
                next!(progress)
                f[m, j, i] = @variable(model, base_name="f[$m,$j,$i]", lower_bound=0, upper_bound=1)
                next!(progress)
            end
            if generate_timeout_sec != Inf64 && generate_time > generate_timeout_sec
                cancel(progress)
                return
            end
        end
    end

    alpha::Int = 1

    if generate_timeout_sec != Inf64 && generate_time > generate_timeout_sec
        cancel(progress)
        return
    end

    for i in 2:n
        if has_edge(graph, i, alpha)
            for m in 2:n
                generate_time += @elapsed begin
                    fix(f[m, i, alpha], 0, force=true)
                end
                next!(progress)
                generate_time += @elapsed begin
                    fix(variable_by_name(model, "y[$i,$alpha]"), 0, force=true)
                end
                next!(progress)
                if generate_timeout_sec != Inf64 && generate_time > generate_timeout_sec
                    cancel(progress)
                    return
                end
            end
        end
    end



    #@constraint(model, source_flow, sum(sum(f[m, alpha, j] for j in 2:n if has_edge(graph, alpha, j)) for m in 2:n) == k)
    generate_time += @elapsed @constraint(model, one_real_root, sum(variable_by_name(model, "y[$alpha,$j]") for j in 2:n if has_edge(graph, alpha, j)) == 1)
    next!(progress)
    for i in 1:n
        for j in 1:n
            if has_edge(graph, i, j)
                for m in 2:n
                    generate_time += @elapsed begin
                        @constraint(model, f[m, i, j] <= variable_by_name(model, "y[$i,$j]"))
                    end
                    next!(progress)
                    if generate_timeout_sec != Inf64 && generate_time > generate_timeout_sec
                        cancel(progress)
                        return
                    end
                end
            end
        end
    end
    for m in 2:n
        generate_time += @elapsed begin
            @constraint(model, sum(f[m, i, m] for i in 1:n if has_edge(graph, i, m))
                            - sum(f[m, m, j] for j in 2:n if has_edge(graph, m, j)) == z[m])
        end
        next!(progress)
        if generate_timeout_sec != Inf64 && generate_time > generate_timeout_sec
            cancel(progress)
            return
        end
        generate_time += @elapsed begin
            @constraint(model, sum(f[m, alpha, j] for j in 2:n if has_edge(graph, alpha, j)) == z[m])
        end
        next!(progress)
        if generate_timeout_sec != Inf64 && generate_time > generate_timeout_sec
            cancel(progress)
            return
        end
        for v in 2:n
            if m != v
                generate_time += @elapsed begin 
                    @constraint(model, sum(f[m, i, v] for i in 1:n if has_edge(graph, i, v))
                                    - sum(f[m, v, j] for j in 2:n if has_edge(graph, v, j)) == 0)
                end
                next!(progress)
            end
            if generate_timeout_sec != Inf64 && generate_time > generate_timeout_sec
                cancel(progress)
                return
            end
        end
    end

    @constraint(model, treesize, sum(z[i] for i in 2:n) == k)
    next!(progress)
    finish!(progress)
end
