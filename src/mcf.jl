
using JuMP
using LightGraphs, SimpleWeightedGraphs


function multi_commodity_flow!(model,
                                graph :: SimpleWeightedGraph,
                                k :: Int;
                                generate_timeout_sec::Float64 = Inf64,
                                isdebug::Bool=false,
                                total_memory::Union{Int,UInt,Float64}=Base.Sys.total_memory(),
                                memory_limit_ratio::Union{Rational{Int},Rational{UInt},Float64}=1//20)
    n::Int = nv(graph)
    nes::Int = ne(graph)
    es::Vector{SimpleWeightedEdge} = collect(edges(graph))
    sort!(es, by=dst)
    sort!(es, by=src, alg=MergeSort)

    progress::Progress = Progress( 
        (n-1) * 4               # z[2:n], Bin, >= 0, <= 1
        + 2 * 3 * (n-1) * nes   # f[m,i,j] f[m,j,i], <= 0, >= 1
        + (n-1)                 # fix(y[i,alpha])
        + (n-1) * (n-1)         # fix(f[m,i,alpha])
        + 1                     # one_real_root
        + nes * (n-1) * 2       # f[m,i,j] <= y[i,j], f[m,j,i] <= y[j,i]
        + n-1                   # sum(f[m,m,j])- sum(f[m,j,m]) == z[m]
        + n-1                   # sum(f[m,alpha,j]) == z[m]
        + (n-1) * (n-1) - n     # sum(f[m,i,v]) - sum(f[m,v,j]) == 0
        + 1                     # treesize
        ,
                                    dt=0.1,
                                    barglyphs=BarGlyphs('[', '#', ['-','~','+','*','=','>','#'], ' ', ']'),
                                    enabled=isdebug)
    generate_time::Float64 = 0

    generate_time += @elapsed begin
        @variable(model, z[2:n], Bin, lower_bound=0, upper_bound=1)
    end
    ProgressMeter.update!(progress, (n-1) * 4)
    f::Dict{Tuple{Int, Int, Int},VariableRef} = Dict()
    for m in 2:n
        for e in es
            i::Int = src(e)
            j::Int = dst(e)
            generate_time += @elapsed begin
                f[m, i, j] = @variable(model, base_name="f[$m,$i,$j]", lower_bound=0, upper_bound=1)
                ProgressMeter.update!(progress, 3)
                f[m, j, i] = @variable(model, base_name="f[$m,$j,$i]", lower_bound=0, upper_bound=1)
                ProgressMeter.update!(progress, 3)
            end
            if (generate_timeout_sec != Inf64 && generate_time > generate_timeout_sec) || (Base.Sys.free_memory()//total_memory <= memory_limit_ratio)
                cancel(progress)
                return
            end
        end
    end

    alpha::Int = 1

    if (generate_timeout_sec != Inf64 && generate_time > generate_timeout_sec) || (Base.Sys.free_memory()//total_memory <= memory_limit_ratio)
        cancel(progress)
        return
    end

    for i in 2:n
        if has_edge(graph, i, alpha)
            generate_time += @elapsed begin
                fix(variable_by_name(model, "y[$i,$alpha]"), 0, force=true)
            end
            next!(progress)
            for m in 2:n
                generate_time += @elapsed begin
                    fix(f[m, i, alpha], 0, force=true)
                end
                next!(progress)
                if (generate_timeout_sec != Inf64 && generate_time > generate_timeout_sec) || (Base.Sys.free_memory()//total_memory <= memory_limit_ratio)
                    cancel(progress)
                    return
                end
            end
        end
    end



    #@constraint(model, source_flow, sum(sum(f[m, alpha, j] for j in 2:n if has_edge(graph, alpha, j)) for m in 2:n) == k)
    generate_time += @elapsed @constraint(model, one_real_root, sum(variable_by_name(model, "y[$alpha,$j]") for j in 2:n if has_edge(graph, alpha, j)) == 1)
    next!(progress)
    for e in es
        i :: Int = src(e)
        j :: Int = dst(e)
        for m in 2:n
            generate_time += @elapsed begin
                @constraint(model, f[m, i, j] <= variable_by_name(model, "y[$i,$j]"))
            end
            next!(progress)
            generate_time += @elapsed begin
                @constraint(model, f[m, j, i] <= variable_by_name(model, "y[$j,$i]"))
            end
            next!(progress)
            if (generate_timeout_sec != Inf64 && generate_time > generate_timeout_sec) || (Base.Sys.free_memory()//total_memory <= memory_limit_ratio)
                cancel(progress)
                return
            end
        end
    end
    for m in 2:n
        generate_time += @elapsed begin
            @constraint(model, sum(f[m, i, m] for i in 1:n if has_edge(graph, i, m))
                            - sum(f[m, m, j] for j in 2:n if has_edge(graph, m, j)) == z[m])
        end
        next!(progress)
        if (generate_timeout_sec != Inf64 && generate_time > generate_timeout_sec) || (Base.Sys.free_memory()//total_memory <= memory_limit_ratio)
            cancel(progress)
            return
        end
        generate_time += @elapsed begin
            @constraint(model, sum(f[m, alpha, j] for j in 2:n if has_edge(graph, alpha, j)) == z[m])
        end
        next!(progress)
        if (generate_timeout_sec != Inf64 && generate_time > generate_timeout_sec) || (Base.Sys.free_memory()//total_memory <= memory_limit_ratio)
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
            if (generate_timeout_sec != Inf64 && generate_time > generate_timeout_sec) || (Base.Sys.free_memory()//total_memory <= memory_limit_ratio)
                cancel(progress)
                return
            end
        end
    end

    @constraint(model, treesize, sum(z[i] for i in 2:n) == k)
    next!(progress)
    finish!(progress)
end
