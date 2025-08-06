

"""
    transformFSimplecase(ex::Union{Float64,Int64,Expr,Symbol})

 transforms expressions of the right hand side of differential equations and zero-crossing functions to personalized ones that use caching to a form that can be used in the [`createT`](@ref) function. The right hand side of the equations should be a number or a q[i] term.
# Example:   
```jldoctest
using QuantizedSystemSolver
ex=:(q[2])
newEx=QuantizedSystemSolver.transformFSimplecase(ex);


# output

:(createT(q[2], cache[1])) 

```
"""  
function transformFSimplecase(ex::Union{Float64,Int64,Expr,Symbol})
    #  it s easier and healthier to leave the big case alone in one prewalk (below) when returning the size of the distributed cache
    ex=Expr(:call, :createT,ex)# case where rhs of eq is a number needed to change to taylor (in cache) to be used for qss ||  #case where rhs is q[i]...return cache that contains it
    cachexpr = Expr(:ref, :cache)   # add cache[1] to createT(ex,....)
    push!(cachexpr.args,1)
    push!( ex.args, cachexpr)    
    return ex
end 

is_binary(expr, op) = expr isa Expr &&
                      expr.head == :call &&
                      expr.args[1] == op &&
                      length(expr.args) == 3 &&
                      !(expr.args[2] isa Number && expr.args[3] isa Number) &&  # Prevent number+number
                      !(expr.args[2] isa Int || expr.args[3] isa Int)  # Avoid symbol+int transformations


is_unary(expr, op)  = expr isa Expr && expr.head == :call && expr.args[1] == op && length(expr.args) == 2

next_cache!(tracker::Ref{Int}) = (tracker[] += 1; Expr(:ref, :cache, tracker[]))

# --- Main transform ---
"""
    transformF(ex::Expr)

 transforms expressions of the right hand side of differential equations and zero-crossing functions to personalized ones that use caching to a form that can be used in functions like  [`addT`](@ref), [`subT`](@ref), [`mulT`](@ref), [`muladdT`](@ref). The right hand side of the equations can be any form of expression.
# Example:   
```jldoctest
using QuantizedSystemSolver
ex=:(q[2] - 2.0*q[1]*q[2],1)
newEx=QuantizedSystemSolver.transformF(ex);


# output

:((subT(q[2], mulTT(2.0, q[1], q[2], cache[2], cache[3]), cache[1]), 3))

```
"""
 function transformF(ex::Expr)
    tracker = Ref(0)  # For cache index tracking

    prewalk(ex) do x
        if x isa Expr && x.head == :call && x.args[1] isa Symbol # do not transform herlper functions (not base functions)
            sym=x.args[1]
            if !(sym in (:+, :-, :*, :/, :^, :%, :&, :|, :!, :(=),:(==), :!=, :<, :>, :<=, :>=)) && !(isdefined(Base, sym) && getfield(Base, sym) isa Function)
               # @show sym
                return x
            end
        end
        if is_binary(x, :-)
            a, b = x.args[2], x.args[3]
            if is_binary(a, :-)
                x.args = [:subsub, a.args[2], a.args[3], b, next_cache!(tracker)]
            elseif is_binary(a, :+)
                x.args = [:addsub, a.args[2], a.args[3], b, next_cache!(tracker)]
            elseif is_binary(a, :*)
                x.args = [:mulsub, a.args[2], a.args[3], b, next_cache!(tracker)]
            else
                x.args = [:subT, a, b, next_cache!(tracker)]
            end

        elseif is_unary(x, :-)
            x.args = [:negateT, x.args[2], next_cache!(tracker)]

        elseif is_binary(x, :+)
            a, b = x.args[2], x.args[3]
            if is_binary(a, :-)
                x.args = [:subadd, a.args[2], a.args[3], b, next_cache!(tracker)]
            elseif is_binary(a, :*)
                x.args = [:muladdT, a.args[2], a.args[3], b, next_cache!(tracker)]
            else
                x.args = [:addT, a, b, next_cache!(tracker)]
            end

        elseif x isa Expr && x.head == :call && x.args[1] == :+ && 4 ≤ length(x.args) ≤ 9
            x.args[1] = :addT
            push!(x.args, next_cache!(tracker))

        elseif is_binary(x, :*)
            x.args[1] = :mulT
            push!(x.args, next_cache!(tracker))

        elseif x isa Expr && x.head == :call && x.args[1] == :* && 4 ≤ length(x.args) ≤ 7
            x.args[1] = :mulTT
            push!(x.args, next_cache!(tracker))
            push!(x.args, next_cache!(tracker))

        elseif is_binary(x, :/)
            x.args[1] = :divT
            push!(x.args, next_cache!(tracker))

        elseif x isa Expr && x.head == :call && x.args[1] in (:exp, :log, :sqrt, :abs, :atan2)
            push!(x.args, next_cache!(tracker))

        elseif x isa Expr && x.head == :call && x.args[1] == :^
            x.args[1] = :powerT
            push!(x.args, next_cache!(tracker))

        elseif x isa Expr && x.head == :call && x.args[1] in (:cos, :sin, :tan, :atan)
            push!(x.args, next_cache!(tracker))
            push!(x.args, next_cache!(tracker))

        elseif x isa Expr && x.head == :call && x.args[1] in (:acos, :asin)
            push!(x.args, next_cache!(tracker))
            push!(x.args, next_cache!(tracker))
            push!(x.args, next_cache!(tracker))
        end

        return x
    end

    ex.args[2] = tracker[]  # total cache count
    return ex
end

#= function transform_expr!(x, tracker::Base.RefValue{Int})
    if is_binary(x, :-)
        a, b = x.args[2], x.args[3]
        if is_binary(a, :-)
            x.args = [:subsub, a.args[2], a.args[3], b, next_cache!(tracker)]
        elseif is_binary(a, :+)
            x.args = [:addsub, a.args[2], a.args[3], b, next_cache!(tracker)]
        elseif is_binary(a, :*)
            x.args = [:mulsub, a.args[2], a.args[3], b, next_cache!(tracker)]
        else
            x.args = [:subT, a, b, next_cache!(tracker)]
        end

    elseif is_unary(x, :-)
        x.args = [:negateT, x.args[2], next_cache!(tracker)]

    elseif is_binary(x, :+)
        a, b = x.args[2], x.args[3]
        if is_binary(a, :-)
            x.args = [:subadd, a.args[2], a.args[3], b, next_cache!(tracker)]
        elseif is_binary(a, :*)
            x.args = [:muladdT, a.args[2], a.args[3], b, next_cache!(tracker)]
        else
            x.args = [:addT, a, b, next_cache!(tracker)]
        end

    elseif x isa Expr && x.head == :call && x.args[1] == :+ && 4 ≤ length(x.args) ≤ 9
        x.args[1] = :addT
        push!(x.args, next_cache!(tracker))

    elseif is_binary(x, :*)
        x.args[1] = :mulT
        push!(x.args, next_cache!(tracker))

    elseif x isa Expr && x.head == :call && x.args[1] == :* && 4 ≤ length(x.args) ≤ 7
        x.args[1] = :mulTT
        push!(x.args, next_cache!(tracker))
        push!(x.args, next_cache!(tracker))

    elseif is_binary(x, :/)
        x.args[1] = :divT
        push!(x.args, next_cache!(tracker))

    elseif x isa Expr && x.head == :call && x.args[1] in (:exp, :log, :sqrt, :abs, :atan2)
        push!(x.args, next_cache!(tracker))

    elseif x isa Expr && x.head == :call && x.args[1] == :^
        x.args[1] = :powerT
        push!(x.args, next_cache!(tracker))

    elseif x isa Expr && x.head == :call && x.args[1] in (:cos, :sin, :tan, :atan)
        push!(x.args, next_cache!(tracker))
        push!(x.args, next_cache!(tracker))

    elseif x isa Expr && x.head == :call && x.args[1] in (:acos, :asin)
        push!(x.args, next_cache!(tracker))
        push!(x.args, next_cache!(tracker))
        push!(x.args, next_cache!(tracker))
    end

    return x
end

function transformF(ex::Expr)
    tracker = Ref(0)

    # Walk only the first argument of the expression
    ex.args[1] = walk_expr(ex.args[1], false, tracker)

    # Store final cache count in the second argument
    ex.args[2] = tracker[]

    return ex
end


function walk_expr(x, inside_helper::Bool,tracker::Base.RefValue{Int})
    if x isa Expr && x.head == :call
        fname = x.args[1]

        is_helper = fname isa Symbol &&
            !(fname in (:+, :-, :*, :/, :^, :%, :&, :|, :!, :(=), :(==), :!=, :<, :>, :<=, :>=)) &&
            !(isdefined(Base, fname) && getfield(Base, fname) isa Function)

        # Recurse with helper context if needed
        new_args = map(arg -> walk_expr(arg, inside_helper || is_helper, tracker), x.args)
        x = Expr(x.head, new_args...)
    elseif x isa Expr
        x = Expr(x.head, map(arg -> walk_expr(arg, inside_helper, tracker), x.args)...)
    end

    # Apply transformation only if NOT inside a helper call
    return (!inside_helper && x isa Expr) ? transform_expr!(x, tracker) : x
end
 =#