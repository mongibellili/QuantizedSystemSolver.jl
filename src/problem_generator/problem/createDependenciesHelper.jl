
"""
    EventDependencyStruct
 A struct that holds the event dependency information. It has the following fields:
  - `id::Int:` the id of the event
  - `evCont::Vector{Int}:` the index tracking used for HD & HZ. Also it is used to update q,quantum,recomputeNext when x is modified in an event
  - `evDisc::Vector{Int}:` the index tracking used for HD & HZ.
  - `evContRHS::Vector{Int}:` the index tracking used to update other Qs before executing the event

"""
struct EventDependencyStruct
  id::Int
  evCont::Vector{Int} #LHS: index tracking used for HD & HZ. Also it is used to update q,quantum,recomputeNext when x is modified in an event
  evDisc::Vector{Int} #LHS: index tracking used for HD & HZ.
  evContRHS::Vector{Int} #Here we look at the RHS: index tracking used to update other Qs before executing the event
end

function changeEventExprToFirstValue(ex::Expr)
  newEx=postwalk(ex) do a  
      if (a==:t)|| (a isa Expr && a.head == :ref && a.args[1]==:q) 
          outerRef=Expr(:ref)
          push!(outerRef.args,a)
          push!(outerRef.args,:(0))
          a=outerRef
      end
      return a
  end
  newEx
end





function find_assignments(expr::Expr)
    assignments = Expr[]
    if expr.head in [:(=), :+=, :-=, :*=, :/=]
        push!(assignments, expr)
    end
    for arg in expr.args
        if arg isa Expr
            append!(assignments, find_assignments(arg))
        end
    end
    return assignments
end




"""
    changeExprToFirstValue(ex::Expr)

 changes an expression in the form u[1] to an expression in the form u[1][0] inside exact jacobian expressions and inside events, because linear coefficients (a_{ii}) do not have derivatives, and updates in events affect the value of a variable directly and there is no need to update its higher derivatives. It is called by the [`restoreRef`](@ref) function for jacobian expressions, and called by the [`handleEvents`](@ref) function for events.
# Example:
```jldoctest
using QuantizedSystemSolver 
ex=:(q[i - 1]) 
newEx=QuantizedSystemSolver.changeExprToFirstValue(ex)

# output

:((q[i - 1])[0])
```
"""

function changeExprToFirstValue(expr)
    if expr isa Expr
        # Skip already transformed expressions like q[i][0] or q[i][1]
        if expr.head == :ref && expr.args[1] isa Expr && expr.args[1].head == :ref && expr.args[1].args[1] == :q
            return expr  # Already transformed → leave as is
        end

        # Match q[i] exactly
        if expr.head == :ref && expr.args[1] == :q && length(expr.args) == 2
            return Expr(:ref, expr, 0)  # q[i] → q[i][0]
        end

        # Recurse into subexpressions
        new_args = map(changeExprToFirstValue2, expr.args)
        return Expr(expr.head, new_args...)
    else
        return expr
    end
end



# Handles individual index terms
function _idxToString(ref)
    if ref isa Symbol || ref isa Int || ref isa String
        return string(ref)
    elseif ref isa Expr && ref.head == :call && length(ref.args) == 3
        op = ref.args[1]
        lhs = _idxToString(ref.args[2])
        rhs = _idxToString(ref.args[3])
        opname = _operatorName(op)
        return lhs * opname * rhs
    else
        error("Unsupported index expression: $ref")
    end
end

# Maps operators to readable names
function _operatorName(op::Symbol)
    return Dict(
        :+ => "plus",
        :- => "minus",
        :* => "times",
        :/ => "div",
        :^ => "pow"
    )[op]  # fallback is error if op not in Dict
end

function expr_to_flat_name(e)
    if e isa Symbol
        return string(e)

    elseif e isa Expr
        if e.head == :call
            fname = expr_to_flat_name(e.args[1])
            args = join(map(expr_to_flat_name, e.args[2:end]), "_")
            return "$(fname)_$(args)"

        elseif e.head == :vect
            items = join(map(expr_to_flat_name, e.args), "_")
            return "vect_$(items)"

        elseif e.head == :ref
            base = expr_to_flat_name(e.args[1])
            indices = join(map(expr_to_flat_name, e.args[2:end]), "_")
            return "$(base)_$(indices)"

        else
            return replace(string(e), r"[^\w]" => "_")
        end

    else
        return string(e)
    end
end




"""
    restoreRef(coefExpr,symDict)

This function is the opposite of symbolFromRef. After using the symbols in symbolic differentiation, it gets back expressions like p[i+Number] and q[i+Number][0] from symbols diplusNumber and qiplusNumber. Adding a zero to q variables is beacause q is a taylor variable while p is a vector.\n 
  # arguments:
- `coefExpr::Expr`: the expression to be changed
- `symDict::Dict{Symbol,Expr}`: the dictionary to store the translation of symbols of continous and discrete variables (q[i] <-> qi)

# Example:
```jldoctest
using QuantizedSystemSolver

symDict= Dict{Symbol, Expr}(:qi => :(q[i]), :q10 => :(q[10]), :q2 => :(q[2]), :qiminus1 => :(q[i - 1]), :q1 => :(q[1]))
coefExpr=:(1.5qiminus1) 
  newEx=QuantizedSystemSolver.restoreRef(coefExpr, symDict)


# output

:(1.5 * (q[i - 1])[0]) 
```
"""
function restoreRef(coefExpr,symDict)
  newEx=prewalk(coefExpr) do element# 
    if element isa Symbol && !(element in (:+,:-,:*,:/)) && haskey(symDict, element) 
        if String(element)[1] == 'q' 
          element=symDict[element]
          element=changeExprToFirstValue(element)# change q[1] to q[1][0]
        else  #element== :p or any other constant parameter
            element=symDict[element]  
        end
    end
    return element
  end#end postwalk
  newEx
end


function is_custom_function(expr::Expr)
    return expr.args[1] isa Symbol &&
           !(expr.args[1] in (:+, :-, :*, :/, :^, :%, :&, :|, :!, :(=), :(==), :!=, :<, :>, :<=, :>=)) &&
           !(isdefined(Base, expr.args[1]) && getfield(Base, expr.args[1]) isa Function)
end



