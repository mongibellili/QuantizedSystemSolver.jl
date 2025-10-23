
function changeT(expr)
    if expr == :t
        return :(t[0])
    elseif expr isa Expr
        return Expr(expr.head, map(changeT, expr.args)...)
    else
        return expr
    end
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
    if e isa Expr
        if e.head == :call
            fname = expr_to_flat_name(e.args[1])
            args  = join(map(expr_to_flat_name, e.args[2:end]), "_")
            return "$(fname)_$(args)"

        elseif e.head == :vect
            items = join(map(expr_to_flat_name, e.args), "_")
            return "vect_$(items)"

        elseif e.head == :tuple
            items = join(map(expr_to_flat_name, e.args), "_")
            return "tuple_$(items)"

        elseif e.head == :ref
            base    = expr_to_flat_name(e.args[1])
            indices = join(map(expr_to_flat_name, e.args[2:end]), "_")
            return "$(base)_$(indices)"

        else
            return replace(string(e), r"[^\w]" => "_")
        end

    else
        return replace(string(e), r"[^\w]" => "_")
        
    end
end

function restoreRef(coefExpr,symDict)
  newEx=prewalk(coefExpr) do element# 
    if element isa Symbol && !(element in (:+,:-,:*,:/)) && haskey(symDict, element) 
        if String(element)[1] == 'q' 
          element=symDict[element]
          #element=changeExprToFirstValue(element)# change q[1] to q[1][0]
        else  #element== :p or any other constant parameter
            element=symDict[element]  
        end
    end
    return element
  end#end 
  newEx
end


function is_custom_function(expr::Expr)
    # expr is assumed to be a :call
    f = expr.args[1]
    if f isa Symbol
        # builtin operator or Base function?
        if f in (:+, :-, :*, :/, :^, :%, :&, :|, :!, :(=), :(==), :!=, :<, :>, :<=, :>=)
            return false
        elseif isdefined(Base, f) && getfield(Base, f) isa Function
            return false
        else
            return true   # user-defined symbol
        end
    else
        return true       # anything else (lambda, getfield, (g()), etc.) → treat as custom
    end
end




@doc """
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
""" restoreRef