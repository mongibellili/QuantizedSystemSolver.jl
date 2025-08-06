


#A helper struct that holds the final IR. in addition, it used to return the number of `if_statements` and the number of helper functions written by the user. These two numbers are extracted during the normalization
struct probInfo
    ir::ODEFunctionIR
    numZC::Int
    helperFunSymSet::Int64
end


function convert_ints_except_indices(ex,  ::Val{false}) #convert all integers to float64 except indices
  if ex isa Expr
      if ex.head == :ref
          # Process the :ref expression, marking that we are inside a reference
          return Expr(:ref, ex.args[1], map(arg -> convert_ints_except_indices(arg, Val(true)), ex.args[2:end])...)
      else
          # Recursively process other expressions, but do NOT mark as inside :ref
          return Expr(ex.head, map(arg -> convert_ints_except_indices(arg, Val(false)), ex.args)...)
      end
  elseif ex isa Int
      # Convert to Float64 only if we are NOT inside a reference
      return  Float64(ex)
  else
      return ex
  end
end
function convert_ints_except_indices(ex,  ::Val{true})   #change any float put by user inside ref to int
  newEx=postwalk(ex) do element
      if element isa Float64 
        return Int(element)
      end
      return element
    end#end postwalk
  return newEx
end
#= function changeExprToFirstValue2(ex::Expr)
  newEx=postwalk(ex) do a  
      if a isa Expr && a.head == :ref && a.args[1]==:q  
          outerRef=Expr(:ref)
          push!(outerRef.args,a)
          push!(outerRef.args,:(0))
          a=outerRef
      end
      return a
  end
  newEx
end =#
function changeExprToFirstValue2(expr)
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


"""
    changeVarNames_params(ex::Expr,stateVarName::Symbol,discrParamName::Symbol,muteVar::Symbol,param::Dict{Symbol,Union{Float64,Int64,Expr,Symbol}},helperFunSymSet::Set{Symbol})

As the name suggests, this changes the continuous variables names to :q and the discrete variable name to :p and any mute variables to :i. It also plugs the parameters values from a parameter dictionary into the differential equations. The function changeVarNames_params has three methods. One for RHS of equations, one for if-statements when RHS is an expression, and one for if-statements when RHS is a symbol. This is method one. It has an additional symDict::Dict{Symbol,Expr} to collect the translation of symbols of continous and discrete variables (q[i] <-> qi). 

# arguments:
- `ex::Expr`: the expression to be changed
- `stateVarName::Symbol`: the name of the state variable
- `muteVar::Symbol`: the name of the mute variable
- `param::Dict{Symbol,Union{Float64,Int64,Expr,Symbol}}`: the dictionary of parameters


# Example:
```jldoctest
using QuantizedSystemSolver

(ex, stateVarName, discrParamName,muteVar, param) = (:(du[k] = u[k] * u[k - 1] * coef2), :u,:p, :k, Dict{Symbol, Union{Float64, Int64,Expr,Symbol}}(:coef1 => 2.0, :coef2 => 1.5))

  newEx=QuantizedSystemSolver.changeVarNames_params(ex, stateVarName,discrParamName, muteVar, param,Set([:f]))
(newEx, stateVarName, muteVar, param)
# output

(:(du[i] = q[i] * q[i - 1] * 1.5), :u, :k, Dict{Symbol, Union{Float64, Int64, Expr, Symbol}}(:coef1 => 2.0, :coef2 => 1.5))
```
"""
function changeVarNames_params(ex::Expr,stateVarName::Symbol,discrParamName::Symbol,muteVar::Symbol,param::Dict{Symbol,Union{Float64,Int64,Expr,Symbol}},helperFunSymSet::Set{Symbol})
  newEx=postwalk(ex) do element#postwalk to change var names and parameters
      if element isa Symbol   
        element=changeVarNames_params(element,stateVarName,discrParamName,muteVar,param,helperFunSymSet)
      elseif element isa Expr && element.head == :call && element.args[1] isa Symbol 
        sym=element.args[1]
        if !(sym in (:+, :-, :*, :/, :^, :%, :&, :|, :!, :(=),:(==), :!=, :<, :>, :<=, :>=)) && !(isdefined(Base, sym) && getfield(Base, sym) isa Function)
          push!(helperFunSymSet, sym)# collect the helper functions used in the rhs of the equations 
          element=changeExprToFirstValue2(element) # change q[1] to q[1][0]
        end
      elseif element isa Expr && element.head == :ref
        # Process the :ref expression, marking that we are inside a reference
        #return Expr(:ref, element.args[1], map(arg -> changeVarNames_params(arg, stateVarName, discrParamName, muteVar, param, helperFunSymSet), element.args[2:end])...)
      end
      return element
    end#end postwalk
    newEx = convert_ints_except_indices(newEx,Val(false))#convert all integers to float64 except indices
  newEx
end


"""
    changeVarNames_params(element::Symbol,stateVarName::Symbol,discrParamName::Symbol,muteVar::Symbol,param::Dict{Symbol,Union{Float64,Int64,Expr,Symbol}},helperFunSymSet::Set{Symbol})

This is method three of the function changeVarNames_params. It is for if-statements when RHS is a symbol. 
Again, it changes the symbol to :q if it is a continuous variable, to :p if it is a discrete variable, to :i if it is a mute variable, and to its corresponding value if it is a parameter.


"""
function changeVarNames_params(element::Symbol,stateVarName::Symbol,discrParamName::Symbol,muteVar::Symbol,param::Dict{Symbol,Union{Float64,Int64,Expr,Symbol}},helperFunSymSet::Set{Symbol})
          #@show stateVarName, element
        if haskey(param, element)#symbol is a parameter
            if param[element] isa Symbol
              element=param[element]#copy(::Symbol) does not exist
            else
              element=copy(param[element])
            end
          elseif element==stateVarName #symbol is a var
              element=:q 
          elseif element==:discrete || element==discrParamName#symbol is a discr var
              element=:p
          elseif element==muteVar #symbol is a mute var
              element=:i
          end
      return element
end
function changeVarNames_params(element::Number,stateVarName::Symbol,discrParamName::Symbol,muteVar::Symbol,param::Dict{Symbol,Union{Float64,Int64,Expr,Symbol}},helperFunSymSet::Set{Symbol})
   return element
end





"""
    recurse(e::Expr,flattened::Vector{Expr})

Break down a compound condition into basic components.
This function recursively decomposes a condition expression into its basic components, flattening nested logical operators (&&, ||) into a list of expressions. It returns a vector of flattened expressions.
# Example:
```jldoctest
using QuantizedSystemSolver
ex=:(u < 1 || u > 10)
flattened = Expr[]
QuantizedSystemSolver.recurse(ex, flattened)
flattened

# output

2-element Vector{Expr}:
 :(u < 1)
 :(u > 10)
```
"""
function recurse(e::Expr,flattened::Vector{Expr})
    if e isa Expr && e.head in [:&&, :||]
        for arg in e.args
            recurse(arg,flattened)
        end
    else
        push!(flattened, e)
    end
end
"""
    decompose_condition(cond::Expr)

Break down a compound condition into basic components.
Returns a tuple (kind, normalized conditions).
where kind is :or or :and depending on the original condition.
# Example:
```jldoctest
using QuantizedSystemSolver
ex=:(u < 1 || u > 10)
(kind,flattened)=QuantizedSystemSolver.decompose_condition(ex)
(kind,flattened)

# output

(:or, Expr[:(u < 1), :(u > 10)])
```
"""
function decompose_condition(cond::Expr)
    if cond.head in [:&&, :||]
        kind = cond.head == :&& ? :and : :or
        flattened = Expr[]
        recurse(cond,flattened)
        if length(flattened) > 2
            return (:multi, flattened)
        else
            return (kind, flattened)
        end
    else
        return (:single, [cond])
    end
end
"""
    to_zcf(expr::Expr)

Convert a condition expression to zero-crossing form.
This function transforms a condition expression into a zero-crossing form (ZCF) by rearranging the terms. For example, it converts expressions like `A < B` to `B - A ` and `A > B` to `A - B`.
It is used to prepare conditions for further processing in the normalization of IRs.    
# Example:
```jldoctest
using QuantizedSystemSolver
ex1 = :(u > 10)
ex2 = :(u < 1)
zcf1 = QuantizedSystemSolver.to_zcf(ex1)
zcf2 = QuantizedSystemSolver.to_zcf(ex2)
(zcf1, zcf2)
# output

(:(u - 10), :(1 - u))
"""
function to_zcf(expr)
    if expr.head == :call && expr.args[1] in [:>, :>=,:(==)]
        return Expr(:call, :-, expr.args[2], expr.args[3])
    elseif expr.head == :call && expr.args[1] in [:<, :<=]
        return Expr(:call, :-, expr.args[3], expr.args[2])
    end

end

"""
    process_if_condition(cond, stateVarName, discrParamName, param, helperFunSymSet)

Processes an `if` condition within the intermediate representation (IR).

# Arguments
- `cond`: The condition expression to be processed.
- `stateVarName`: The name of the state variable involved in the condition.
- `discrParamName`: The name of the discrete parameter relevant to the condition.
- `param`: Additional parameters required for processing.
- `helperFunSymSet`: A set of helper function symbols used during processing.

# Returns
Returns the processed representation of the `if` condition, potentially transformed for normalization within the IR.

# Notes
This function is intended for internal use in the normalization of IR.
"""
function process_if_condition(cond, stateVarName, discrParamName, param, helperFunSymSet)
    kind, original_conds = decompose_condition(cond)

    new_subconds = [
        changeVarNames_params(
            c,
            stateVarName,
            discrParamName,
            :nothing,
            param,
            helperFunSymSet,
        ) for c in original_conds
    ]
    cond_zcfs = [
            to_zcf(c)
            for c in new_subconds
    ]

    return kind, cond_zcfs,new_subconds
end

"""
    process_if_block(block_expr::Expr, stateVarName, discrParamName, param, helperFunSymSet)

Processes an `if` block (body) expression within the intermediate representation (IR).

# Arguments
- `block_expr::Expr`: The Julia expression representing the `if` block to be processed.
- `stateVarName`: The name of the state variable involved in the block.
- `discrParamName`: The name  of the discrete variable.
- `param`: Additional parameter(s) required for processing the block.
- `helperFunSymSet`: A set of helper function symbols used during processing.

# Returns
- The processed `if` block body (change of variable names).

# Notes
- This function is intended for internal use within the IR normalization pipeline.
"""
function process_if_block(block_expr::Expr, stateVarName, discrParamName, param, helperFunSymSet)
    # Ensure it's a block of statements
    stmts = block_expr.head == :block ? block_expr.args : [block_expr]
    to_delete = Int[]
    for (i, stmt) in enumerate(stmts)
        if stmt isa Expr && stmt.head in [:(=), :+=, :-=, :*=, :/=]
            lhs, rhs = stmt.args[1], stmt.args[2]
            stmt.args[1] = changeVarNames_params(lhs, stateVarName, discrParamName, :nothing, param, helperFunSymSet)
            #@show stmt.args[1]
            if stmt.args[1] isa Symbol && rhs isa Expr && rhs.head in [:call, :ref]
                stmt.args[2] = changeVarNames_params(rhs, stateVarName, discrParamName, :nothing, param, helperFunSymSet)
                param[lhs] = stmt.args[2]
                push!(to_delete, i)
            else
                stmt.args[2] = changeVarNames_params(rhs, stateVarName, discrParamName, :nothing, param, helperFunSymSet)
            end
        else
            stmts[i] = changeVarNames_params(stmt, stateVarName, discrParamName, :nothing, param, helperFunSymSet)
        end
    end

    # Remove unnecessary assignments
    for i in reverse(to_delete)
        splice!(stmts, i)
    end
    return length(stmts) == 1 ? stmts[1] : Expr(:block, stmts...)
end




"""
    process_if_expr(statement,stateVarName,discrParamName,param,helperFunSymSet)

Processes an `IfStatement` within the intermediate representation (IR) of a simple model. 
This function normalizes the given `IfStatement` according to the provided state variable name, discrete variable name, additional parameters, and a set of helper function symbols:
    - Applies variable substitution and body block processing.
    - Decomposes composite conditions (&&, ||).

# Arguments
- `statement::IfStatement`: The `IfStatement` node to be processed.
- `stateVarName`: The name of the state variable relevant to the normalization.
- `discrParamName`: The name of the discretization parameter.
- `param`: Additional parameter(s) required for normalization.
- `helperFunSymSet`: A set of symbols representing helper functions used during normalization.

# Returns
- `Vector{IfStatement}`: A vector of normalized `IfStatement` objects resulting from the processing.

# Notes
This function is typically used as part of the IR normalization pipeline in the `SimpleModelIR` module.
"""
function process_if_expr(statement,stateVarName,discrParamName,param,helperFunSymSet)
    cond = statement.condition
    if !(cond isa Expr && (
            (cond.head == :call && cond.args[1] in [:>, :>=, :<, :<=,:(==)]) ||
            cond.head in [:&&, :||]
        ))
        #@show statement
        return [statement]  # Already normalized
    end
    # === Normalize inequalities: A < B -> B - A > 0 ===
    if cond.head == :call && cond.args[1] in [:<, :<=]
        cond = Expr(:call, :>, cond.args[3], cond.args[2])
    elseif cond.head == :call && cond.args[1] == :(==)
        if_expr = statement.body
        if length(if_expr.args) == 3
            error("Equality condition must not have an else branch. Please contact the developers if this feature is needed.")
        elseif length(if_expr.args) == 2
            push!(if_expr.args, if_expr.args[2] ) # add the then block as else block, because rising or falling should trigger the same action.
        end
    elseif cond.head == :call && cond.args[1] in [:>, :>=]
        # already fine
    elseif cond.head in [:&&, :||]
        # OK for now, we'll handle decomposition later
    else
        error("Unsupported if condition: $cond")
    end
    # === Apply changeVarNames_params to whole condition ===
    kind, cond_zcfs,new_conds = process_if_condition(cond, stateVarName, discrParamName, param, helperFunSymSet)

    # === Process IF body ===
    if_expr = statement.body
    processed_expr = process_if_block(if_expr.args[2], stateVarName, discrParamName, param, helperFunSymSet)
    then_block = processed_expr.head == :block ? processed_expr : Expr(:block, processed_expr)
    if_expr.args[2] = then_block
    has_else = length(if_expr.args) == 3
   # @show has_else
    if has_else
        processed_expr = process_if_block(if_expr.args[3], stateVarName, discrParamName, param, helperFunSymSet)
        else_block = processed_expr.head == :block ? processed_expr : Expr(:block, processed_expr)
        if_expr.args[3] = else_block
    end
    # === Generate new IfStatement IRs ===
    new_if_statements = IfStatement[]
    if kind == :single
        if_expr.args[1] = cond_zcfs[1]
        push!(new_if_statements, IfStatement(cond_zcfs[1], if_expr))

    elseif kind == :and && !has_else
        zcfcond1, zcfcond2 = cond_zcfs
        cond1, cond2 = new_conds
        inner = Expr(:if,deepcopy(cond2), then_block) # deepcopy because inside events we do not want transformation to conditions in taylorEquationConstruction to be reflected back (pass by reference)
        outer = Expr(:if, deepcopy(cond1), inner)
        push!(new_if_statements, IfStatement(zcfcond1, outer))
        #  reverse the order
        inner_rev = Expr(:if, deepcopy(cond1), then_block)
        outer_rev = Expr(:if, deepcopy(cond2), inner_rev)
        push!(new_if_statements, IfStatement(zcfcond2, outer_rev))

    elseif kind == :and && has_else
        zcfcond1, zcfcond2 = cond_zcfs
        cond1, cond2 = new_conds
        inner = Expr(:if, deepcopy(cond2), then_block)
        outer = Expr(:if, deepcopy(cond1), inner, else_block)
        push!(new_if_statements, IfStatement(zcfcond1, outer))
        #  reverse the order
        inner_rev = Expr(:if, deepcopy(cond1), then_block)
        outer_rev = Expr(:if, deepcopy(cond2), inner_rev, else_block)
        push!(new_if_statements, IfStatement(zcfcond2, outer_rev))

    elseif kind == :or && !has_else
        zcfcond1, zcfcond2 = cond_zcfs
        cond1, cond2 = new_conds
        push!(new_if_statements, IfStatement(zcfcond1, Expr(:if, deepcopy(cond2), then_block)))
        push!(new_if_statements, IfStatement(zcfcond2, Expr(:if, deepcopy(cond1), then_block)))
    elseif kind == :or && has_else
        cond1, cond2 = cond_zcfs
        zcfcond1, zcfcond2 = cond_zcfs
        cond1, cond2 = new_conds
        fallback1 = Expr(:if, Expr(:call, :!,  deepcopy(cond2)), else_block)
        push!(new_if_statements, IfStatement(zcfcond1, Expr(:if, deepcopy(cond1), then_block, fallback1)))

        fallback2 = Expr(:if, Expr(:call, :!, deepcopy(cond1)), else_block)
        push!(new_if_statements, IfStatement(zcfcond2, Expr(:if, deepcopy(cond2), then_block, fallback2)))
    elseif kind == :multi
        renamed_cond = changeVarNames_params(statement.condition, stateVarName, discrParamName, :nothing, param, helperFunSymSet)
        if_expr.args[1] = deepcopy(renamed_cond)
        for cond_i in cond_zcfs
            new_expr = Expr(:if, deepcopy(cond_i), if_expr, if_expr)
            push!(new_if_statements, IfStatement(cond_i, new_expr))
        end
    end

    return new_if_statements
end



"""
    normalize_ir(ir, stateVarName::Symbol, discrParamName::Symbol)

Normalizes the intermediate representation (IR) of an ODE function.

# Arguments
- `ir::ODEFunctionIR`: The intermediate representation of the ODE function to be normalized.
- `stateVarName::Symbol`: The symbol representing the state variable in the IR.
- `discrParamName::Symbol`: The symbol representing the name of the discrete variable.

# Returns
- A normalized version of the input `ODEFunctionIR`.

# Description
This function processes the given ODE function IR, normalizing its structure with respect to the specified state variable and discretization parameter.
The normalization includes renaming variables, swapping parameters, and decomposing composite if-statements.
"""
function normalize_ir(ir, stateVarName::Symbol, discrParamName::Symbol)
    param = Dict{Symbol, Union{Float64, Int64, Expr, Symbol}}()
    helperFunSymSet = Set{Symbol}()
    numZC = 0
    normalized_statements = []
    for statement in ir.statements
        if statement isa AssignStatement
            lhs, rhs = statement.lhs, statement.rhs
            if lhs isa Symbol && rhs isa Number
                param[lhs] = Float64(rhs)
            elseif rhs isa Symbol && lhs isa Expr && lhs.head == :tuple
                base_sym = rhs == discrParamName ? :p : rhs == stateVarName ? :q : rhs
                for (i, sym) in enumerate(lhs.args)
                    param[sym] = :($base_sym[$i])
                end
            elseif lhs isa Symbol && rhs isa Expr && (rhs.head == :call || rhs.head == :ref)
                new_rhs = changeVarNames_params(rhs, stateVarName, discrParamName, :nothing, param, helperFunSymSet)
                statement.rhs = new_rhs
                param[lhs] = new_rhs
            elseif lhs isa Symbol && rhs isa Expr && rhs.head in [:vect, :tuple]
                if rhs.head == :vect
                    @warn "Vector literal assigned to `$(lhs)` may cause allocations. Consider using a tuple or pass it through parameters if appropriate."
                end
                new_rhs = changeVarNames_params(rhs, stateVarName, discrParamName, :nothing, param, helperFunSymSet)
                statement.rhs = new_rhs
            elseif lhs isa Expr && lhs.head == :ref
                if lhs.args[2] isa Expr && lhs.args[2].head == :call
                    lhs.args[2] = Int(eval(changeVarNames_params(lhs.args[2], stateVarName, discrParamName, :nothing, param, helperFunSymSet)))
                end
                if rhs isa Expr && rhs.head != :vect
                    new_rhs = changeVarNames_params(rhs, stateVarName, discrParamName, :nothing, param, helperFunSymSet)
                    statement.rhs = new_rhs
                    if haskey(param, lhs.args[2])
                        lhs.args[2] = param[lhs.args[2]] 
                    end
                elseif rhs isa Symbol && haskey(param, rhs)
                    statement.rhs = param[rhs]
                end
            end
             push!(normalized_statements, statement)
        elseif statement isa ForStatement
            muteVar = statement.var
            b = statement.start
            niter = statement.stop
            for body_statement in statement.body              
                if body_statement isa AssignStatement
                    body_statement.rhs = changeVarNames_params(body_statement.rhs, stateVarName, discrParamName, muteVar, param, helperFunSymSet)
                end
            end
            if !(b isa Int64)
                statement.start = Int(eval(changeVarNames_params(b, stateVarName, discrParamName, :nothing, param, helperFunSymSet)))
            end
            if !(niter isa Int64)
                statement.stop = Int(eval(changeVarNames_params(niter, stateVarName, discrParamName, :nothing, param, helperFunSymSet)))
            end
             push!(normalized_statements, statement)
        elseif statement isa IfStatement
           
           # cond = statement.condition
           # @show "before normalize if", statement
            new_ifs = process_if_expr(statement, stateVarName, discrParamName, param, helperFunSymSet)
           # @show new_ifs
            numZC += length(new_ifs)
            append!(normalized_statements, new_ifs) 
        else
            if statement isa ExprStatement && statement.expr.head == :function
                # put function name in param so that call to it can be replaced
                param[statement.expr.args[1].args[1]]=:f_
            end
             
             append!(normalized_statements, [statement])
        end
    end
    ir.statements = normalized_statements
    numHelperF = length(helperFunSymSet)
    return probInfo(ir,numZC, numHelperF)
end


