




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
        new_args = map(changeExprToFirstValue, expr.args)
        return Expr(expr.head, new_args...)
    else
        return expr
    end
end


"""
    rename_and_swap(ex::Expr,stateVarName::Symbol,discrParamName::Symbol,muteVar::Symbol,stack::SymbolTableStack,helperFunSymSet::Set{Symbol})

As the name suggests, this changes the continuous variables names to :q and the discrete variable name to :p and any mute variables to :i. It also plugs the parameters values from a parameter dictionary into the differential equations. The function rename_and_swap has three methods. One for RHS of equations, one for if-statements when RHS is an expression, and one for if-statements when RHS is a symbol. This is method one. It has an additional symDict::Dict{Symbol,Expr} to collect the translation of symbols of continous and discrete variables (q[i] <-> qi). 

# arguments:
- `ex::Expr`: the expression to be changed
- `stateVarName::Symbol`: the name of the state variable
- `muteVar::Symbol`: the name of the mute variable
- `stack::SymbolTableStack`: contains a dictionary of parameters


# Example:
```jldoctest
using QuantizedSystemSolver

(ex, stateVarName, discrParamName,muteVar, stack) = (:(du[k] = u[k] * u[k - 1] * coef2), :u,:p, :k, Dict{Symbol, Union{Float64, Int64,Expr,Symbol}}(:coef1 => 2.0, :coef2 => 1.5))

  newEx=QuantizedSystemSolver.rename_and_swap(ex, stateVarName,discrParamName, muteVar, stack,Set([:f]))
(newEx, stateVarName, muteVar, stack)
# output

(:(du[i] = q[i] * q[i - 1] * 1.5), :u, :k, Dict{Symbol, Union{Float64, Int64, Expr, Symbol}}(:coef1 => 2.0, :coef2 => 1.5))
```
"""


function is_custom_function(expr::Expr;custom=true)
    # expr is assumed to be a :call
    expr.head == :call || error("Expected Expr with head :call, got $(expr.head). Please report this bug")
    f = expr.args[1]
    if f isa Symbol
        # builtin operator or Base function?
        if f in SKIP_SYMBOLS # (:+, :-, :*, :/, :^, :%, :&, :|, :!, :(=), :(==), :!=, :<, :>, :<=, :>=)
            return false
        elseif custom && isdefined(Base, f) && getfield(Base, f) isa Function
            return false
        else
            return true   # user-defined symbol
        end
    else
        return true       # anything else (lambda, getfield, (g()), etc.) → treat as custom
    end
end


function rename_and_swap(ex::Expr,stateVarName::Symbol,discrParamName::Symbol,muteVar::Symbol,stack::SymbolTableStack,helperFunSymSet::Set{Symbol})
  newEx=postwalk(ex) do element#postwalk to change var names and parameters
      if element isa Symbol   
        element=rename_and_swap(element,stateVarName,discrParamName,muteVar,stack,helperFunSymSet)
      elseif element isa Expr && element.head == :call #&& element.args[1] isa Symbol 
   #=         sym=element.args[1]
        if !(sym in SKIP_SYMBOLS) && !(isdefined(Base, sym) && getfield(Base, sym) isa Function)  =#
        if is_custom_function(element)
            if element.args[1] isa Symbol
                 push!(helperFunSymSet, element.args[1])# collect the helper functions used in the rhs of the equations 
            else
                 push!(helperFunSymSet, Symbol(element.args[1]) )
            end
          element=changeExprToFirstValue(element) # change q[1] to q[1][0]
        end
      elseif element isa Expr && element.head == :ref  #change [q[1], 6.6, q[1]][1] to q[1] 
            container = element.args[1]
            idx = element.args[2]

            if container isa Expr && container.head == :vect && idx isa Int
                return container.args[idx]
            end
        # Process the :ref expression, marking that we are inside a reference
        #return Expr(:ref, element.args[1], map(arg -> rename_and_swap(arg, stateVarName, discrParamName, muteVar, stack, helperFunSymSet), element.args[2:end])...)
      end
      return element
    end#end postwalk
    newEx = convert_ints_except_indices(newEx,Val(false))#convert all integers to float64 except indices
  newEx
end


"""
    rename_and_swap(element::Symbol,stateVarName::Symbol,discrParamName::Symbol,muteVar::Symbol,stack::SymbolTableStack,helperFunSymSet::Set{Symbol})

This is method three of the function rename_and_swap. It is for if-statements when RHS is a symbol. 
Again, it changes the symbol to :q if it is a continuous variable, to :p if it is a discrete variable, to :i if it is a mute variable, and to its corresponding value if it is a parameter.


"""
function rename_and_swap(element::Symbol,stateVarName::Symbol,discrParamName::Symbol,muteVar::Symbol,stack::SymbolTableStack,helperFunSymSet::Set{Symbol})
    
    if element in SKIP_SYMBOLS || isdefined(Base, element)#symbol is sign or in base like cos...
    elseif element==stateVarName #symbol is a var
        element=:q 
    elseif element==:discrete || element==discrParamName#symbol is a discr var: p can be used both for simple param passing and for discrete variables (the word discrete if in future we want to separate them)
        element=:p
    elseif element==muteVar #symbol is a mute var
        element=:i
    else
        entry = lookup(stack, element) # lookup the symbol in the stack
        if !isnothing(entry) && entry.safe_to_inline
            val = entry.value
            element = (val isa Expr) ?  copy(val) : val  # Copy if Expr to avoid IR side effects
        end
    end
    return element
end
function rename_and_swap(element::Union{Number,String,Char,Bool},stateVarName::Symbol,discrParamName::Symbol,muteVar::Symbol,stack::SymbolTableStack,helperFunSymSet::Set{Symbol})
   return element
end


function contains_vect_or_custom(ex; skip_self=false)
    found = false
    postwalk(ex) do node
        if node isa Expr
            if (node.head in [:vect, :vcat]) || (node.head == :call && is_custom_function(node))
                # if skip_self=true, only mark found if node is not the root
                if !(skip_self && node === ex)
                    found = true
                end
            end
        end
        return node
    end
    return found
end
# --- Decide AssignAction based on mode, user flags, RHS ---
function decide_action(rhs; mode::InlineMode, user_inline::Bool=false, user_noinline::Bool=false)
    if mode == MANUAL
        return user_inline ? RegisterSwapRemove() : RegisterNoSwapKeep()
    elseif mode == AUTO
        # Hybrid AUTO: scalar custom functions are inlined with warning
        if user_inline
            return RegisterSwapRemove()
        elseif user_noinline
            return RegisterNoSwapKeep()
        elseif rhs isa Expr
            if rhs.head == :call
                if is_custom_function(rhs)
                    return RegisterSwapRemoveWarn()  # inline with warning
                elseif !contains_vect_or_custom(rhs)
                    return RegisterSwapRemove()      # safe built-in op
                else
                    return RegisterNoSwapKeep()
                end
            elseif rhs.head in [:vect, :vcat]
                return RegisterNoSwapKeep()          # never inline vectors
            elseif rhs.head in [:tuple, :curly]
                return RegisterNoSwapKeep()          # never inline tuples/curly
            elseif rhs.head == :ref  && (rhs.args[1]==:q || rhs.args[1]==:p)
                 return RegisterSwapRemove() 
            else
                return RegisterNoSwapKeep()
            end
        else
            return RegisterSwapRemove()              # literal or symbol
        end
    else # FULL
        return user_noinline ? RegisterNoSwapKeep() : RegisterSwapRemove()
    end
end

# --- Main handler ---
function handleAssignStatement!(statement::AssignStatement,
                                stateVarName::Symbol, discrParamName::Symbol,
                                stack::SymbolTableStack, helperFunSymSet::Set{Symbol},
                                muteVar::Symbol,
                                mode::InlineMode;
                                user_inline::Bool=false,
                                user_noinline::Bool=false)

    currentTable = peek(stack)
    lhs, rhs = statement.lhs, statement.rhs

    # rename / swap RHS first
    rhs = rename_and_swap(rhs, stateVarName, discrParamName, muteVar, stack, helperFunSymSet)
    statement.rhs = rhs

    # Decide policy action
    action = decide_action(rhs; mode=mode, user_inline=user_inline, user_noinline=user_noinline)

    # Emit warning if Hybrid AUTO inlined custom function
    if action.warn
        @warn "Inlined custom function $(rhs). If it returns a vector/matrix, this may cause allocations."
    end

    # --- Handle LHS shape ---
    if lhs isa Symbol
        if action.register
            add_symbol!(currentTable, lhs, rhs, action.swap)
        end
        statement.keep_assignment = action.keep

    elseif lhs isa Expr
        if lhs.head == :tuple
            # Tuple destructuring: per element
            keep = action.keep
            base_sym = rhs == discrParamName ? :p : rhs == stateVarName ? :q : rhs
            for (i, sym) in enumerate(lhs.args)
                if sym isa Symbol
                    if action.register
                        add_symbol!(currentTable, sym, :($base_sym[$i]), action.swap)
                    end
                else
                    keep = true
                    if sym isa Expr && sym.head == :ref
                        sym.args[1] = rename_and_swap_lhs(sym.args[1], stateVarName, discrParamName, muteVar, stack, helperFunSymSet)
                    end
                end
            end
            statement.keep_assignment = keep

        elseif lhs.head == :ref
            # Cannot register ref LHS; just rename indices
            lhs.args[1] = rename_and_swap_lhs(lhs.args[1], stateVarName, discrParamName, muteVar, stack, helperFunSymSet)
            if lhs.args[2] isa Expr && lhs.args[2].head == :call
                lhs.args[2] = Int(eval(rename_and_swap(lhs.args[2], stateVarName, discrParamName, muteVar, stack, helperFunSymSet)))
            end
            statement.keep_assignment = true

        else
            # Fallback: e.g., obj.field = rhs
            statement.keep_assignment = true
        end
    end

    return nothing
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
function normalize_ir(ir_statements,stack, stateVarName::Symbol, discrParamName::Symbol,muteVar::Symbol,inline_mode::InlineMode)
    helperFunSymSet = Set{Symbol}()
    numZC = 0
    normalized_statements = Vector{AbstractODEStatement}()
    for statement in ir_statements
        if statement isa AssignStatement
            handleAssignStatement!(statement, stateVarName, discrParamName, stack, helperFunSymSet,muteVar,inline_mode) # notice that muteVar is passed to handleAssignStatement! so that it is not swapped
            push!(normalized_statements, statement)
        elseif statement isa ForStatement
            localSymTable = SymbolTable()
            push_scope!(stack, localSymTable)
            muteVar = statement.var
            b = statement.start
            niter = statement.stop
            if !(b isa Int64)
                statement.start = Int(eval(rename_and_swap(b, stateVarName, discrParamName, :nothing, stack, helperFunSymSet)))
            end
            if !(niter isa Int64)
                statement.stop = Int(eval(rename_and_swap(niter, stateVarName, discrParamName, :nothing, stack, helperFunSymSet)))
            end
            # recursive call to normalize the body of the for loop: no need to get the return because here all we want is the swapping: 
            # notice that inside normalize, only process_if_expr gets a new expr...so, an if statment under a for statment won't change (which what we want because it s not an event)
            normalize_ir(statement.body,stack, stateVarName, discrParamName,muteVar,inline_mode) 

            push!(normalized_statements, statement)
            pop_scope!(stack)
        elseif statement isa IfStatement
            localSymTable = SymbolTable()
            push_scope!(stack, localSymTable)
            new_ifs = process_if_expr(statement, stateVarName, discrParamName, stack, helperFunSymSet)
            numZC += length(new_ifs)
            append!(normalized_statements, new_ifs) 
            pop_scope!(stack)
        elseif statement isa ExprStatement && statement.expr.head == :function
                # put function name in stack so that call to it can be replaced. This is a closure function and not allowed inside RuntimeGeneratedFunctions, so i pass it to integrator as :f_
                currentTable = peek(stack)
                add_symbol!(currentTable, statement.expr.args[1].args[1], :f_)
                push!(normalized_statements, statement)
        elseif statement isa ExprStatement && statement.expr.head == :macrocall && statement.expr.args[1] == Symbol("@inline") && statement.expr.args[3] isa Expr && statement.expr.args[3].head == :(=)
            lhs = statement.expr.args[3].args[1]
            rhs = statement.expr.args[3].args[2]

            # create an AssignStatement from the @inline. no need to handle @no_inline because it goes under the default behavior: else below
            assignment_from_macro = AssignStatement(lhs, rhs)
            handleAssignStatement!(assignment_from_macro, stateVarName, discrParamName, stack, helperFunSymSet,muteVar,inline_mode,user_inline=true) # true means user_inline
          
            push!(normalized_statements, assignment_from_macro)
        elseif statement isa ExprStatement && statement.expr.head == :macrocall && statement.expr.args[1] == Symbol("@noinline") && statement.expr.args[3] isa Expr && statement.expr.args[3].head == :(=)
            lhs = statement.expr.args[3].args[1]
            rhs = statement.expr.args[3].args[2]

            # create an AssignStatement from the @inline. no need to handle @no_inline because it goes under the default behavior: else below
            assignment_from_macro = AssignStatement(lhs, rhs)
            handleAssignStatement!(assignment_from_macro, stateVarName, discrParamName, stack, helperFunSymSet,muteVar,inline_mode,user_noinline=true) # true means user_inline
          
            push!(normalized_statements, assignment_from_macro)
        else #ExprStatement
            new_expr = rename_and_swap(statement.expr, stateVarName, discrParamName, muteVar, stack, helperFunSymSet) # rhs = f(q[i], p[i]) or rhs=q[i]
            statement.expr = new_expr
            push!(normalized_statements, statement)
        end
    end

    #@show ir_statements
    numHelperF = length(helperFunSymSet)
    #return probInfo(numZC, numHelperF) # return the number of zero-crossings and the number of helper functions
    return normalized_statements,numZC, numHelperF
end


