

function rename_and_swap_lhs(element::Symbol,stateVarName::Symbol,discrParamName::Symbol,muteVar::Symbol,stack::SymbolTableStack,helperFunSymSet::Set{Symbol})
    if element==stateVarName #symbol is a var
        element=:q 
    elseif element==:discrete || element==discrParamName#symbol is a discr var: p can be used both for simple param passing and for discrete variables (the word discrete if in future we want to separate them)
        element=:p
    elseif element==muteVar #symbol is a mute var
        element=:i
    else
        entry = lookup(stack, element) # lookup the symbol in the stack
        if !isnothing(entry) && entry.safe_to_inline
            val = entry.value
            if val isa Expr 
                if val.head == :ref # if symbol in lhs of event to be replaced, then allow only refs
                    element = copy(val)  # Copy if Expr to avoid IR side effects
                end
            else #
                element = val
            end
        end
    end
    return element
end


function rename_and_swap_lhs(ex::Expr,stateVarName::Symbol,discrParamName::Symbol,muteVar::Symbol,stack::SymbolTableStack,helperFunSymSet::Set{Symbol})
  newEx=postwalk(ex) do element#postwalk to change var names and parameters
      if element isa Symbol   
        element=rename_and_swap(element,stateVarName,discrParamName,muteVar,stack,helperFunSymSet)
      end
      return element
    end#end postwalk
    #newEx = convert_ints_except_indices(newEx,Val(false))#convert all integers to float64 except indices
  newEx
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
        if expr.args[3] == 0.0 # if expr.args[3]==0.0, no need for substraction
            return  expr.args[2]
        else
            return Expr(:call, :-, expr.args[2], expr.args[3])
        end

    elseif expr.head == :call && expr.args[1] in [:<, :<=]
        return Expr(:call, :-, expr.args[3], expr.args[2])
    end

end

"""
    process_if_condition(cond, stateVarName, discrParamName, stack, helperFunSymSet)

Processes an `if` condition within the intermediate representation (IR).

# Arguments
- `cond`: The condition expression to be processed.
- `stateVarName`: The name of the state variable involved in the condition.
- `discrParamName`: The name of the discrete parameter relevant to the condition.
- `stack`: Additional parameters required for processing.
- `helperFunSymSet`: A set of helper function symbols used during processing.

# Returns
Returns the processed representation of the `if` condition, potentially transformed for normalization within the IR.

# Notes
This function is intended for internal use in the normalization of IR.
"""
function process_if_condition(cond, stateVarName, discrParamName, stack, helperFunSymSet)
    kind, original_conds = decompose_condition(cond)

    new_subconds = [
        rename_and_swap(
            c,
            stateVarName,
            discrParamName,
            :nothing,
            stack,
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
    process_if_block(block_expr::Expr, stateVarName, discrParamName, stack, helperFunSymSet)

Processes an `if` block (body) expression within the intermediate representation (IR).

# Arguments
- `block_expr::Expr`: The Julia expression representing the `if` block to be processed.
- `stateVarName`: The name of the state variable involved in the block.
- `discrParamName`: The name  of the discrete variable.
- `stack`: Additional parameter(s) required for processing the block.
- `helperFunSymSet`: A set of helper function symbols used during processing.

# Returns
- The processed `if` block body (change of variable names).

# Notes
- This function is intended for internal use within the IR normalization pipeline.
"""
function process_if_block(block_expr::Expr, stateVarName, discrParamName, stack, helperFunSymSet)
    # Ensure it's a block of statements
    stmts = block_expr.head == :block ? block_expr.args : [block_expr]
    to_delete = Int[]
    for (i, stmt) in enumerate(stmts)
        if stmt isa Expr && stmt.head in [:(=), :+=, :-=, :*=, :/=]
            lhs, rhs = stmt.args[1], stmt.args[2]
            
            stmt.args[1] = rename_and_swap_lhs(lhs, stateVarName, discrParamName, :nothing, stack, helperFunSymSet)  # even if lhs is a symbol, it still goes through the swapping first
            #@show stmt.args[1]
            if stmt.args[1] isa Symbol && (!(rhs isa Expr) || rhs isa Expr && rhs.head in [:call, :ref] )   # we cannot use/write  "if lhs" because stmt.args[1] can be a ref while lhs still symbol and this assignment would be removed which is bad
                stmt.args[2] = rename_and_swap(rhs, stateVarName, discrParamName, :nothing, stack, helperFunSymSet)
                currentTable = peek(stack)
                add_symbol!(currentTable, lhs, stmt.args[2])# safe_to_inline already set to true by Default
                if stmt.args[1] !=:t push!(to_delete, i) end # do not delete time assignments
            else
                stmt.args[2] = rename_and_swap(rhs, stateVarName, discrParamName, :nothing, stack, helperFunSymSet)
            end
            
        else # If it's not an assignment
            stmts[i] = rename_and_swap(stmt, stateVarName, discrParamName, :nothing, stack, helperFunSymSet)
        end
    end

    # Remove unnecessary assignments
    for i in reverse(to_delete)
        splice!(stmts, i)
    end
    return length(stmts) == 1 ? stmts[1] : Expr(:block, stmts...)
end




"""
    process_if_expr(statement,stateVarName,discrParamName,stack,helperFunSymSet)

Processes an `IfStatement` within the intermediate representation (IR) of a simple model. 
This function normalizes the given `IfStatement` according to the provided state variable name, discrete variable name, additional parameters, and a set of helper function symbols:
    - Applies variable substitution and body block processing.
    - Decomposes composite conditions (&&, ||).

# Arguments
- `statement::IfStatement`: The `IfStatement` node to be processed.
- `stateVarName`: The name of the state variable relevant to the normalization.
- `discrParamName`: The name of the discretization parameter.
- `stack`: Additional parameter(s) required for normalization.
- `helperFunSymSet`: A set of symbols representing helper functions used during normalization.

# Returns
- `Vector{IfStatement}`: A vector of normalized `IfStatement` objects resulting from the processing.

# Notes
This function is typically used as part of the IR normalization pipeline in the `SimpleModelIR` module.
"""
function expand_comparison(ex::Expr)
    # Check it's a comparison with exactly 5 args
    if ex.head == :comparison && length(ex.args) == 5
        a, op1, b, op2, c = ex.args
        # Only expand if both operators are the same and are <, >, <=, or >=
        if  op1 in (:<, :>, :<=, :>=) &&  op2  in (:<, :>, :<=, :>=)
            expanded_expr=Expr(:&&,
                Expr(:call, op1, a, b),
                Expr(:call, op2, b, c)
            )
            return expanded_expr
        end
    end
    return ex # leave untouched if not matching
end


function process_if_expr(statement,stateVarName,discrParamName,stack,helperFunSymSet)
    cond = statement.condition
    if !(cond isa Expr && ((cond.head == :call && cond.args[1] in [:>, :>=, :<, :<=,:(==)]) || cond.head in [:&&, :||] || (cond.head == :comparison)))
        return [statement]  # Already normalized  (there is ! not here)
    end
    if cond.head == :comparison
        cond = expand_comparison(cond)
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
    # === Apply rename_and_swap to whole condition ===
    kind, cond_zcfs,new_conds = process_if_condition(cond, stateVarName, discrParamName, stack, helperFunSymSet)
    #@show kind, cond_zcfs, new_conds

    # === Process IF body ===
    if_expr = statement.body
    processed_expr = process_if_block(if_expr.args[2], stateVarName, discrParamName, stack, helperFunSymSet)
    then_block = processed_expr.head == :block ? processed_expr : Expr(:block, processed_expr)
    if_expr.args[2] = then_block
    has_else = length(if_expr.args) == 3
    if has_else
        processed_expr = process_if_block(if_expr.args[3], stateVarName, discrParamName, stack, helperFunSymSet)
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
        push!(new_if_statements, IfStatement(zcfcond1, Expr(:if, deepcopy(cond1), then_block)))
        push!(new_if_statements, IfStatement(zcfcond2, Expr(:if, deepcopy(cond2), then_block)))
    elseif kind == :or && has_else
        cond1, cond2 = cond_zcfs
        zcfcond1, zcfcond2 = cond_zcfs
        cond1, cond2 = new_conds
        fallback1 = Expr(:if, Expr(:call, :!,  deepcopy(cond2)), else_block)
        push!(new_if_statements, IfStatement(zcfcond1, Expr(:if, deepcopy(cond1), then_block, fallback1)))

        fallback2 = Expr(:if, Expr(:call, :!, deepcopy(cond1)), else_block)
        push!(new_if_statements, IfStatement(zcfcond2, Expr(:if, deepcopy(cond2), then_block, fallback2)))
    elseif kind == :multi
        renamed_cond = rename_and_swap(statement.condition, stateVarName, discrParamName, :nothing, stack, helperFunSymSet)
        if_expr.args[1] = deepcopy(renamed_cond)
        for cond_i in cond_zcfs
            new_expr = Expr(:if, deepcopy(cond_i), if_expr, if_expr)
            push!(new_if_statements, IfStatement(cond_i, new_expr))
        end
    end
    return new_if_statements
end

