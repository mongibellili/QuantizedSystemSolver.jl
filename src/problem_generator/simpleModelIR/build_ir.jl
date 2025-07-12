# ================================
# IR Building
# ================================
function parse_expr_tree(expr::Expr)
    statements = AbstractODEStatement[]
    for subexpr in expr.args
        if subexpr isa Expr && subexpr.head == :(=)
            push!(statements, AssignStatement(subexpr.args[1], subexpr.args[2]))
            
        elseif subexpr isa Expr && subexpr.head == :if
            condition = subexpr.args[1]
            then_branch = subexpr.args[2]
            else_branch = length(subexpr.args) == 3 ? subexpr.args[3] : nothing

            if_expr = isnothing(else_branch) ?
                Expr(:if, condition, then_branch) :
                Expr(:if, condition, then_branch, else_branch)

            push!(statements, IfStatement(condition, if_expr))


        elseif subexpr isa Expr && subexpr.head == :while
           # condition = subexpr.args[1]
           # body_exprs = subexpr.args[2] isa Expr && subexpr.args[2].head == :block ? subexpr.args[2].args : [subexpr.args[2]]
           # body_statements = parse_expr_tree(Expr(:block, body_exprs...))
           # push!(statements, WhileStatement(condition, body_statements))
        
        elseif subexpr isa Expr && subexpr.head == :for
            loop_var, loop_iter = subexpr.args[1].args
            # Determine loop type
            if loop_iter isa Expr && loop_iter.head == :call && loop_iter.args[1] == :(:)
                args = loop_iter.args[2:end]
                if length(args) == 2
                    loop_start, loop_end = args
                    loop_statement = nothing
                    loop_type = :range
                elseif length(args) == 3
                    #loop_start, loop_statement, loop_end = args
                    #loop_type = :range_with_statement
                    error("Unsupported colon expression: $loop_iter")
                else
                    error("Unsupported colon expression: $loop_iter")
                end
            else
                error("Unsupported loop iterator format: $loop_iter")
            end
            local_statements = parse_expr_tree(Expr(:block, (subexpr.args[2] isa Expr && subexpr.args[2].head == :block ?
                subexpr.args[2].args : [subexpr.args[2]])...))
            push!(statements, ForStatement(loop_var, loop_start, loop_end, local_statements, loop_statement, loop_iter, loop_type))
        #elseif subexpr isa Expr && subexpr.head == :function
        elseif subexpr isa Expr
        # generic fallback – should come **last**
            push!(statements, ExprStatement(subexpr))
        else
            push!(statements, ExprStatement(Expr(:quote, subexpr)))
        end
    end
    return statements
end




"""
    build_ir(expr::Expr)

Converts a symbolic expression into an ODEFunctionIR representation.
# Arguments
- `expr::Expr`: The symbolic expression to be converted into an ODEFunctionIR.
# Returns
- `ODEFunctionIR`: An intermediate representation of the ODE function, encapsulating the parsed
"""
function build_ir(expr::Expr)
    Base.remove_linenums!(expr)
    statements = parse_expr_tree(expr)
    return ODEFunctionIR(statements)
end




"""
    problem_to_normalized_ir(expr::Expr, stateVarName::Symbol, discrParamName::Symbol)

Converts a symbolic problem expression into a normalized intermediate representation (IR).

# Arguments
- `expr::Expr`: The symbolic expression representing the problem to be normalized.
- `stateVarName::Symbol`: The symbol representing the state variable in the problem.
- `discrParamName::Symbol`: The symbol representing the discrete variable in the problem.

# Returns
- `probInfo`: A structure containing the normalized IR and associated problem information.

# Description
This function processes the input symbolic expression, extracting relevant information and transforming it into a normalized IR suitable for further analysis or code generation. 
It uses the provided state variable and discretization parameter names to correctly interpret the structure of the problem.
this process is delegated to build_ir and normalize_ir functions.
"""
function problem_to_normalized_ir(expr::Expr, stateVarName::Symbol, discrParamName::Symbol)
    ir = build_ir(expr) # build IR
    return normalize_ir(ir, stateVarName, discrParamName)
end



#= function process_if_expr(statement,stateVarName,discrParamName,param,helperFunSymSet)
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
    elseif cond.head == :call && cond.args[1] in [:>, :>=,:(==)]
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
      #=   inner = Expr(:if,Expr(:call, :>, deepcopy(cond2),0), then_block) # deepcopy because inside events we do not want transformation to conditions in taylorEquationConstruction to be reflected back (pass by reference)
        outer = Expr(:if, Expr(:call, :>, deepcopy(cond1),0), inner) =#

        inner = Expr(:if,Expr(:call, :>, deepcopy(cond2),0), then_block) # deepcopy because inside events we do not want transformation to conditions in taylorEquationConstruction to be reflected back (pass by reference)
        outer = Expr(:if, Expr(:call, :>, deepcopy(cond1),0), inner)
        push!(new_if_statements, IfStatement(zcfcond1, outer))
        #  reverse the order
        inner_rev = Expr(:if, Expr(:call, :>, deepcopy(cond1),0), then_block)
        outer_rev = Expr(:if, Expr(:call, :>, deepcopy(cond2),0), inner_rev)
        push!(new_if_statements, IfStatement(zcfcond2, outer_rev))

    elseif kind == :and && has_else
        cond1, cond2 = cond_zcfs
        inner = Expr(:if, Expr(:call, :>, deepcopy(cond2),0), then_block)
        outer = Expr(:if, Expr(:call, :>, deepcopy(cond1),0), inner, else_block)
        push!(new_if_statements, IfStatement(cond1, outer))
        #  reverse the order
        inner_rev = Expr(:if, Expr(:call, :>, deepcopy(cond1),0), then_block)
        outer_rev = Expr(:if, Expr(:call, :>, deepcopy(cond2),0), inner_rev, else_block)
        push!(new_if_statements, IfStatement(cond2, outer_rev))

    elseif kind == :or && !has_else
        for c in cond_zcfs
            push!(new_if_statements, IfStatement(c, Expr(:if, Expr(:call, :>, deepcopy(c),0), then_block)))
        end
    elseif kind == :or && has_else
        c1, c2 = cond_zcfs
        fallback1 = Expr(:if, Expr(:call, :>, Expr(:call, :-, 0, deepcopy(c2)), 0), else_block)
        push!(new_if_statements, IfStatement(c1, Expr(:if, Expr(:call, :>, deepcopy(c1),0), then_block, fallback1)))

        fallback2 = Expr(:if, Expr(:call, :>, Expr(:call, :-, 0, deepcopy(c1)), 0), else_block)
        push!(new_if_statements, IfStatement(c2, Expr(:if, Expr(:call, :>, deepcopy(c2),0), then_block, fallback2)))
    elseif kind == :multi
        renamed_cond = changeVarNames_params(statement.condition, stateVarName, discrParamName, :nothing, param, helperFunSymSet)
        if_expr.args[1] = deepcopy(renamed_cond)
        for cond_i in cond_zcfs
            new_expr = Expr(:if, deepcopy(cond_i), if_expr, if_expr)
            push!(new_if_statements, IfStatement(cond_i, new_expr))
        end
    end

    return new_if_statements
end =#
