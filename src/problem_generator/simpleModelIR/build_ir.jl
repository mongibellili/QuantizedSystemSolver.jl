# ================================
# IR Building
# ================================
function parse_expr_tree(expr::Expr)
    statements = AbstractODEStatement[]
    for subexpr in expr.args
        if subexpr isa Expr && subexpr.head in [:(=), :+=, :-=, :*=, :/=]
            lhs = subexpr.args[1]
            rhs = subexpr.args[2]
            if lhs isa Expr && lhs.head == :tuple && rhs isa Expr && rhs.head == :tuple
                length(lhs.args) == length(rhs.args) || 
                    error("lhs and rhs must have the same length in $lhs = $rhs")

                for (l, r) in zip(lhs.args, rhs.args)
                    push!(statements, AssignStatement(l, r))
                end
            else
                # regular assignment
                if subexpr.head == :+=
                    rhs = Expr(:call, :+, lhs, rhs)
                elseif subexpr.head == :-=
                    rhs = Expr(:call, :-, lhs, rhs)
                elseif subexpr.head == :*=
                    rhs = Expr(:call, :*, lhs, rhs)
                elseif subexpr.head == :/=
                    rhs = Expr(:call, :/, lhs, rhs)
                end
                push!(statements, AssignStatement(lhs, rhs))
            end
     

  
        elseif subexpr isa Expr && subexpr.head == :if
            condition = subexpr.args[1]
            then_branch = subexpr.args[2]
            else_branch = length(subexpr.args) == 3 ? subexpr.args[3] : nothing
            if_expr = isnothing(else_branch) ?
                Expr(:if, condition, then_branch) :
                Expr(:if, condition, then_branch, else_branch)
            push!(statements, IfStatement(condition, if_expr))
        elseif subexpr isa Expr && subexpr.head == :while
            error("While loops are not supported in the current IR parser.")
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
                    error("Unsupported colon expression: $loop_iter")
                else
                    error("Unsupported colon expression: $loop_iter")
                end
            else
                error("Unsupported loop iterator format: $loop_iter")
            end
            local_statements = parse_expr_tree(Expr(:block, (subexpr.args[2] isa Expr && subexpr.args[2].head == :block ? subexpr.args[2].args : [subexpr.args[2]])...))
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
function problem_to_normalized_ir(expr::Expr, stateVarName::Symbol, discrParamName::Symbol, inline_mode::InlineMode)
    ir = build_ir(expr) # build IR
    stack = SymbolTableStack() # initialize the stack with one table
    normalized_ir_statements,numZC, numHelperF = normalize_ir(ir.statements,stack, stateVarName, discrParamName,:nothing,inline_mode) # normalize IR
    ir.statements=normalized_ir_statements
    return probInfo(ir,  numZC, numHelperF)
end

