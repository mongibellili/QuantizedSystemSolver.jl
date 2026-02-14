

to_expr(a::AssignStatement) = Expr(:(=), a.lhs, a.rhs)


to_expr(a::IfStatement)     = Expr(:if, a.condition, to_expr(a.body), nothing)
to_expr(a::ForStatement)    = Expr(:for,
                                   Expr(:in, a.var, Expr(:colon, a.start, a.stop)),
                                   to_expr(a.body))
to_expr(a::WhileStatement)  = Expr(:while, a.condition, to_expr(a.body))
to_expr(a::ExprStatement)   = a.expr
# Fallback for things already expressions or basic values
to_expr(x::Expr)    = x
to_expr(x::Symbol)  = x
to_expr(x::Number)  = x
to_expr(x) = Meta.parse(string(x))



"""
    createExactJacFun(finalJacFunctionCode::Expr,Exactjac:: Dict{Expr,Union{Float64,Int,Symbol,Expr}},funName::Symbol,f) 

 constructs the exact jacobian entries as a function from the existing dictionary Exactjac resulted from extractJacDep and extractJacDepLoop functions. \n

  From the input dictionary exactJac, we see that the key is always a expression that hold a tuple of size 2: the first element is going to be the index 'i' and the second element if the index 'j'. The value corresponding to this key, which is the exact jacobian entry, is put in a cache. i.e the function maps the keys of the dictionary to their values using an 'if-statement.. This approach does not depend on the size of the problem.
    # Example:   
```jldoctest
using QuantizedSystemSolver
exacteJacExpr=Dict{Expr,Union{Float64,Int,Symbol,Expr}}(:((1, 1)) => :(-2.0 * (q[2])[0]), :((1, 2)) => :(1 - 2.0 * (q[1])[0]),:(((2, 9), i - 1)) => :((q[i])[0]), :(((2, 9), i)) => :((q[i - 1])[0]),:((10, 10)) => -1, :((10, 1)) => 1);
exactJac=QuantizedSystemSolver.createExactJacFun(:(),exacteJacExpr,:f,0);
exactJac

# output

:(function exactJacf(q::Vector{Taylor0}, p::Vector{Float64}, cache::AbstractVector{Float64}, i::Int, j::Int, t::Float64, f_)
      (if i == 0
              return nothing
          elseif i == 1 && j == 1
              cache[1] = -2.0 * (q[2])[0]
              return nothing
          elseif i == 10 && j == 10
              cache[1] = -1
              return nothing
          elseif 2 <= i <= 9 && j == i - 1
              cache[1] = (q[i])[0]
              return nothing
          elseif i == 1 && j == 2
              cache[1] = 1 - 2.0 * (q[1])[0]
              return nothing
          elseif 2 <= i <= 9 && j == i
              cache[1] = (q[i - 1])[0]
              return nothing
          elseif i == 10 && j == 1
              cache[1] = 1
              return nothing
          end,)
  end)
```
"""


function createExactJacFun(modelHelperCode::Expr,Exactjac:: Dict{Expr,ScopedEquation},funName::Symbol) 
    finalJacFunctionCode=copy(modelHelperCode)#start by kept helper assignments. copy to avoid changing the original modelHelperCode (needed in diffeqfunction)

    branches = []
    #push!(branches, :(i == 0) => :(return nothing))  # first branch

    for (key_expr, scoped) in Exactjac
        cond_expr = nothing
        body_exprs = []

        # Condition: either (i == X && j == Y) or (range <= i <= range && j == Y)
        if key_expr.args[1] isa Int
            cond_expr = :((i == $(key_expr.args[1])) && (j == $(key_expr.args[2])))
        elseif key_expr.args[1] isa Expr
            start_val = key_expr.args[1].args[1]
            stop_val  = key_expr.args[1].args[2]
            cond_expr = :(($start_val <= i <= $stop_val) && (j == $(key_expr.args[2])))
        else
            error("Unsupported key type in Exactjac: $key_expr")
        end

        # Add helper assignments if any
        for k in scoped.helperAssignments
            if !(k isa AssignStatement) || k.keep_assignment
                push!(body_exprs, to_expr(k))
            end
        end

        # Main assignment: cache[1] = (rhs)
        push!(body_exprs, Expr(:(=), Expr(:ref, :cache, 1), to_expr(scoped.eqs_RHS)))
        push!(body_exprs, :(return nothing))

        push!(branches, cond_expr => Expr(:block, body_exprs...))
    end

    # Else branch: error if no match
    #else_branch = :(error("Invalid i,j: ", i, ", ", j))

    ex = ex = :(nothing)#else_branch
    for (cond, body) in reverse(branches)
        ex = Expr(:if, cond, body, ex)
    end

    push!(finalJacFunctionCode.args, ex)
 


  
  Base.remove_linenums!(finalJacFunctionCode)
  def1=Dict{Symbol,Any}() 
  def1[:head] = :function
  def1[:name] = Symbol(:exactJac,funName)  
  #def1[:args] = [:(q::Vector{Taylor0}),:(p::Vector{Float64}),:(cache::AbstractVector{Float64}),:(i::Int),:(j::Int),:(t::Float64),:(f_)]
  def1[:args] = [:(q::Vector{Taylor0}),:(p),:(cache::AbstractVector{Float64}),:(i::Int),:(j::Int),:(t::Float64),:(f_)]
  def1[:body] = finalJacFunctionCode
  functioncode1=combinedef(def1)
end
 


"""
    build_flat_if(varsym::Symbol, branches::Vector{Any})

Builds an AST of the form:

    if varsym == 1
        branch1
        return nothing
    elseif varsym == 2
        branch2
        return nothing
    ...
    else
        nothing
    end
"""
function build_flat_if(varsym::Symbol, branches_vect::Vector)
    pairs = []
    for (idx, rhs) in enumerate(branches_vect)
        cond_expr = :( $varsym == $idx )
        body_expr = Expr(:block, to_expr(rhs), :(return nothing))
        push!(pairs, cond_expr => body_expr)
    end

    # default else branch
    ex = :(nothing)
    for (cond, body) in reverse(pairs)
        ex = Expr(:if, cond, body, ex)
    end
    return ex
end


function createEqFun(modelHelperCode::Expr,equs::Dict{Union{Int,Symbol,Expr},ScopedEquation},zcequs::Vector{Expr},eventequs::Vector{Expr},fname::Symbol)# where{F}
    finalFunctionCode=copy(modelHelperCode)#start by kept helper assignments. copy to avoid changing the original modelHelperCode (needed in exactJacfunction)
    branches = []

    for (condkey, scoped) in equs
        cond_expr = nothing
        body_exprs = []

        if condkey isa Int
            cond_expr = :(i == $condkey)
        elseif condkey isa Expr && condkey.head == :tuple  # i.e. :(2:3)
            start_val = condkey.args[1]
            stop_val  = condkey.args[2]
            cond_expr = :($start_val <= i <= $stop_val)
        elseif condkey isa Symbol   # check and remove this line later, symbol never visited because varnum gets converted to int in normalize_ir.jl line 170
            cond_expr = :(i == $condkey)
        else
            error("Unsupported key type: $condkey")
        end

        # Add helper assignments if any
        for k in scoped.helperAssignments
            if !(k isa AssignStatement) || k.keep_assignment
                push!(body_exprs, to_expr(k))
            end
        end

        # Add main equation RHS
        push!(body_exprs, to_expr(scoped.eqs_RHS))

        push!(body_exprs, :(return nothing))

        push!(branches, cond_expr => Expr(:block, body_exprs...))
    end

    # Build flat if/elseif chain
    ex = :(nothing)  # default else
    for (cond, body) in reverse(branches)
        ex = Expr(:if, cond, body, ex)
    end
    push!(finalFunctionCode.args, ex)

    ############## ZCF ###################
    if !isempty(zcequs)
        myex2 = build_flat_if(:zc, zcequs)
        push!(finalFunctionCode.args, myex2)
    end

    ############## events ################
    if !isempty(eventequs)
        myex3 = build_flat_if(:ev, eventequs)
        push!(finalFunctionCode.args, myex3)
    end
 
  Base.remove_linenums!(finalFunctionCode)
  def=Dict{Symbol,Any}()
  def[:head] = :function
  def[:name] = fname  
  #def[:args] = [:(i::Int),:(zc::Int),:(ev::Int),:(q::Vector{Taylor0}),:(p::Vector{Float64}), :(t::Taylor0),:(cache::Vector{Taylor0}),:(f_)]  
  def[:args] = [:(i::Int),:(zc::Int),:(ev::Int),:(q::Vector{Taylor0}),:(p), :(t::Taylor0),:(cache::Vector{Taylor0}),:(f_)]  
  def[:body] = finalFunctionCode 
  functioncode=combinedef(def)
 # @show functioncode

end
