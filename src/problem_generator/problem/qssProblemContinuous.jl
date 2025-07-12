
"""
    odeProblemFunc(ir::ODEFunctionIR,::Val{T},::Val{D},::Val{0},initCond::Vector{Float64},discrVars::Union{Vector{EM}, Tuple{Vararg{EM}}},preProcessData::PreProcessData,jac_mode ::Symbol) where {T,D,EM} 

This function continues building a continuous problem. it receives an expression and useful info from the main interface. it calls the transform function from the taylorEquationConstruction.jl file to change the AST of all operations to personlized ones and update the needed cache size. It also construct via helper functions the Exact jacobian function, the jacobian dependecy and the state-derivative dependency (opposite of jacobian) as vectors. Finally, it groups all differential equations in one function, and constructs a continous problem from the qssProblemDefinition.jl file.
# Arguments
- `odeExprs::Expr`: The expression of the whole user code in the function defining the problem with names modified and parameters plugged in.
- `Val{T}`: the dimensions of the system of differential equations.    
- `Val{0}`: No zero-crossing functions. pure continous problem.   
- `Val{0}`: No events functions. pure continous problem.  
- `initConditions::Vector{Float64}`: No zero-crossing functions. pure continous problem.  
- `du::Symbol`: to distinguish the start of a differential equations.  
- `symDict::Dict{Symbol,Expr}`: maps a reference expression to a symbol (qi->q[i]).  
- `tspan::Tuple{Float64, Float64}`: stores the initial time and final time of the simulation.  
- `prbName::Symbol`: The problem name as chosen by the user to be carried to the solution for displaying purposes.  
"""

function odeProblemFunc(ir::ODEFunctionIR,::Val{T},::Val{D},::Val{0},initCond::Vector{Float64},discrVars::Union{Vector{EM}, Tuple{Vararg{EM}}},preProcessData::PreProcessData,jac_mode ::Symbol) where {T,D,EM} 
    du=preProcessData.du
    tspan= preProcessData.tspan
    fname= preProcessData.prbName
    mod= preProcessData.mod
    is_top_level= preProcessData.is_top_level
    numHelperFunCalls=preProcessData.numHelperFunCalls
    symDict=Dict{Symbol,Expr}()
    equs=Dict{Union{Int,Expr},Union{Int,Symbol,Expr}}()
    jac = Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}()# set used because do not want to re-insert an existing varNum
    dD = Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}() # Not needed for continuous. like jac but for discrete variables: kept for continuous probl to have common function signature
    exacteJacExpr = Dict{Expr,Union{Float64,Int,Symbol,Expr}}()

    num_cache_equs=1#cachesize
    functionEx=Expr(:block) #function to be inserted inside the main function
    otherCode=Expr(:block) #other code to be inserted inside the main function
    for statement in ir.statements
        if statement isa AssignStatement
            lhs = statement.lhs
            rhs = statement.rhs

            if lhs isa Expr && lhs.head == :ref && lhs.args[1] == du
                varNum = lhs.args[2]
                if rhs isa Number || rhs isa Symbol
                    equs[varNum] = transformFSimplecase(rhs)
                elseif rhs isa Expr && rhs.head == :ref
                    rhs = extractJacDepNormal(varNum, rhs, jac, exacteJacExpr, jac_mode, symDict,dD)
                    equs[varNum] = transformFSimplecase(rhs)
                elseif rhs isa Expr && rhs.head == :call && is_custom_function(rhs)
                    rhs = extractJacDepNormal(varNum, rhs, jac, exacteJacExpr, jac_mode, symDict,dD)
                    equs[varNum] = transformFSimplecase(rhs)
                else
                    rhs = extractJacDepNormal(varNum, rhs, jac, exacteJacExpr, jac_mode, symDict,dD)
                    temp = (transformF(:($(rhs), 1))).args[2]
                    num_cache_equs = max(num_cache_equs, temp)
                    equs[varNum] = rhs
                end

            elseif lhs isa Symbol && rhs isa Expr && rhs.head in [:call, :ref]
                # Already processed in normalize_ir
                continue
            elseif lhs isa Expr && lhs.head == :tuple && rhs isa Symbol
                # Already processed
                continue
            elseif lhs isa Symbol && rhs isa Number
                continue
            else
                push!(otherCode.args, Expr(:(=), lhs, rhs))
            end

        elseif statement isa ForStatement
            b = statement.start
            niter = statement.stop
            specRHS = statement.body[1] isa AssignStatement ? statement.body[1].rhs : nothing  # Simplified assumption

            if specRHS isa Number || specRHS isa Symbol
                equs[:(($b, $niter))] = transformFSimplecase(specRHS)
            elseif specRHS isa Expr && specRHS.head == :ref
                specRHS = extractJacDepLoop(b, niter, specRHS, jac, exacteJacExpr, jac_mode, symDict,dD)
                equs[:(($b, $niter))] = transformFSimplecase(specRHS)
            elseif specRHS isa Expr && specRHS.head == :call && is_custom_function(specRHS)
                specRHS = extractJacDepLoop(b, niter, specRHS, jac, exacteJacExpr, jac_mode, symDict,dD)
                equs[:(($b, $niter))] = transformFSimplecase(specRHS)
            else
                specRHS = extractJacDepLoop(b, niter, specRHS, jac, exacteJacExpr, jac_mode, symDict,dD)
                temp = (transformF(:($(specRHS), 1))).args[2]
                num_cache_equs = max(num_cache_equs, temp)
                equs[:(($b, $niter))] = specRHS
            end

        elseif statement isa ExprStatement && statement.expr.head == :function
            push!(functionEx.args, statement.expr)

        elseif statement isa ExprStatement
            push!(otherCode.args, statement.expr)
        end
    end

      #println("after second pass, otherCode=$otherCode , equs=$equs, jac=$jac, dD=$dD, exacteJacExpr=$exacteJacExpr, zcequs=$zcequs, eventequs=$eventequs, ZCjac=$ZCjac, SZ=$SZ, dZ=$dZ")
    if length(functionEx.args)==0
        closurefunc=0
        diffEqfunctionExpression=createContEqFun(otherCode,equs,fname,closurefunc)# diff equations before this are stored in a dict:: now we have a giant function that holds all diff equations
        exactJacfunction=createExactJacFun(otherCode,exacteJacExpr,fname,closurefunc) 
    elseif length(functionEx.args)==1
        if jac_mode==:symbolic
            @warn("symbolic jacobian mode is not advised with helper functions. Explicitly use jac_mode = :approximate in ODEProblem.")   
        end
        closurefunc=@RuntimeGeneratedFunction(functionEx.args[1]) 
        diffEqfunctionExpression=createContEqFun(otherCode,equs,fname,closurefunc)# diff equations before this are stored in a dict:: now we have a giant function that holds all diff equations
        exactJacfunction=createExactJacFun(otherCode,exacteJacExpr,fname,closurefunc)
    else
        error("Error: Currently only one helper function is allowed inside the model function. Consider moving them to the top level.")
    end
 
    if numHelperFunCalls>length(functionEx.args)# It is considered a user mistake the case when the user has defined a helper function inside the model function and one helperF outside and only called the outside helperF  (either user made a typo or forgot to remove).
        if jac_mode==:symbolic
            @warn("symbolic jacobian mode is not advised with helper functions. Explicitly use jac_mode = :approximate in ODEProblem.")   
        end
        if is_top_level                     
            RuntimeGeneratedFunctions.init(Main)
            exactJacfunctionF=RuntimeGeneratedFunction(mod, mod,exactJacfunction)
            diffEqfunctionF=RuntimeGeneratedFunction(mod, mod,diffEqfunctionExpression)   
        else
            @warn("Simulation is not running on same level with helper functions. This may hurt performance. Consider placing the solve function on top level, or placing the helper function inside the model.")
            exactJacfunctionF1 = Base.eval(mod, exactJacfunction)
            exactJacfunctionF = (args...) -> Base.invokelatest(exactJacfunctionF1, args...)
            diffEqfunctionF1 = Base.eval(mod, diffEqfunctionExpression)
            diffEqfunctionF = (args...) -> Base.invokelatest(diffEqfunctionF1, args...)
        end
    else
        exactJacfunctionF=@RuntimeGeneratedFunction(exactJacfunction)
        diffEqfunctionF=@RuntimeGeneratedFunction(diffEqfunctionExpression) 
    end

    jacVect=createJacVect(jac,Val(T))
    SDVect=createSDVect(jac,Val(T))
   
    myodeProblem = ODEContProblem(fname,Val(jac_mode),Val(T),Val(D),Val(0),Val(num_cache_equs),initCond, collect(discrVars), jacVect   ,diffEqfunctionF,SDVect,exactJacfunctionF,tspan,[closurefunc])
end

