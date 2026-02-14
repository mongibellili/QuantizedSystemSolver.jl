


"""
    process_du_rhs(rhs::Union{Number,Symbol,Expr}, b::Int, niter::Int, 
                   localHelperAssignments::Vector{AbstractODEStatement}, 
                   equs::Dict{Union{Int,Symbol,Expr},ScopedEquation},
                   exactJacExpr::Dict{Expr,ScopedEquation},
                   symDict::Dict{Symbol,Expr},
                   jac::Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},
                   dD::Dict{Union{Int,Symbol,Expr},Set{Union{Int,Symbol,Expr}}},
                   jac_mode::Symbol)

Process the right-hand side (RHS) of a differential equation term.

# Arguments
- `rhs::Union{Number,Symbol,Expr}`: The right-hand side expression to process
- `b::Int`: Block index parameter if equation inside a loop or normal index if not inside a loop
- `niter::Int`: Number of iterations in case there is a loop
- `localHelperAssignments::Vector{AbstractODEStatement}`: Vector of helper assignment statements to this equation
- `equs::Dict{Union{Int,Symbol,Expr},ScopedEquation}`: Dictionary of scoped equations
- `exactJacExpr::Dict{Expr,ScopedEquation}`: Dictionary of exact Jacobian expressions
- `symDict::Dict{Symbol,Expr}`: Symbol to expression mapping dictionary
- `jac::Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}`: Jacobian dependency dictionary
- `dD::Dict{Union{Int,Symbol,Expr},Set{Union{Int,Symbol,Expr}}}`: Derivative dependency dictionary
- `jac_mode::Symbol`: Jacobian computation mode

# Returns
Processed representation of the differential equation term RHS.
"""
function process_du_rhs(rhs::Union{Number,Symbol,Expr},b::Int,niter::Int,localHelperAssignments::Vector{AbstractODEStatement}, equs::Dict{Union{Int,Symbol,Expr},ScopedEquation},exactJacExpr :: Dict{Expr,ScopedEquation},symDict::Dict{Symbol,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},dD :: Dict{Union{Int,Symbol,Expr},Set{Union{Int,Symbol,Expr}}},jac_mode::Symbol)
    newCacheSize = 1 # default cache size
    if rhs isa Number || rhs isa Symbol
        rhs=transformFSimplecase(rhs)
    elseif rhs isa Expr && (rhs.head == :ref  || (rhs.head == :call && is_custom_function(rhs)))
        jacset=extractJacDep(b,niter, rhs, jac,  dD)
        if jac_mode==:symbolic extractJacExpression(b,niter, rhs, jacset, exactJacExpr, localHelperAssignments,symDict) end
        rhs=transformFSimplecase(rhs)
    else
        
        jacset=extractJacDep(b,niter, rhs, jac,  dD)
        if jac_mode==:symbolic extractJacExpression(b,niter, rhs, jacset, exactJacExpr, localHelperAssignments,symDict) end
        newCacheSize = (transformF!(:($(rhs), 1)))# note the expression contains rhs and 1 ...different from (:$(rhs),1)
    end
    if b==-1
        equs[niter] =  ScopedEquation(localHelperAssignments,rhs) # each equation is stored in a dict with its lhs (varNum or for loop index) as key and its rhs and helper assignments as value. 
    else
        equs[:(($b, $niter))] =  ScopedEquation(localHelperAssignments,rhs)
    end
    return newCacheSize
end


"""
    odeProblemFunc(ir::ODEFunctionIR,::Val{T},::Val{D},::Val{Z},initCond::Vector{Float64},discrVars::Union{Vector{EM}, Tuple{Vararg{EM}}},preProcessData::PreProcessData,jac_mode ::Symbol) where {T,D,Z,EM} 

        continues building a discrete problem. 
    It receives an expression and useful info from the main interface. It calls the transform function from the taylorEquationConstruction.jl file to change the AST of all operations to personlized ones and update the needed cache size. It also construct via helper functions the Exact jacobian function, the jacobian dependecy (jac) and the state-derivative dependency (SD:opposite of jacobian), the state to zero-crossing dependency (SZ) and events to derivative and zero-crossing (HD and HZ) as vectors. Finally, it groups all differential equations and events in one function, and constructs a discrete problem from the qssProblemDefinition.jl file.
"""
function odeProblemFunc(ir::ODEFunctionIR,::Val{T},::Val{D},::Val{Z},initCond::Vector{Float64},discrVars::Union{Vector{EM}, Tuple{Vararg{EM}}},preProcessData::PreProcessData,jac_mode ::Symbol) where {T,D,Z,EM} 
    # used EM to account for when problem contains if-statements whithout discrete vars
    #println("after first pass, odeexpr=$odeExprs")
    du=preProcessData.du
    tspan= preProcessData.tspan
    fname= preProcessData.prbName
    mod= preProcessData.mod 
    is_top_level= preProcessData.is_top_level
    numHelperFunCalls=preProcessData.numHelperFunCalls
    symDict=Dict{Symbol,Expr}()
    equs=Dict{Union{Int,Symbol,Expr},ScopedEquation}()
    jac = Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}()# inside-set used because do not want to re-insert an existing varNum
    dD =  Dict{Union{Int,Symbol,Expr},Set{Union{Int,Symbol,Expr}}}() # like jac but for discrete variables (and also dependency in the other direction)
    exactJacExpr = Dict{Expr,ScopedEquation}()
    zcequs=Vector{Expr}()#vect to collect the conditions of if-statements (zero-crossing)
    eventequs=Vector{Expr}()#vect to collect events (effect)
    ZCjac=Vector{Vector{Int}}() # the dependency of zero crossing to variables
    ZCFCounter=0
    SZ= Dict{Int,Set{Int}}() # opposite of ZCjac: the effect of a continuous variable on which zc
    dZ= Dict{Int,Set{Int}}()  # the effect of a discrete variable on which zc
    evsArr = EventDependencyStruct[]
    num_cache_equs=1#cachesize
    functionEx=Expr(:block) #function to be inserted inside the main function
    modelHelperCode=Expr(:block) #other code to be inserted inside the main function
    helperAssignments=AbstractODEStatement[] # vect of helper assignments to be used in model scope. kept empty because modelHelperCode already does that. it is used merely to match function signature in handling other local scopes
    for statement in ir.statements
        if statement isa AssignStatement
            lhs = statement.lhs
            rhs = statement.rhs
            keep_assignment = statement.keep_assignment         
            if lhs isa Expr && lhs.head == :ref && lhs.args[1] == du
                varNum = lhs.args[2]
                newCacheSize=process_du_rhs(rhs,-1,varNum,helperAssignments, equs,exactJacExpr,symDict,jac,dD,jac_mode)
                num_cache_equs = max(num_cache_equs, newCacheSize)
            else
                if keep_assignment# If the assignment is kept, add it to modelHelperCode
                    push!(modelHelperCode.args, Expr(:(=), lhs, rhs))
                end
            end
        elseif statement isa ForStatement
            localHelperAssignments=AbstractODEStatement[] 
            b = statement.start
            niter = statement.stop
            for bodyElement in statement.body
                if bodyElement isa AssignStatement
                    lhs = bodyElement.lhs
                    rhs = bodyElement.rhs
                    keep_assignment = bodyElement.keep_assignment
                    if lhs isa Expr && lhs.head == :ref && lhs.args[1] == du
                        newCacheSize=process_du_rhs(rhs,b,niter,localHelperAssignments, equs,exactJacExpr,symDict,jac,dD,jac_mode)
                        num_cache_equs = max(num_cache_equs, newCacheSize)
                        break # currently only one diff eq inside one for loop
                    else
                        if keep_assignment
                            push!(localHelperAssignments, bodyElement)
                        end
                    end
                else
                    push!(localHelperAssignments, bodyElement)
                end
            end

        elseif statement isa IfStatement    
            zcf = statement.condition
            ZCFCounter += 1
            extractZCJacDep(ZCFCounter, statement.condition, ZCjac, SZ, dZ)
            if zcf isa Symbol || (zcf isa Expr && zcf.head == :ref)
                push!(zcequs, transformFSimplecase(zcf))
            else
                newCacheSize = (transformF!(:($(zcf), 1)))#.args[2]
                num_cache_equs = max(num_cache_equs, newCacheSize)
                push!(zcequs, zcf)
            end
            handleEvents(statement.body, eventequs, length(zcequs), evsArr)
        elseif statement isa ExprStatement && statement.expr.head == :function
            push!(functionEx.args, statement.expr)
        elseif statement isa ExprStatement
            push!(modelHelperCode.args, statement.expr)
        end
    end
    closurefunc=0
    if length(functionEx.args)==0   # helper functions inside the model
    elseif length(functionEx.args)==1
        closurefunc=@RuntimeGeneratedFunction(functionEx.args[1]) 
        jac_mode==:symbolic && @warn("symbolic jacobian mode is not advised with helper functions. use jac_mode = :approximate in ODEProblem.")   
    else
        error("Error: Currently only one helper function is allowed inside the model function. Consider moving them to the top level.")
    end
    if jac_mode==:symbolic        
        exactJacfunction=createExactJacFun(modelHelperCode,exactJacExpr,fname) 
    else
        exactJacfunction=createExactJacFun(Expr(:block),exactJacExpr,fname) 
    end
    diffEqfunctionExpression=createEqFun(modelHelperCode,equs,zcequs,eventequs,fname)   # diff equations before this are stored in a dict:: now we have a giant function that holds all diff equations


    if numHelperFunCalls>length(functionEx.args)# helper functions outside the model. the case when the user has defined a helper function inside the model function and one helperF outside and only called the outside helperF is considered a user mistake (either user made a typo or forgot to remove).
        jac_mode==:symbolic && @warn("symbolic jacobian mode is not advised with helper functions. use jac_mode = :approximate in ODEProblem.")  
        if is_top_level                    
            RuntimeGeneratedFunctions.init(Main)
            exactJacfunctionF=RuntimeGeneratedFunction(mod, mod,exactJacfunction)
            diffEqfunctionF=RuntimeGeneratedFunction(mod, mod,diffEqfunctionExpression)   
        else
            @warn("Simulation is not running on top level with helper functions. This may hurt performance. Consider placing your code on top level, or placing the helper function inside the model.")
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
    dDVect =createdDVect(dD,Val(D))

    SZVect=createSZVect(SZ,Val(T))
    # temporary dependencies to be used to determine HD and HZ...determine HD: event-->Derivative   && determine HZ:Event-->ZCfunction....influence of events on derivatives and zcfunctions:
    #an event is a discrteVar change or a cont Var change. So HD=HD1 UNION HD2  (same for HZ=HZ1 UNION HZ2)
    # (1) through a discrete Var: 
    #  ============================
    # HD1=Hd-->dD  where Hd  comes from the eventDependecies Struct. 
    # HZ1=Hd-->dZ  where Hd  comes from the eventDependecies Struct. 
    HZ1HD1=createDependencyToEventsDiscr(dDVect,dZ,evsArr) 
    # (2) through a continous Var: 
    # ==============================
     #  HD2=Hs-->sD where  Hs comes from the eventDependecies Struct.(sd already created) 
     #  HZ2=Hs-->sZ where  Hs comes from the eventDependecies Struct.(sZ already created) 
    HZ2HD2=createDependencyToEventsCont(SDVect,SZ,evsArr) 
    ##########UNION##############
    HZ=unionDependency(HZ1HD1[1],HZ2HD2[1])
    HD=unionDependency(HZ1HD1[2],HZ2HD2[2])
    # mapFun=createMapFun(jac,fname)
    # mapFunF=RuntimeGeneratedFunction(mapFun)
    myodeProblem = ODEDiscProblem(fname,Val(jac_mode),Val(T),Val(D),Val(Z),Val(num_cache_equs),initCond, collect(discrVars), jacVect ,ZCjac  ,diffEqfunctionF, evsArr,SDVect,HZ,HD,SZVect,exactJacfunctionF,tspan,[closurefunc])
end

