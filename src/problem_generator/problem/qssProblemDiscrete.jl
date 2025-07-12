
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
    equs=Dict{Union{Int,Expr},Union{Int,Symbol,Expr}}()
    jac = Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}()# set used because do not want to re-insert an existing varNum
    dD = Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}() # like jac but for discrete variables
    exacteJacExpr = Dict{Expr,Union{Float64,Int,Symbol,Expr}}()
    zcequs=Vector{Expr}()#vect to collect the conditions of if-statements (zero-crossing)
    eventequs=Vector{Expr}()#vect to collect events (effect)
    ZCjac=Vector{Vector{Int}}() # the dependency of zero crossing to variables
    ZCFCounter=0
    SZ= Dict{Int,Set{Int}}() # opposite of ZCjac: the effect of a continuous variable on which zc
    dZ= Dict{Int,Set{Int}}()  # the effect of a discrete variable on which zc
    evsArr = EventDependencyStruct[]
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
                    rhs = extractJacDepNormal(varNum, rhs, jac, exacteJacExpr, jac_mode, symDict, dD)
                    equs[varNum] = transformFSimplecase(rhs)
                elseif rhs isa Expr && rhs.head == :call && is_custom_function(rhs)
                    rhs = extractJacDepNormal(varNum, rhs, jac, exacteJacExpr, jac_mode, symDict, dD)
                    equs[varNum] = transformFSimplecase(rhs)
                else
                    rhs = extractJacDepNormal(varNum, rhs, jac, exacteJacExpr, jac_mode, symDict, dD)
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
                specRHS = extractJacDepLoop(b, niter, specRHS, jac, exacteJacExpr, jac_mode, symDict, dD)
                equs[:(($b, $niter))] = transformFSimplecase(specRHS)
            elseif specRHS isa Expr && specRHS.head == :call && is_custom_function(specRHS)
                specRHS = extractJacDepLoop(b, niter, specRHS, jac, exacteJacExpr, jac_mode, symDict, dD)
                equs[:(($b, $niter))] = transformFSimplecase(specRHS)
            else
                specRHS = extractJacDepLoop(b, niter, specRHS, jac, exacteJacExpr, jac_mode, symDict, dD)
                temp = (transformF(:($(specRHS), 1))).args[2]
                num_cache_equs = max(num_cache_equs, temp)
                equs[:(($b, $niter))] = specRHS
            end

        elseif statement isa IfStatement
            
            zcf = statement.condition
            ZCFCounter += 1
            extractZCJacDep(ZCFCounter, statement.condition, ZCjac, SZ, dZ)

            if zcf isa Expr && zcf.head == :ref
                push!(zcequs, transformFSimplecase(zcf))
            else
                temp = (transformF(:($(zcf), 1))).args[2]
                num_cache_equs = max(num_cache_equs, temp)
                push!(zcequs, zcf)
            end

            # You may want to rewrite `handleEvents` to take IR-form if-bodies
           # @show statement.body
            handleEvents(statement.body, eventequs, length(zcequs), evsArr)

        elseif statement isa ExprStatement && statement.expr.head == :function
            push!(functionEx.args, statement.expr)
        elseif statement isa ExprStatement
            push!(otherCode.args, statement.expr)
        end
    end

    #@show eventequs

      #println("after second pass, otherCode=$otherCode , equs=$equs, jac=$jac, dD=$dD, exacteJacExpr=$exacteJacExpr, zcequs=$zcequs, eventequs=$eventequs, ZCjac=$ZCjac, SZ=$SZ, dZ=$dZ")
    if length(functionEx.args)==0
        closurefunc=0
        diffEqfunctionExpression=createDiscEqFun(otherCode,equs,zcequs,eventequs,fname,closurefunc)# diff equations before this are stored in a dict:: now we have a giant function that holds all diff equations
        exactJacfunction=createExactJacFun(otherCode,exacteJacExpr,fname,closurefunc) 
    elseif length(functionEx.args)==1
        if jac_mode==:symbolic
            @warn("symbolic jacobian mode is not advised with helper functions. Explicitly use jac_mode = :approximate in ODEProblem.")   
        end
        closurefunc=@RuntimeGeneratedFunction(functionEx.args[1]) 
        diffEqfunctionExpression=createDiscEqFun(otherCode,equs,zcequs,eventequs,fname,closurefunc)# diff equations before this are stored in a dict:: now we have a giant function that holds all diff equations
        exactJacfunction=createExactJacFun(otherCode,exacteJacExpr,fname,closurefunc)
    else
        error("Error: Currently only one helper function is allowed inside the model function. Consider moving them to the top level.")
    end
 
    if numHelperFunCalls>length(functionEx.args)# the case when the user has defined a helper function inside the model function and one helperF outside and only called the outside helperF is considered a user mistake (either user made a typo or forgot to remove).
        
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

