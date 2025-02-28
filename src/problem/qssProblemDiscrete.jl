
"""
    NLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{D},::Val{Z},initCond::Vector{Float64},du::Symbol,tspan::Tuple{Float64, Float64},discrVars::Union{Vector{EM}, Tuple{Vararg{EM}}},prbName::Symbol) where {T,D,Z,EM} 

        continues building a discrete problem. 
    It receives an expression and useful info from the main interface. It calls the transform function from the taylorEquationConstruction.jl file to change the AST of all operations to personlized ones and update the needed cache size. It also construct via helper functions the Exact jacobian function, the jacobian dependecy (jac) and the state-derivative dependency (SD:opposite of jacobian), the state to zero-crossing dependency (SZ) and events to derivative and zero-crossing (HD and HZ) as vectors. Finally, it groups all differential equations and events in one function, and constructs a discrete problem from the qssProblemDefinition.jl file.
"""
function NLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{D},::Val{Z},initCond::Vector{Float64},du::Symbol,tspan::Tuple{Float64, Float64},discrVars::Union{Vector{EM}, Tuple{Vararg{EM}}},prbName::Symbol) where {T,D,Z,EM} 
    # used EM to account for when problem contains if-statements whithout discrete vars
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
    functionEx=quote end #function to be inserted inside the main function
    otherCode=quote end #other code to be inserted inside the main function
    for argI in odeExprs.args
        if argI isa Expr &&  argI.head == :(=) && argI.args[1]== :discrete # old approach use the keyword discrete to define discrete variables
            discrVars = Vector{Float64}(argI.args[2].args)
         #only diff eqs: du[]= number/one ref/call  
        elseif argI isa Expr &&  argI.head == :(=)  && argI.args[1] isa Expr && argI.args[1].head == :ref && argI.args[1].args[1]==du
            rhs=argI.args[2];varNum=argI.args[1].args[2] # varnum is the order/index of variable
            if rhs isa Number || rhs isa Symbol # rhs of equ =number  or symbol
                equs[varNum]=:($((transformFSimplecase(:($(rhs))))))
            elseif rhs isa Expr && rhs.head==:ref 
                extractJacDepNormalDiscrete(varNum,rhs,jac,exacteJacExpr ,symDict,dD )
                equs[varNum ]=:($((transformFSimplecase(:($(rhs))))))
            else #rhs head==call              
                extractJacDepNormalDiscrete(varNum,rhs,jac,exacteJacExpr ,symDict,dD )
                temp=(transformF(:($(rhs),1))).args[2]  #number of caches distibuted   
                if num_cache_equs<temp 
                        num_cache_equs=temp
                end 
                equs[varNum]=rhs
            end  
        elseif argI isa Expr && argI.head==:function
            push!(functionEx.args,argI)
        elseif @capture(argI, for counter_ in b_:niter_ loopbody__ end)
             specRHS=loopbody[1].args[2] 
             extractJacDepLoopDiscrete(b,niter,specRHS,jac,exacteJacExpr,symDict,dD ) 
             temp=(transformF(:($(specRHS),1))).args[2]
                if num_cache_equs<temp 
                    num_cache_equs=temp
                end 
               equs[:(($b,$niter))]=specRHS            
        elseif argI isa Expr && argI.head==:if   
           zcf=argI.args[1].args[2]
           ZCFCounter+=1
           extractZCJacDepNormal(ZCFCounter,zcf,ZCjac ,SZ ,dZ ) 
            if zcf.head==:ref  #if one_Var
                  push!(zcequs,(transformFSimplecase(:($(zcf)))))    
            else # if whole expre ops with many vars            
                temp=:($((transformF(:($(zcf),1))).args[2]))   #number of caches distibuted, given 1 place holder for ex.args[2] to be filled inside and returned
                if num_cache_equs<temp 
                      num_cache_equs=temp
                end 
                push!(zcequs,(zcf))      
            end     
            handleEvents(argI,eventequs,length(zcequs),evsArr) #extract events
        elseif argI.args[1] isa Symbol && argI.args[2] isa Number  # already changed in main and plugged in diff eqs
        elseif  argI.args[1] isa Symbol && argI.args[2] isa Expr && (argI.args[2].head==:call || argI.args[2].head==:ref) #already changed in main and plugged in diff eqs
        else# keep any other code written by user
            push!(otherCode.args,argI)
        end #end cases 
 
    end #end for #

    fname= prbName # problem name received as outside function name as written by user
    if odeExprs.args[1] isa Expr && odeExprs.args[1].args[2] isa Expr && odeExprs.args[1].args[2].head == :tuple#old usage: user has to enter problem info as a string in a tuple
        fname= Symbol(odeExprs.args[1].args[2].args[1])
    end

    if length(functionEx.args)>1
        #closurefunc=()->nothing
        closurefunc=@RuntimeGeneratedFunction(functionEx.args[2]) 
        diffEqfunction=createDiscEqFun(otherCode,equs,zcequs,eventequs,fname,closurefunc)# diff equations before this are stored in a dict:: now we have a giant function that holds all diff equations
        exactJacfunction=createExactJacFun(otherCode,exacteJacExpr,fname,closurefunc)    
    else
        closurefunc=0
        diffEqfunction=createDiscEqFun(otherCode,equs,zcequs,eventequs,fname,closurefunc)# diff equations before this are stored in a dict:: now we have a giant function that holds all diff equations
        exactJacfunction=createExactJacFun(otherCode,exacteJacExpr,fname,closurefunc)
    end
    exactJacfunctionF=@RuntimeGeneratedFunction(exactJacfunction)
    functioncodeF=@RuntimeGeneratedFunction(diffEqfunction)

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
    # mapFunF=@RuntimeGeneratedFunction(mapFun)
    myodeProblem = NLODEDiscProblemSpan(fname,Val(1),Val(T),Val(D),Val(Z),Val(num_cache_equs),initCond, collect(discrVars), jacVect ,ZCjac  ,functioncodeF, evsArr,SDVect,HZ,HD,SZVect,exactJacfunctionF,tspan,[closurefunc])
end

#old interface without tspan
function NLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{D},::Val{Z},initCond::Vector{Float64},du::Symbol)where {T,D,Z}
    discrVars=Vector{Float64}()
    tspan = (0.0,1.0)
    prbName=:_
    probspan=NLodeProblemFunc(odeExprs,Val(T),Val(D),Val(Z),initCond,du,tspan,discrVars,prbName)  
    myodeProblem = NLODEDiscProblem(probspan.prname,Val(1),Val(T),Val(D),Val(Z),probspan.cacheSize,probspan.initConditions, probspan.discreteVars, probspan.jac ,probspan.ZCjac  ,probspan.eqs, probspan.eventDependencies,probspan.SD,probspan.HZ,probspan.HD,probspan.SZ,probspan.exactJac, probspan.closureFuncs)
end



