
"""
    NLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{D},::Val{Z},initCond::Vector{Float64},du::Symbol,symDict::Dict{Symbol,Expr},tspan::Tuple{Float64, Float64},discrVars::Vector{EM},prbName::Symbol) where {T,D,Z,EM}

        continues building a discrete problem. 
    It receives an expression and useful info from the main interface. It calls the transform function from the taylorEquationConstruction.jl file to change the AST of all operations to personlized ones and update the needed cache size. It also construct via helper functions the Exact jacobian function, the jacobian dependecy (jac) and the state-derivative dependency (SD:opposite of jacobian), the state to zero-crossing dependency (SZ) and events to derivative and zero-crossing (HD and HZ) as vectors. Finally, it groups all differential equations and events in one function, and constructs a discrete problem from the qssProblemDefinition.jl file.
"""
function NLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{D},::Val{Z},initCond::Vector{Float64},du::Symbol,symDict::Dict{Symbol,Expr},tspan::Tuple{Float64, Float64},discrVars::Vector{EM},prbName::Symbol) where {T,D,Z,EM} # used EM to account for when problem contains if-statements whithout discrete vars
    if VERBOSE println("discrete nlodeprobfun  T D Z= $T $D $Z") end
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
    for argI in odeExprs.args
        if argI isa Expr &&  argI.head == :(=) && argI.args[1]== :discrete # old approach use the keyword discrete to define discrete variables
            discrVars = Vector{Float64}(argI.args[2].args)
         #only diff eqs: du[]= number/one ref/call  
        elseif argI isa Expr &&  argI.head == :(=)  && argI.args[1] isa Expr && argI.args[1].head == :ref && argI.args[1].args[1]==du
            y=argI.args[1];rhs=argI.args[2]
            varNum=y.args[2] # order of variable
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
            ##################################################################################################################
            #                                                      events     
            ##################################################################################################################                  
            # each 'if-statmets' has 2 events (arg[2]=posEv and arg[3]=NegEv) each pos or neg event has a function...later i can try one event for zc
            if length(argI.args)==2  #if user only wrote the positive evnt, here I added the negative event wich does nothing
                nothingexpr = quote nothing end # neg dummy event 
                push!(argI.args, nothingexpr)
                Base.remove_linenums!(argI.args[3])
            end
            #pos event
            newPosEventExprToFunc=changeExprToFirstValue(argI.args[2])  #change u[1] to u[1][0]  # pos ev can't be a symbol ...later maybe add check anyway
            push!(eventequs,newPosEventExprToFunc) 
            #neg eve
            if argI.args[3].args[1] isa Expr # argI.args[3] isa expr and is either an equation or :block symbol end ...so lets check .args[1]
                newNegEventExprToFunc=changeExprToFirstValue(argI.args[3])
                push!(eventequs,newNegEventExprToFunc) 
            else
                push!(eventequs,argI.args[3]) #symbol nothing
            end
            #after constructing the equations we move to dependencies: we need to change A[n] to An so that they become symbols
            posEvExp =  argI.args[2]
            negEvExp =  argI.args[3]
            indexPosEv = 2 * length(zcequs) - 1 # store events in order
            indexNegEv = 2 * length(zcequs)   
              #------------------pos Event--------------------#
            posEv_disArrLHS= Vector{Int}()  
            posEv_conArrLHS= Vector{Int}() 
            posEv_conArrRHS=Vector{Int}()    #to be used inside intgrator to updateOtherQs (intgrateState) before executing the event there is no discArrRHS because d is not changing overtime to be updated      
            for j = 1:length(posEvExp.args)  # j coressponds the number of statements under one posEvent
                if (posEvExp.args[j]  isa Expr &&  posEvExp.args[j].head == :(=)) 
                        poslhs=posEvExp.args[j].args[1];posrhs=posEvExp.args[j].args[2]
                    if (poslhs  isa Expr &&  poslhs.head == :ref && (poslhs.args[1]==:q || poslhs.args[1]==:d))    
                       if poslhs.args[1]==:q
                            push!(posEv_conArrLHS,poslhs.args[2])
                        else # lhs is a disc var 
                            push!(posEv_disArrLHS,poslhs.args[2])
                        end
                        postwalk(posrhs) do a   #
                            if a isa Expr && a.head == :ref && a.args[1]==:q# 
                                push!(posEv_conArrRHS,  (a.args[2]))  #                    
                            end
                            return a 
                        end
                    end
                end
            end
            #------------------neg Event--------------------#
            negEv_disArrLHS= Vector{Int}()#
            negEv_conArrLHS= Vector{Int}()# 
            negEv_conArrRHS=Vector{Int}()#to be used inside intgrator to updateOtherQs (intgrateState) before executing the event there is no discArrRHS because d is not changing overtime to be updated      
            if negEvExp.args[1] != :nothing
                for j = 1:length(negEvExp.args)  # j coressponds the number of statements under one negEvent
                    neglhs=negEvExp.args[j].args[1];negrhs=negEvExp.args[j].args[1]
                    if (neglhs  isa Expr &&  neglhs.head == :ref && (neglhs.args[1]==:q || neglhs.args[1]==:d))    
                        if neglhs.args[1]==:q
                            push!(negEv_conArrLHS,neglhs.args[2])
                        else # lhs is a disc var 
                            push!(negEv_disArrLHS,neglhs.args[2])
                        end
                        postwalk(negrhs) do a   #
                            if a isa Expr && a.head == :ref && a.args[1]==:q# 
                                push!(negEv_conArrRHS,  (a.args[2]))  #                    
                            end
                            return a 
                        end
                    end
                end
            end 
            structposEvent = EventDependencyStruct(indexPosEv, posEv_conArrLHS, posEv_disArrLHS,posEv_conArrRHS) # posEv_conArr is vect 
            push!(evsArr, structposEvent)
            structnegEvent = EventDependencyStruct(indexNegEv, negEv_conArrLHS, negEv_disArrLHS,negEv_conArrRHS)
            push!(evsArr, structnegEvent)
        end #end cases 
 
    end #end for #
    allEpxpr=Expr(:block)
    ##############diffEqua###############
    s="if i==0 return nothing\n"  # :i is the mute var
    for elmt in equs
        Base.remove_linenums!(elmt[1])
        Base.remove_linenums!(elmt[2])
        if elmt[1] isa Int
            s*="elseif i==$(elmt[1]) $(elmt[2]) ;return nothing\n"
        end
        if elmt[1] isa Expr
            s*="elseif $(elmt[1].args[1])<=i<=$(elmt[1].args[2]) $(elmt[2]) ;return nothing\n"
        end
    end
    s*=" end "
    myex1=Meta.parse(s)
    push!(allEpxpr.args,myex1)
     ##############ZCF###################
    if length(zcequs)>0
        s="if zc==1  $(zcequs[1]) ;return nothing"
        for i=2:length(zcequs)
            s*= " elseif zc==$i $(zcequs[i]) ;return nothing"
        end
        s*= " end "
        myex2=Meta.parse(s)
        push!(allEpxpr.args,myex2)
    end
     #############events#################
    if length(eventequs)>0
        s= "if ev==1  $(eventequs[1]) ;return nothing"
        for i=2:length(eventequs)
            s*= " elseif ev==$i $(eventequs[i]) ;return nothing"
        end
        s*= " end "
        myex3=Meta.parse(s)
        push!(allEpxpr.args,myex3)
    end
    fname= prbName # problem name received as outside function name as written by user
    if odeExprs.args[1] isa Expr && odeExprs.args[1].args[2] isa Expr && odeExprs.args[1].args[2].head == :tuple#user has to enter problem info in a tuple
        fname= odeExprs.args[1].args[2].args[1]
    end
    Base.remove_linenums!(allEpxpr)
    def=Dict{Symbol,Any}()
    def[:head] = :function
    def[:name] = fname  
    def[:args] = [:(i::Int),:(zc::Int),:(ev::Int),:(q::Vector{Taylor0}),:(d::Vector{Float64}), :(t::Taylor0),:(cache::Vector{Taylor0})]
    def[:body] = allEpxpr 
    functioncode=combinedef(def)
   # @show functioncode
    functioncodeF=@RuntimeGeneratedFunction(functioncode)
    jacVect=createJacVect(jac,Val(T))
    SDVect=createSDVect(jac,Val(T))
    dDVect =createdDvect(dD,Val(D))
    SZvect=createSZvect(SZ,Val(T))
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
    exactJacfunction=createExactJacFun(exacteJacExpr,fname)
    exactJacfunctionF=@RuntimeGeneratedFunction(exactJacfunction)
    if VERBOSE println("discrete problem created") end
    myodeProblem = NLODEDiscProblemSpan(fname,Val(1),Val(T),Val(D),Val(Z),Val(num_cache_equs),initCond, discrVars, jacVect ,ZCjac  ,functioncodeF, evsArr,SDVect,HZ,HD,SZvect,exactJacfunctionF,tspan)
end

#old interface without tspan
function NLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{D},::Val{Z},initCond::Vector{Float64},du::Symbol,symDict::Dict{Symbol,Expr})where {T,D,Z}
    discrVars=Vector{Float64}()
    tspan = (0.0,1.0)
    prbName=:_
    probspan=NLodeProblemFunc(odeExprs,Val(T),Val(D),Val(Z),initCond,du,symDict,tspan,discrVars,prbName)  
    myodeProblem = NLODEDiscProblem(probspan.prname,Val(1),Val(T),Val(D),Val(Z),probspan.cacheSize,probspan.initConditions, probspan.discreteVars, probspan.jac ,probspan.ZCjac  ,probspan.eqs, probspan.eventDependencies,probspan.SD,probspan.HZ,probspan.HD,probspan.SZ,probspan.exactJac)
end



