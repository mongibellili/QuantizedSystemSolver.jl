"""
    ODEProblem(f::Function, u::Vector{Float64}, tspan::Tuple{A,B}, p::Union{Vector{EM}, Tuple{Vararg{EM}}}) where{EM,A<:Union{Float64, Int64},B<:Union{Float64, Int64}}

Creates an ODE problem with the given function `f`, initial conditions `u`, parameters `p`, and time span `tspan`.

# Arguments
- `f::Function`: The function defining the ODE.
- `u::Vector{Float64}`: The initial conditions.
- `p::Vector{Float64}`: The parameters.
- `tspan::Tuple{Float64,Float64}`: The time span for the ODE.

# Returns
- An ODE problem.
"""
function ODEProblem(f::Function, u::Vector{Float64}, tspan::Tuple{A,B}, p::Union{Vector{EM}, Tuple{Vararg{EM}}}) where{EM,A<:Union{Float64, Int64},B<:Union{Float64, Int64}}
   #=  uu=deepcopy(u)
    pp=deepcopy(p)
    f(uu, uu, pp, tspan) =#

    bt = stacktrace()
    is_top_level=false
    for i in 2:min(5, length(bt))
        if String(bt[i].func) == "top-level scope"
            if String(bt[i-1].func) == "ODEProblem"
                is_top_level=true
            end
        end
    end


    odeestring = @code_string f(u, u, p, tspan)
    odeex = Meta.parse(odeestring)
    Base.remove_linenums!(odeex)
    odeExprs = odeex.args[2] #function body
    stateVarName = odeex.args[1].args[3] # Symbol(:u) or other
    discrParamName = odeex.args[1].args[4] # Symbol(:u) or other
    problemName = odeex.args[1].args[1]
    probInfo = prepareInfo(odeExprs, stateVarName,discrParamName) # replace symbols and params, extract info about size, symbols, initconds
    probSize = length(u)
    discSize = length(p) 
    if EM==Any
        p=Float64[]#if user enters an empty vector... not force user to enter a type and keep using the type float64
    end
    zcSize = probInfo.numZC
    numHelperF=probInfo.helperFunSymSet
    #symDict = probInfo.symDict
    du = odeex.args[1].args[2] # Symbol(:du)
    #u = Float64.(u)  # Convert each element to Float64
    tspan = Float64.(tspan)  # Convert each element in tuple
    p=Float64.(p)
    mod = typeof(f).name.module # Get the module where the function is defined
    preProcessData=PreProcessData(du,tspan,problemName,mod,is_top_level,numHelperF)
    NLodeProblemFunc(odeExprs, Val(probSize), Val(discSize), Val(zcSize), u,  p,preProcessData) # returns prob
end

"""
    ODEProblem(f::Function, u::Vector{Float64}, tspan::Tuple{A,B}) where{A<:Union{Float64, Int64},B<:Union{Float64, Int64}}

Creates an ODE problem with the given function `f`, initial conditions `u`, and time span `tspan`.
# Arguments
- `f::Function`: The function defining the ODE.
- `u::Vector{Float64}`: The initial conditions.
- `tspan::Tuple{Float64,Float64}`: The time span for the ODE.
# Returns
- An ODE problem.
"""
function ODEProblem(f::Function, u::Vector{Float64}, tspan::Tuple{A,B}) where{A<:Union{Float64, Int64},B<:Union{Float64, Int64}}
    p=Float64[]
    ODEProblem(f, u, tspan, p)
end

struct probInfo #helper struct to return number of "if_statements" and a dictionary of symbols from prepareInfo
    numZC::Int
    helperFunSymSet::Int64
end
"""
    prepareInfo(odeExprs::Expr,stateVarName::Symbol,discrParamName::Symbol) 

Prepares information about the ODE problem by replacing symbols and parameters, and extracting information about size, symbols, and initial conditions.

# Arguments
- `odeExprs::Expr`: The expressions defining the ODE.
- `stateVarName::Symbol`: The Symbol representing the continuous variables.

# Returns
- A `probInfo` struct containing the number of zero-crossings (`numZC`) and a dictionary of symbols and expressions (`symDict`).
"""
function prepareInfo(odeExprs::Expr,stateVarName::Symbol,discrParamName::Symbol) # replace symbols and params , extract info about sizes,symbols,initconds
    param=Dict{Symbol,Union{Float64,Int64,Expr,Symbol}}()
    helperFunSymSet=Set{Symbol}()
    numZC=0
    for argI in odeExprs.args
        if argI isa Expr &&  argI.head == :(=) #&& @capture(argI, y_ = rhs_)  
            lhs=argI.args[1];rhs=argI.args[2]
            if lhs isa Symbol && rhs isa Number #params: fill dict of param
                #rhs=Float64(rhs)
                param[lhs]=Float64(rhs)
            elseif rhs isa Symbol && lhs isa Expr && lhs.head==:tuple #params: fill dict of param
                if rhs==discrParamName
                  rhs=:p
                elseif rhs==stateVarName
                    rhs=:q
                end
                for k=1: length(lhs.args)
                    param[lhs.args[k]]=:($rhs[$k])
                end
            elseif lhs isa Symbol && rhs isa Expr && (rhs.head==:call || rhs.head==:ref) #params=epression fill dict of param
                argI.args[2]=changeVarNames_params(rhs,stateVarName,discrParamName,:nothing,param,helperFunSymSet)
                param[lhs]=argI.args[2]
            elseif lhs isa Expr && lhs.head == :ref 
                if lhs.args[2] isa Expr && lhs.args[2].head==:call   # left hand side
                    lhs.args[2]=Int(eval(changeVarNames_params(lhs.args[2],stateVarName,discrParamName,:nothing,param,helperFunSymSet)))
                end
                if (rhs isa Expr && rhs.head !=:vect)#&& rhs.head==:call or ref # a diff equa not in a loop
                    argI.args[2]=changeVarNames_params(rhs,stateVarName,discrParamName,:nothing,param,helperFunSymSet)
                    if haskey(param, lhs.args[2])#Symbol in LHS is a parameter
                        lhs.args[2]=copy(param[lhs.args[2]]) 
                    end
                   
                elseif rhs isa Symbol
                    if haskey(param, rhs)#Symbol is a parameter 
                        argI.args[2]=copy(param[rhs]) 
                    end
                #elseif rhs 
                end
                #the following 3 commented lines are the beginning of the implementation of allowing du[expression]
               #=  
                if lhs.args[2] isa Expr && lhs.args[2].head==:call
                    lhs.args[2]=changeVarNames_params(lhs.args[2],stateVarName,discrParamName,:nothing,param,helperFunSymSet)
                end =#
           
            #= else
                argI.args[2]=changeVarNames_params(rhs,stateVarName,discrParamName,:nothing,param,helperFunSymSet) # multiple params on lhs and rhs contains symbols =#
            end
        elseif argI isa Expr && argI.head==:function 
            param[argI.args[1].args[1]]=:f_  #add counter if to allow more than one function
        elseif @capture(argI, for muteVar_ in b_:niter_ loopbody__ end)
            argI.args[2]=changeVarNames_params(loopbody[1],stateVarName,discrParamName,muteVar,param,helperFunSymSet) 
            #the following 5 commented lines are the beginning of the implementation of allowing many du inside the loop
           #=  tempVect=:()
            for lpbody in loopbody
                push!(tempVect.args,changeVarNames_params(lpbody,stateVarName,discrParamName,muteVar,param,helperFunSymSet))
            end
            argI.args[2]=tempVect =#
            if !(argI.args[1].args[2].args[2] isa Int64)#if b is not a number
                argI.args[1].args[2].args[2]=Int(eval(changeVarNames_params(argI.args[1].args[2].args[2],stateVarName,discrParamName,:nothing,param,helperFunSymSet)))#
            end
            if !(argI.args[1].args[2].args[3] isa Int64)#if niter is not a number
                argI.args[1].args[2].args[3]=Int(eval(changeVarNames_params(argI.args[1].args[2].args[3],stateVarName,discrParamName,:nothing,param,helperFunSymSet)))#
            end
        elseif argI isa Expr && argI.head==:if
            numZC+=1
            (length(argI.args)!=3 && length(argI.args)!=2) && error("use format if A>B C else D or if A>B C")
            # !(argI.args[1] isa Expr && argI.args[1].head==:call && argI.args[1].args[1]==:> && (argI.args[1].args[3]==0||argI.args[1].args[3]==0.0)) && error("use the format 'if a>0: change if a>b to if a-b>0")
            #   !(argI.args[1].args[2] isa Expr) && error("LHS of >  must be be an expression!")
            zcf_LHS=argI.args[1].args[2]
            zcf_RHS=argI.args[1].args[3]
            if argI.args[1].args[1] in [:>,:>=]
            elseif argI.args[1].args[1] in [:<,:<=]
                        argI.args[1].args[1]=:>
                zcf_RHS=argI.args[1].args[2]
                zcf_LHS=argI.args[1].args[3]
                
            else
                error("use :>,:>= or :<,:<= in if statement")
            end
               
            
            zcf=:()
            zcf.head=:call
            zcf.args=[:-,zcf_LHS,zcf_RHS]
            #argI.args[1].args[2]=zcf
            argI.args[1].args[2]=changeVarNames_params(zcf,stateVarName,discrParamName,:nothing,param,helperFunSymSet)#zcf
            # name changes have to be per block
            to_delete = Int[] # delete the lines that are no longer needed in the body of the if statement
            for (i, argii) in enumerate(argI.args[2].args)  # body of the if statement
                if argii isa Expr &&  argii.head in [:(=),:+=, :-=, :*=, :/=] 
                    bodyif_lhs=argii.args[1];bodyif_rhs=argii.args[2]
                    argii.args[1]=changeVarNames_params(bodyif_lhs,stateVarName,discrParamName,:nothing,param,helperFunSymSet) # first change the lhs of the if statement: symbol to u[1] or p[1] for example
                    if argii.args[1] isa Symbol && bodyif_rhs isa Expr && (bodyif_rhs.head==:call || bodyif_rhs.head==:ref) #params=epression fill dict of param: #if there is a param that holds an expression, register it in the param so that it is changed inside the event so that dependency is tracked 
                        argii.args[2]=changeVarNames_params(bodyif_rhs,stateVarName,discrParamName,:nothing,param,helperFunSymSet) 
                        param[bodyif_lhs]=argii.args[2]#copy(argii.args[2])
                        push!(to_delete, i)
                    else
                        argii.args[2]=changeVarNames_params(bodyif_rhs,stateVarName,discrParamName,:nothing,param,helperFunSymSet) #posEv 
                    end
                else #
                    argI.args[2].args[i]=changeVarNames_params(argI.args[2].args[i],stateVarName,discrParamName,:nothing,param,helperFunSymSet) 
                    
                end
            end

            for i in reverse(to_delete)
                splice!(argI.args[2].args, i)
            end
            
            if length(argI.args)==3 #user used if a b else c
                to_delete = Int[]
                for (i, argii) in enumerate(argI.args[3].args)  # body of the if statement
                    if argii isa Expr &&  argii.head in [:(=),:+=, :-=, :*=, :/=] 
                        bodyif_lhs=argii.args[1];bodyif_rhs=argii.args[2]
                        argii.args[1]=changeVarNames_params(bodyif_lhs,stateVarName,discrParamName,:nothing,param,helperFunSymSet) # first change the lhs of the if statement: symbol to u[1] or p[1] for example
                        if argii.args[1] isa Symbol && bodyif_rhs isa Expr && (bodyif_rhs.head==:call || bodyif_rhs.head==:ref) #params=epression fill dict of param: #if there is a param that holds an expression, register it in the param so that it is changed inside the event so that dependency is tracked 
                            argii.args[2]=changeVarNames_params(bodyif_rhs,stateVarName,discrParamName,:nothing,param,helperFunSymSet) 
                            param[bodyif_lhs]=argii.args[2]
                            push!(to_delete, i)
                        else
                            argii.args[2]=changeVarNames_params(bodyif_rhs,stateVarName,discrParamName,:nothing,param,helperFunSymSet) #posEv 
                        end
                    else #
                        argI.args[3].args[i]=changeVarNames_params(argI.args[3].args[i],stateVarName,discrParamName,:nothing,param,helperFunSymSet) 
                    end
                end
              
                for i in reverse(to_delete)
                    splice!(argI.args[3].args, i)
                end

            end
        #else
            #argI=changeVarNames_params(argI,stateVarName,discrParamName,:nothing,param,helperFunSymSet)
        end#end cases of argI
    end#end for argI in args

    numHelperF=length(helperFunSymSet)
    p=probInfo(numZC,numHelperF)
    
end#end function



