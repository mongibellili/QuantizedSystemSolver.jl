"""
    ODEProblem(f::Function, u::Vector{Float64}, tspan::Tuple{Float64,Float64}, p::Vector{EM}) where {EM}

Creates an ODE problem with the given function `f`, initial conditions `u`, parameters `p`, and time span `tspan`.

# Arguments
- `f::Function`: The function defining the ODE.
- `u::Vector{Float64}`: The initial conditions.
- `p::Vector{Float64}`: The parameters.
- `tspan::Tuple{Float64,Float64}`: The time span for the ODE.

# Returns
- An ODE problem.
"""
function ODEProblem(f::Function, u::Vector{Float64}, tspan::Tuple{Float64,Float64}, p::Vector{EM}) where {EM}# T used because for empty p=[] type is any
    odeestring = @code_string f(u, u, p, tspan)
    odeex = Meta.parse(odeestring)
    Base.remove_linenums!(odeex)
    if VERBOSE println("starting prob parsing...") end
    odeExprs = odeex.args[2]
    contVarSymbole = odeex.args[1].args[3] # Symbol(:u) or other
    problemName = odeex.args[1].args[1]
    probInfo = prepareInfo(odeExprs, contVarSymbole) # replace symbols and params, extract info about size, symbols, initconds
    probSize = length(u)
    discSize = length(p) 
    if EM==Any
        p=[0.0] #empty vector of Float64 to not force user to enter a type and keep using the fields of type float64
    end
    zcSize = probInfo.numZC
    symDict = probInfo.symDict
    du = odeex.args[1].args[2] # Symbol(:du)
    NLodeProblemFunc(odeExprs, Val(probSize), Val(discSize), Val(zcSize), u, du, symDict, tspan, p, problemName) # returns prob
end

"""
    ODEProblem(f::Function, u::Vector{Float64}, tspan::Tuple{Float64,Float64})

Creates an ODE problem with the given function `f`, initial conditions `u`, and time span `tspan`.

# Arguments
- `f::Function`: The function defining the ODE.
- `u::Vector{Float64}`: The initial conditions.
- `tspan::Tuple{Float64,Float64}`: The time span for the ODE.

# Returns
- An ODE problem.
"""
function ODEProblem(f::Function, u::Vector{Float64}, tspan::Tuple{Float64,Float64})
    odeestring = @code_string f(u, u,u, tspan) 
    odeex = Meta.parse(odeestring)
    #odeex = @code_expr f(u, u, tspan) 
    #odeestring =code_string(f, (Function,Expr))
    
    Base.remove_linenums!(odeex)
    if VERBOSE println("starting prob parsing...") end
    odeExprs = odeex.args[2]
    contVarSymbole = odeex.args[1].args[3] # Symbol(:u) or other
    problemName = odeex.args[1].args[1]
    probInfo = prepareInfo(odeExprs, contVarSymbole) # replace symbols and params, extract info about size, symbols, initconds
    probSize = length(u)
    discSize = 0
    zcSize = probInfo.numZC
    symDict = probInfo.symDict
    du = odeex.args[1].args[2] # Symbol(:du)
    p=[0.0]
    NLodeProblemFunc(odeExprs, Val(probSize), Val(discSize), Val(zcSize), u, du, symDict, tspan,p, problemName) # returns prob
end

struct probInfo #helper struct to return stuff from prepareInfo
    numZC::Int
    symDict::Dict{Symbol,Expr}
end
"""
    prepareInfo(x::Expr,stateVarName::Symbol) 

Prepares information about the ODE problem by replacing symbols and parameters, and extracting information about size, symbols, and initial conditions.

# Arguments
- `x::Expr`: The expressions defining the ODE.
- `stateVarName::Symbol`: The symbol representing the continuous variables.

# Returns
- A `probInfo` struct containing the number of zero-crossings (`numZC`) and a dictionary of symbols and expressions (`symDict`).
"""
function prepareInfo(x::Expr,stateVarName::Symbol) # replace symbols and params , extract info about sizes,symbols,initconds
    param=Dict{Symbol,Union{Float64,Expr}}()
    symDict=Dict{Symbol,Expr}()
    numZC=0
    for argI in x.args
        if argI isa Expr &&  argI.head == :(=) #&& @capture(argI, y_ = rhs_) 
            y=argI.args[1];rhs=argI.args[2]
            if y isa Symbol && rhs isa Number #params: fill dict of param
                param[y]=rhs
            elseif y isa Symbol && rhs isa Expr && (rhs.head==:call || rhs.head==:ref) #params=epression fill dict of param
                    argI.args[2]=changeVarNames_params(rhs,stateVarName,:nothing,param,symDict)
                    param[y]=argI.args[2]
            elseif y isa Expr && y.head == :ref && (rhs isa Expr && rhs.head !=:vect)#&& rhs.head==:call or ref # a diff equa not in a loop
                argI.args[2]=changeVarNames_params(rhs,stateVarName,:nothing,param,symDict)
            elseif y isa Expr && y.head == :ref && rhs isa Symbol
                if haskey(param, rhs)#symbol is a parameter
                    argI.args[2]=copy(param[rhs]) 
                end
            end
        elseif @capture(argI, for var_ in b_:niter_ loopbody__ end)
            argI.args[2]=changeVarNames_params(loopbody[1],stateVarName,var,param,symDict)
        elseif argI isa Expr && argI.head==:if
            numZC+=1
            (length(argI.args)!=3 && length(argI.args)!=2) && error("use format if A>0 B else C or if A>0 B")
            !(argI.args[1] isa Expr && argI.args[1].head==:call && argI.args[1].args[1]==:> && (argI.args[1].args[3]==0||argI.args[1].args[3]==0.0)) && error("use the format 'if a>0: change if a>b to if a-b>0")
           #   !(argI.args[1].args[2] isa Expr) && error("LHS of >  must be be an expression!")
              argI.args[1].args[2]=changeVarNames_params(argI.args[1].args[2],stateVarName,:nothing,param)#zcf
              # name changes have to be per block
              if length(argI.args)==2 #user used if a b
                argI.args[2]=changeVarNames_params(argI.args[2],stateVarName,:nothing,param) #posEv
              elseif length(argI.args)==3 #user used if a b else c
                argI.args[2]=changeVarNames_params(argI.args[2],stateVarName,:nothing,param)#posEv
                argI.args[3]=changeVarNames_params(argI.args[3],stateVarName,:nothing,param)#negEv
              end
        end#end cases of argI
    end#end for argI in args
    p=probInfo(numZC,symDict)
end#end function
