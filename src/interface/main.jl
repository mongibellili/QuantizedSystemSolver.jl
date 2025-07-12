




"""
    ODEProblem(f::Function, u::Vector{Float64}, tspan::Tuple{A,B}, p::Union{Vector{EM}, Tuple{Vararg{EM}}};jac_mode ::Symbol= :symbolic) where{EM,A<:Union{Float64, Int64},B<:Union{Float64, Int64}}

Creates an ODE problem with the given function `f`, initial conditions `u`, parameters `p`, and time span `tspan`.

# Arguments
- `f::Function`: The function defining the ODE.
- `u::Vector{Float64}`: The initial conditions.
- `p::Vector{Float64}`: The parameters.
- `tspan::Tuple{Float64,Float64}`: The time span for the ODE.

# Returns
- An ODE problem.
"""
function ODEProblem(f::Function, u::Vector{Float64}, tspan::Tuple{A,B}, p::Union{Vector{EM}, Tuple{Vararg{EM}}};jac_mode ::Symbol= :symbolic) where{EM,A<:Union{Float64, Int64},B<:Union{Float64, Int64}}
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
    probSize = length(u)
    discSize = length(p) 
    if EM==Any
        p=Float64[]#if user enters an empty vector... not force user to enter a type and keep using the type float64
    end
    #ir=build_ir(odeExprs) #build IR
    #probInfo =normalize_ir(ir, stateVarName, discrParamName) # normalize IR
    probInfo=problem_to_normalized_ir(odeExprs, stateVarName, discrParamName) # combine build and normalize IR
    ir= probInfo.ir # get the ir from probInfo
    zcSize = probInfo.numZC
    numHelperF=probInfo.helperFunSymSet
    #symDict = probInfo.symDict
    du = odeex.args[1].args[2] # Symbol(:du)
    #u = Float64.(u)  # Convert each element to Float64
    tspan = Float64.(tspan)  # Convert each element in tuple
    p=Float64.(p)
    mod = typeof(f).name.module # Get the module where the function is defined
    preProcessData=PreProcessData(du,tspan,problemName,mod,is_top_level,numHelperF)
   # @show zcSize
    #odeProblemFunc(odeExprs, Val(probSize), Val(discSize), Val(zcSize), u,  p,preProcessData,jac_mode) # returns prob
    odeProblemFunc(ir, Val(probSize), Val(discSize), Val(zcSize), u,  p,preProcessData,jac_mode) # returns prob
end

"""
    ODEProblem(f::Function, u::Vector{Float64}, tspan::Tuple{A,B};jac_mode ::Symbol= :symbolic) where{A<:Union{Float64, Int64},B<:Union{Float64, Int64}}

Creates an ODE problem with the given function `f`, initial conditions `u`, and time span `tspan`.
# Arguments
- `f::Function`: The function defining the ODE.
- `u::Vector{Float64}`: The initial conditions.
- `tspan::Tuple{Float64,Float64}`: The time span for the ODE.
# Returns
- An ODE problem.
"""
function ODEProblem(f::Function, u::Vector{Float64}, tspan::Tuple{A,B};jac_mode ::Symbol= :symbolic) where{A<:Union{Float64, Int64},B<:Union{Float64, Int64}}
    p=Float64[]
    ODEProblem(f, u, tspan, p,;jac_mode =jac_mode)
end

struct probInfo #helper struct to return number of "if_statements" and a dictionary of symbols from prepareInfo
    numZC::Int
    helperFunSymSet::Int64
end
