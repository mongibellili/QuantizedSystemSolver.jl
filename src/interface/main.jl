




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
function ODEProblem(f::Function, u::Vector{C}, tspan::Tuple{A,B}, p::Union{Vector{EM}, Tuple{Vararg{EM}}};jac_mode ::Symbol= :symbolic) where{EM,A<:Union{Float64, Int64},B<:Union{Float64, Int64},C<:Union{Float64, Int64}}
    bt = stacktrace()
    is_top_level=false
    for i in 2:min(5, length(bt))
        if String(bt[i].func) == "top-level scope"
            if String(bt[i-1].func) == "ODEProblem"
                is_top_level=true
            end
        end
    end
   odeString=""
    for m in methods(f).ms
        if length(m.sig.parameters) == 5   
            odeString = @code_string f(u, u, p, tspan)
            break
        elseif length(m.sig.parameters) == 4
            odeString = @code_string f(u, u, tspan)
            break
        end
        error("Function $(f) must accept exactly 3 arguments (du, u, t) or 4 arguments (du, u, p, t).") 
    end
    odeex = Meta.parse(odeString)
    Base.remove_linenums!(odeex)
    stateVarName = odeex.args[1].args[3] # Symbol(:u) or other
     if length(odeex.args[1].args)==4
        discrParamName = odeex.args[1].args[4] # Symbol(:p) or other
     else
        discrParamName = :p # if no parameters are given, we create a dummy symbol
     end
    problemName = odeex.args[1].args[1]
    odeBodyExprs = odeex.args[2] #function body
    probSize = length(u)
    discParamSize = length(p) # parameters may or may not be discrete.
 
    probInfo=problem_to_normalized_ir(odeBodyExprs, stateVarName, discrParamName) # combine build and normalize IR
    ir= probInfo.ir # get the ir from probInfo
    #@show ir
    zcSize = probInfo.numZC
    numHelperF=probInfo.helperFunSymSet
    du = odeex.args[1].args[2] # Symbol(:du)
    if !(u isa Vector{Float64}) u = Float64.(u) end
    tspan = Float64.(tspan)  # Convert each element in tuple

    mod = typeof(f).name.module # Get the module where the function is defined
    preProcessData=PreProcessData(du,tspan,problemName,mod,is_top_level,numHelperF)
    odeProblemFunc(ir, Val(probSize), Val(discParamSize), Val(zcSize), u,  p,preProcessData,jac_mode) # returns prob
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
function ODEProblem(f::Function, u::Vector{C}, tspan::Tuple{A,B};jac_mode ::Symbol= :symbolic) where{A<:Union{Float64, Int64},B<:Union{Float64, Int64},C<:Union{Float64, Int64}}
    p=Float64[]
    ODEProblem(f, u, tspan, p,;jac_mode =jac_mode)
end

struct probInfo #helper struct to return number of "if_statements" and a dictionary of symbols from prepareInfo
    numZC::Int
    helperFunSymSet::Int64
end
