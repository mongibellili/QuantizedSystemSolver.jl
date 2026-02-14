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
function IR(f::Function, u::Vector{C}, tspan::Tuple{A,B};inline_mode::InlineMode=FULL, p::Union{Vector{EM}, Tuple{Vararg{EM}}}=Float64[]) where {EM,A<:Union{Float64, Int64},B<:Union{Float64, Int64},C<:Union{Float64, Int64}}
    odeString=""
    for m in methods(f).ms
        if length(m.sig.parameters) == 5   
            odeString = @code_string f(u, u, p, tspan)
            break
        elseif length(m.sig.parameters) == 4
            odeString = @code_string f(u, u, tspan)
            break
        end
    end
    odeex = Meta.parse(odeString)
    Base.remove_linenums!(odeex)

    stateVarName = odeex.args[1].args[3]
    discrParamName = length(odeex.args[1].args)==4 ? odeex.args[1].args[4] : :p
    problemName = odeex.args[1].args[1]
    du = odeex.args[1].args[2]
    odeBodyExprs = odeex.args[2]

    probInfo = problem_to_normalized_ir(odeBodyExprs, stateVarName, discrParamName,inline_mode)

    return (ir = probInfo.ir,
            numZC = probInfo.numZC,
            helperFunSymSet = probInfo.helperFunSymSet,
            du = du,
            problemName = problemName)
end


function ODEProblem(f::Function, u::Vector{C}, tspan::Tuple{A,B}, p::Union{Vector{EM}, Tuple{Vararg{EM}}};inline_mode::InlineMode=FULL, jac_mode::Symbol=:symbolic,verbose::Bool=false) where {EM,A<:Union{Float64, Int64},B<:Union{Float64, Int64},C<:Union{Float64, Int64}}
    
    bt = stacktrace()
    is_top_level=false
    for i in 2:min(5, length(bt))
        if String(bt[i].func) == "top-level scope"
            if String(bt[i-1].func) == "ODEProblem"
                is_top_level=true
            end
        end
    end

    irInfo = IR(f, u, tspan, p=p, inline_mode=inline_mode) # helper function to get the IR and other info. separation is also used to allow to call this from user space (i.e IR is public API for debugging)

    ir = irInfo.ir
    verbose && @show ir
    zcSize = irInfo.numZC
    numHelperF = irInfo.helperFunSymSet
    du = irInfo.du
    problemName = irInfo.problemName
 

    if !(u isa Vector{Float64}) u = Float64.(u) end
    tspan = Float64.(tspan)
    probSize = length(u)
    discParamSize = length(p)
    mod = typeof(f).name.module

    preProcessData = PreProcessData(du, tspan, problemName, mod, is_top_level, numHelperF)

    return odeProblemFunc(ir, Val(probSize), Val(discParamSize), Val(zcSize), u, p, preProcessData, jac_mode)
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
function ODEProblem(f::Function, u::Vector{C}, tspan::Tuple{A,B};inline_mode::InlineMode=FULL,jac_mode ::Symbol= :symbolic,verbose::Bool=false) where{A<:Union{Float64, Int64},B<:Union{Float64, Int64},C<:Union{Float64, Int64}}
    p=Float64[]
    ODEProblem(f, u, tspan, p;inline_mode=inline_mode,jac_mode =jac_mode, verbose=verbose)
end

struct probInfo #helper struct to return number of "if_statements" and a dictionary of symbols from prepareInfo
    numZC::Int
    helperFunSymSet::Int64
end



# dummy macro to be used in IR generation. this is needed to allow users to inline or no_inline some expressions. the custom parsing module will act on the existence of these macros to decide whether to inline or not the expressions.
macro _inline(expr)
end
macro _noinline(expr)
end