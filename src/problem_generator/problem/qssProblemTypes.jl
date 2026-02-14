# with any equation, carry its helper assignments 
struct ScopedEquation
    helperAssignments::Vector{AbstractODEStatement}
    eqs_RHS::Union{Number,Symbol,Expr}  # later change to expr only because eqs_RHS hold rhs (expr) after transformation
end
"""
    EventDependencyStruct
 A struct that holds the event dependency information. It has the following fields:
  - `id::Int:` the id of the event
  - `evCont::Vector{Int}:` the index tracking used for HD & HZ. Also it is used to update q,quantum,recomputeNext when x is modified in an event
  - `evDisc::Vector{Int}:` the index tracking used for HD & HZ.
  - `evContRHS::Vector{Int}:` the index tracking used to update other Qs before executing the event

"""
struct EventDependencyStruct
  id::Int
  evCont::Vector{Int} #LHS: index tracking used for HD & HZ. Also it is used to update q,quantum,recomputeNext when x is modified in an event
  evDisc::Vector{Int} #LHS: index tracking used for HD & HZ.
  evContRHS::Vector{Int} #Here we look at the RHS: index tracking used to update other Qs before executing the event
end
# A structure to hold data required for preprocessing steps.
struct PreProcessData 
  du::Symbol
  tspan::Tuple{Float64, Float64}
  prbName::Symbol
  mod::Module
  is_top_level::Bool
  numHelperFunCalls::Int
end




"""
    ODEDiscProblem{JACMODE,T,D,Z,CS,F,JAC,CLS}
A struct that holds the Problem of a system of ODEs with a set of events with tspan. It has the following fields:  
  -`prname::Symbol  `
  -`prtype::Val{JACMODE} ` 
  -`num_cont_vars::Val{T}  `
  -`num_discr_vars::Val{D}  `
  -`num_zero_cross_func::Val{Z} ` 
  -`cacheSize::Val{CS}  `
  -`initConditions::Vector{Float64} `  
  -`discreteVars::Vector{Float64} `   
  -`jac::Vector{Vector{Int}}`#Jacobian dependency..I have a der and I want to know which vars affect it...opposite of SD  
  -`ZCjac::Vector{Vector{Int}}` # to update other Qs before checking ZCfunction  
  -`equations::Function`#function that holds all ODEs  
  -`eventDependencies::Vector{EventDependencyStruct} ` 
  -`SD::Vector{Vector{Int}}`#  I have a var and I want the der that are affected by it  
  -`HZ::Vector{Vector{Int}}`#  an ev occured and I want the ZC that are affected by it  
  -`HD::Vector{Vector{Int}}`#  an ev occured and I want the der that are affected by it  
  -`SZ::Vector{Vector{Int}}`#  I have a var and I want the ZC that are affected by it  
  -`exactJac::Function `#used only in the implicit integration: linear approximation  
  -`tspan::Tuple{Float64, Float64}`# This field variable did not exist in the original ODEDiscProblem as this simulation time should part of the problem. However, to match the differentialEqation.jl interface, the tspan is added to the definition of the problem.
  -`closureFuncs::Vector{F}` # function that holds closure function inside system defined by user
"""
mutable  struct ODEDiscProblem{JACMODE,T,D,Z,CS,F,JAC,CLS}<: ODEProblemData{JACMODE,T,D,Z,CS,F,JAC,CLS} 
  prname::Symbol
  prtype::Val{JACMODE}
  num_cont_vars::Val{T}
  num_discr_vars::Val{D}
  num_zero_cross_func::Val{Z} 
  cacheSize::Val{CS} 
  initConditions::Vector{Float64}  
  discreteVars#::Vector{Any}  
  jac::Vector{Vector{Int}}#Jacobian dependency..I have a der and I want to know which vars affect it...opposite of SD
  ZCjac::Vector{Vector{Int}} # to update other Qs before checking ZCfunction
  equations::F#function that holds all ODEs
  eventDependencies::Vector{EventDependencyStruct}# 
  SD::Vector{Vector{Int}}#  I have a var and I want the der that are affected by it
  HZ::Vector{Vector{Int}}#  an ev occured and I want the ZC that are affected by it
  HD::Vector{Vector{Int}}#  an ev occured and I want the der that are affected by it
  SZ::Vector{Vector{Int}}#  I have a var and I want the ZC that are affected by it
  exactJac::JAC #used only in the implicit integration: linear approximation
  tspan::Tuple{Float64, Float64}
  closureFuncs::Vector{CLS} # 
end


