

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
    ODEContProblem{F,PRTYPE,T,D,Z,CS}
A struct that holds the continuous problem with tspan. It has the following fields:  
  - `prname`: The name of the problem  
  - `prtype`: The type of the problem  
  - `a`: The size of the problem   
  - `c`: The number of discrete vars 
  - `b`: The number of zero crossing functions 
  - `cacheSize`: The size of the cache  
  - `initConditions`: The initial conditions of the problem  
  - `discreteVars` # to match the differentialEqation.jl interface that wants the parameter p to be part of the problem
  - `eqs`: The function that holds all the ODEs  
  - `jac`: The Jacobian dependency  
  - `SD`: The state derivative dependency  
  - `exactJac`: The exact Jacobian function  
  - `tspan::Tuple{Float64, Float64}`:  This field variable did not exist in the original ODEContProblem as this simulation time should part of the problem. However, to match the differentialEqation.jl interface, the tspan is added to the definition of the problem.
  - `closureFuncs::Vector{F}` # function that holds closure function inside system defined by user
"""
struct ODEContProblem{F,PRTYPE,T,D,Z,CS}<: ODEProblemData{F,PRTYPE,T,D,Z,CS} 
  prname::Symbol
  prtype::Val{PRTYPE}
  a::Val{T}
  c::Val{D}
  b::Val{Z}
  cacheSize::Val{CS}
  initConditions::Vector{Float64}  
  discreteVars::Vector{Float64}  
  jac::Vector{Vector{Int}}#Jacobian dependency..I have a der and I want to know which vars affect it...opposite of SD

  eqs::Function#function that holds all ODEs

  SD::Vector{Vector{Int}}#  I have a var and I want the der that are affected by it

  exactJac::Function #used only in the implicit integration: linear approximation
  tspan::Tuple{Float64, Float64}
  closureFuncs::Vector{F} # 
end



"""
    ODEDiscProblem{F,PRTYPE,T,D,Z,CS}
A struct that holds the Problem of a system of ODEs with a set of events with tspan. It has the following fields:  
  -`prname::Symbol  `
  -`prtype::Val{PRTYPE} ` 
  -`a::Val{T}  `
  -`c::Val{D}  `
  -`b::Val{Z} ` 
  -`cacheSize::Val{CS}  `
  -`initConditions::Vector{Float64} `  
  -`discreteVars::Vector{Float64} `   
  -`jac::Vector{Vector{Int}}`#Jacobian dependency..I have a der and I want to know which vars affect it...opposite of SD  
  -`ZCjac::Vector{Vector{Int}}` # to update other Qs before checking ZCfunction  
  -`eqs::Function`#function that holds all ODEs  
  -`eventDependencies::Vector{EventDependencyStruct} ` 
  -`SD::Vector{Vector{Int}}`#  I have a var and I want the der that are affected by it  
  -`HZ::Vector{Vector{Int}}`#  an ev occured and I want the ZC that are affected by it  
  -`HD::Vector{Vector{Int}}`#  an ev occured and I want the der that are affected by it  
  -`SZ::Vector{Vector{Int}}`#  I have a var and I want the ZC that are affected by it  
  -`exactJac::Function `#used only in the implicit integration: linear approximation  
  -`tspan::Tuple{Float64, Float64}`# This field variable did not exist in the original ODEDiscProblem as this simulation time should part of the problem. However, to match the differentialEqation.jl interface, the tspan is added to the definition of the problem.
  -`closureFuncs::Vector{F}` # function that holds closure function inside system defined by user
"""
struct ODEDiscProblem{F,PRTYPE,T,D,Z,CS}<: ODEProblemData{F,PRTYPE,T,D,Z,CS} 
  prname::Symbol
  prtype::Val{PRTYPE}
  a::Val{T}
  c::Val{D}
  b::Val{Z}
  cacheSize::Val{CS}
  initConditions::Vector{Float64}  
  discreteVars::Vector{Float64}  
  jac::Vector{Vector{Int}}#Jacobian dependency..I have a der and I want to know which vars affect it...opposite of SD
  ZCjac::Vector{Vector{Int}} # to update other Qs before checking ZCfunction
  eqs::Function#function that holds all ODEs
  eventDependencies::Vector{EventDependencyStruct}# 
  SD::Vector{Vector{Int}}#  I have a var and I want the der that are affected by it
  HZ::Vector{Vector{Int}}#  an ev occured and I want the ZC that are affected by it
  HD::Vector{Vector{Int}}#  an ev occured and I want the der that are affected by it
  SZ::Vector{Vector{Int}}#  I have a var and I want the ZC that are affected by it
  exactJac::Function #used only in the implicit integration: linear approximation
  tspan::Tuple{Float64, Float64}
  closureFuncs::Vector{F} # 
end


