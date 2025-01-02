"""
    NLODEContProblem{F,PRTYPE,T,D,Z,CS}
A struct that holds the continuous problem. It has the following fields:  
  - `prname`: The name of the problem  
  - `prtype`: The type of the problem  
  - `a`: The size of the problem  
  - `c`: The number of discrete events  
  - `b`: The number of zero crossing functions  
  - `cacheSize`: The size of the cache  
  - `initConditions`: The initial conditions of the problem  
  - `discreteVars` # to match the differentialEqation.jl interface that wants the parameter p to be part of the problem
  - `eqs`: The function that holds all the ODEs  
  - `jac`: The Jacobian dependency  
  - `SD`: The state derivative dependency  
  - `exactJac`: The exact Jacobian function  
  - `closureFuncs::Vector{F}` # function that holds closure function inside system defined by user
"""
struct NLODEContProblem{F,PRTYPE,T,D,Z,CS}<: NLODEProblem{F,PRTYPE,T,D,Z,CS} 
  prname::Symbol # problem name used to distinguish printed results
  prtype::Val{PRTYPE} # problem type: not used but created in case in the future we want to handle problems differently
  a::Val{T} #problem size based on number of vars: T is used not a: 'a' is a mute var
  c::Val{D} #number of discrete events=2*ZCF: Z is used not c: 'c' is a mute var
  b::Val{Z} #number of Zero crossing functions (ZCF) based on number of 'if statements': Z is used not b: 'b' is a mute var
  cacheSize::Val{CS}# CS= cache size is used  : 'cacheSize' is a mute var
  initConditions::Vector{Float64}  # 
  discreteVars::Vector{Float64} # to match the differentialEqation.jl interface that wants the parameter p to be part of the problem
  eqs::Function#function that holds all ODEs
  jac::Vector{Vector{Int}}#Jacobian dependency..I have a der and I want to know which vars affect it...opposite of SD...is a vect for direct method (later @resumable..closure..for saved method)
  SD::Vector{Vector{Int}}#  I have a var and I want the der that are affected by it
  exactJac::Function  # used only in the implicit intgration
  closureFuncs::Vector{F} # function that holds closure function inside system defined by user
end
"""
    NLODEContProblemSpan{F,PRTYPE,T,D,Z,CS}
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
  - `tspan::Tuple{Float64, Float64}`:  This field variable did not exist in the original NLODEContProblem as this simulation time should part of the problem. However, to match the differentialEqation.jl interface, the tspan is added to the definition of the problem.
  - `closureFuncs::Vector{F}` # function that holds closure function inside system defined by user
"""
struct NLODEContProblemSpan{F,PRTYPE,T,D,Z,CS}<: NLODEProblem{F,PRTYPE,T,D,Z,CS} 
  prname::Symbol # problem name used to distinguish printed results
  prtype::Val{PRTYPE} # problem type: not used but created in case in the future we want to handle problems differently
  a::Val{T} #problem size based on number of vars: T is used not a: 'a' is a mute var
  c::Val{D} #number of discrete events=2*ZCF: Z is used not c: 'c' is a mute var
  b::Val{Z} #number of Zero crossing functions (ZCF) based on number of 'if statements': Z is used not b: 'b' is a mute var
  cacheSize::Val{CS}# CS= cache size is used  : 'cacheSize' is a mute var
  initConditions::Vector{Float64}  # 
  discreteVars::Vector{Float64} # to match the differentialEqation.jl interface that wants the parameter p to be part of the problem
  eqs::Function#function that holds all ODEs
  jac::Vector{Vector{Int}}#Jacobian dependency..I have a der and I want to know which vars affect it...opposite of SD...is a vect for direct method (later @resumable..closure..for saved method)
  SD::Vector{Vector{Int}}#  I have a var and I want the der that are affected by it
  exactJac::Function  # used only in the implicit intgration
  tspan::Tuple{Float64, Float64}
  closureFuncs::Vector{F} # 
end
struct NLODEContAdvProblemSpan{F,PRTYPE,T,D,Z,CS}<: NLODEProblem{F,PRTYPE,T,D,Z,CS} 
  prname::Symbol # problem name used to distinguish printed results
  prtype::Val{PRTYPE} # problem type: not used but created in case in the future we want to handle problems differently
  a::Val{T} #problem size based on number of vars: T is used not a: 'a' is a mute var
  c::Val{D} #number of discrete events=2*ZCF: Z is used not c: 'c' is a mute var
  b::Val{Z} #number of Zero crossing functions (ZCF) based on number of 'if statements': Z is used not b: 'b' is a mute var
  cacheSize::Val{CS}# CS= cache size is used  : 'cacheSize' is a mute var
  initConditions::Vector{Float64}  # 
  discreteVars::Vector{Float64} # to match the differentialEqation.jl interface that wants the parameter p to be part of the problem
  eqs::Function#function that holds all ODEs
  jac::Vector{Vector{Int}}#Jacobian dependency..I have a der and I want to know which vars affect it...opposite of SD...is a vect for direct method (later @resumable..closure..for saved method)
  SD::Vector{Vector{Int}}#  I have a var and I want the der that are affected by it
  exactJac::Function  # used only in the implicit intgration
  tspan::Tuple{Float64, Float64}
  closureFuncs::Vector{F} # 
end

"""
    NLODEDiscProblem{F,PRTYPE,T,D,Z,CS}
A struct that holds the Problem of a system of ODEs with a set of events. It has the following fields:  
  - `prname::Symbol`  
  - `prtype::Val{PRTYPE}`  
  - `a::Val{T} ` 
  - `c::Val{D}  `
  - `b::Val{Z}`  
  - `cacheSize::Val{CS}  `
  - `initConditions::Vector{Float64} `  
  - `discreteVars::Vector{Float64}`    
  - `jac::Vector{Vector{Int}}`#Jacobian dependency..I have a der and I want to know which vars affect it...opposite of SD  
  - `ZCjac::Vector{Vector{Int}}` # to update other Qs before checking ZCfunction  
  - `eqs::Function`#function that holds all ODEs  
  - `eventDependencies::Vector{EventDependencyStruct}`  
  - `SD::Vector{Vector{Int}}`#  I have a var and I want the der that are affected by it  
  - `HZ::Vector{Vector{Int}}`#  an ev occured and I want the ZC that are affected by it  
  - `HD::Vector{Vector{Int}}`#  an ev occured and I want the der that are affected by it  
  - `SZ::Vector{Vector{Int}}`#  I have a var and I want the ZC that are affected by it  
  - `exactJac::Function` #used only in the implicit integration: linear approximation  
  - `closureFuncs::Vector{F}` # function that holds closure function inside system defined by user
"""
struct NLODEDiscProblem{F,PRTYPE,T,D,Z,CS}<: NLODEProblem{F,PRTYPE,T,D,Z,CS} 
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
  closureFuncs::Vector{F} # where {F<:Function}#function that holds all ODEs 
end
"""
    NLODEDiscProblemSpan{F,PRTYPE,T,D,Z,CS}
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
  -`tspan::Tuple{Float64, Float64}`# This field variable did not exist in the original NLODEDiscProblem as this simulation time should part of the problem. However, to match the differentialEqation.jl interface, the tspan is added to the definition of the problem.
  -`closureFuncs::Vector{F}` # function that holds closure function inside system defined by user
"""
struct NLODEDiscProblemSpan{F,PRTYPE,T,D,Z,CS}<: NLODEProblem{F,PRTYPE,T,D,Z,CS} 
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


