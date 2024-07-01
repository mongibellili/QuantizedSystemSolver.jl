# Functions



```@docs
 solve(prob::NLODEProblem{PRTYPE,T,Z,D,CS},tspan::Tuple{Float64, Float64};sparsity::Val{Sparsity}=Val(false)::Float64,saveat=1e-9::Float64::Float64,abstol=1e-4::Float64,reltol=1e-3::Float64,maxErr=Inf::Float64) where {PRTYPE,T,Z,D,CS,Sparsity}    
```