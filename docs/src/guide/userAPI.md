# User API


## Problem definition

```@docs
ODEProblem(f::Function, u::Vector{Float64}, tspan::Tuple{Float64,Float64}, p::Vector{EM}) where {EM}
```
```@docs
ODEProblem(f::Function, u::Vector{Float64}, tspan::Tuple{Float64,Float64})
```

```@docs 
NLodeProblem(odeExprs) 
```


## The solve function:

```@docs
solve(prob::NLODEProblem{PRTYPE,T,D,Z,CS},al::QSSAlgorithm{SolverType, OrderType};sparsity::Val{Sparsity}=Val(false),saveat=Inf::Float64,abstol=1e-4::Float64,reltol=1e-3::Float64,maxErr=Inf::Float64,maxiters=Int(1e7)::Int) where{PRTYPE,T,D,Z,CS,SolverType,OrderType,Sparsity}     
```

```@docs
solve(prob::NLODEProblem{PRTYPE,T,D,Z,CS},al::QSSAlgorithm{SolverType, OrderType},tspan::Tuple{Float64, Float64};sparsity::Val{Sparsity}=Val(false),saveat=Inf::Float64,abstol=1e-4::Float64,reltol=1e-3::Float64,maxErr=Inf::Float64,maxiters=Int(1e7)::Int) where{PRTYPE,T,D,Z,CS,SolverType,OrderType,Sparsity}     
```

## Algorithms


```@docs
 qss1()  
```
```@docs
 qss2()  
```

```@docs
 liqss1()  
```
```@docs
 liqss2()  
```

```@docs
 nmliqss1()  
```
```@docs
 nmliqss2()  
```

## Solution

### Query a solution


```@docs
QuantizedSystemSolver.evaluateSol(sol::Sol{T,O}, index::Int, t::Float64) where {T,O}
```
```@docs
QuantizedSystemSolver.solInterpolated(sol::Sol{T,O}, t::Float64) where {T,O}
```
```@docs
QuantizedSystemSolver.solInterpolated(sol::Sol{T,O},index::Int,step::Float64) where {T,O}
```

```@docs
QuantizedSystemSolver.show(io::IO, a::Stats)
```

```@docs
QuantizedSystemSolver.Stats
```
### Error with respect to an analytic or reference solution

!!! warning "Finding the error"

    To find the error, an interpolated solution has to be used as shown in [`solInterpolated`](@ref). The reference solution (using saveat for example) and the interpolated solution (passed step to interpolated function) must have the same tspan and the same step size.
    
```@docs
getError(sol::Sol{T,O},index::Int,f::Function) where{T,O}
```
```@docs
getAverageError(sol::Sol{T,O},f::Vector{Function}) where{T,O}
```

```@docs
getErrorByRefs(sol::Sol{T,O},index::Int,solRef::Vector{Any}) where{T,O}
```
```@docs
getAverageErrorByRefs(sol::Sol{T,O},solRef::Vector{Any}) where{T,O}
```


### Plot the solution

```@docs
plot_Sol(sol::Sol{T,O},xvars::Int...;note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool,marker=:circle::Symbol,title="") where{T,O}
```

```@docs
plot(sol::Sol{T,O};idxs=Int[]::Vector{Int},note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool,marker=:circle::Symbol,title="") where{T,O}
``` 


```@docs
save_Sol(sol::Sol{T,O},xvars::Int...;note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool) where{T,O}
```


```@docs
save_SolSum(sol::Sol{T,O},xvars::Int...;interp=0.0001,note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool) where{T,O}
```

```@docs
plot_SolSum(sol::Sol{T,O},xvars::Int...;interp=0.0001,note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool) where{T,O}
```