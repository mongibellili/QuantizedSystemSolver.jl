# User API


## Problem definition

```@docs
ODEProblem(f::Function, u::Vector{Float64}, tspan::Tuple{A,B}, p::Union{Vector{EM}, Tuple{Vararg{EM}}};jac_mode ::Symbol= :symbolic) where{EM,A<:Union{Float64, Int64},B<:Union{Float64, Int64}}
```
```@docs
ODEProblem(f::Function, u::Vector{Float64}, tspan::Tuple{A,B};jac_mode ::Symbol= :symbolic) where{A<:Union{Float64, Int64},B<:Union{Float64, Int64}}
```

 


## The solve function:

```@docs
solve(prob::ODEProblemData{F,PRTYPE,T,D,Z,CS},al::QSSAlgorithm{SolverType, OrderType};detection::Val{M}=Val(2),saveat::Float64=Inf,abstol::Float64=1e-3,reltol::Float64=1e-3,maxErr::Float64=1.0,maxiters::Int=Int(1e7),verbose=false::Bool) where{F,PRTYPE,T,D,Z,CS,SolverType,OrderType,M}        
```


```@docs
solve(prob::ODEProblemData{F,PRTYPE,T,D,Z,CS},al::QSSAlgorithm{SolverType, OrderType},tspan::Tuple{Float64, Float64};detection::Val{M}=Val(2),saveat::Float64=Inf,abstol::Float64=1e-3,reltol::Float64=1e-3,maxErr::Float64=1.0,maxiters::Int=Int(1e7),verbose=false::Bool) where{F,PRTYPE,T,D,Z,CS,SolverType,OrderType,M} 
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
QuantizedSystemSolver.evaluateSol(sol::Sol{T,O},index::Int,t::Float64) where {T,O}
```
```@docs
QuantizedSystemSolver.solInterpolated(sol::Sol{T,O},step::Float64) where {T,O}
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
plot_SolSum(sol::Sol{T,O},xvars::Int...;interp=0.0001,note=" "::String,xlims::Tuple{Float64, Float64}=(0.0,0.0),ylims::Tuple{Float64, Float64}=(0.0,0.0),legend::Bool=true) where{T,O}
```

