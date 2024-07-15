
# Application Programming Interface

```@docs
NLodeProblem(odeExprs) 
```
```@docs
solve(prob::NLODEProblem{PRTYPE,T,Z,D,CS},al::QSSAlgorithm{SolverType, OrderType},tspan::Tuple{Float64, Float64};sparsity::Val{Sparsity}=Val(false),saveat=1e-9::Float64,abstol=1e-4::Float64,reltol=1e-3::Float64,maxErr=Inf::Float64,maxiters=10000000) where{PRTYPE,T,Z,D,CS,SolverType,OrderType,Sparsity}     
```
```@docs
plot_Sol(sol::Sol{T,O},xvars::Int...;note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool,marker=:circle::Symbol,title="") where{T,O}
```

```@docs
save_Sol(sol::Sol{T,O},xvars::Int...;note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool) where{T,O}
```
```@docs
getError(sol::Sol{T,O},index::Int,f::Function) where{T,O}
```
```@docs
getAverageErrorByRefs(sol::Sol{T,O},solRef::Vector{Any}) where{T,O}
```
