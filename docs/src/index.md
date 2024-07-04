# Quantized System Solver

```@contents
Pages = ["Taylor0.md", "Problem.md","Algorithm.md","examples.md"]
Depth = 1
```
## Algorithms:
-  [`qss1()`](@ref),[`qss2()`](@ref)



```@docs
 solve(prob::NLODEProblem{PRTYPE,T,Z,D,CS},tspan::Tuple{Float64, Float64};sparsity::Val{Sparsity}=Val(false)::Float64,saveat=1e-9::Float64::Float64,abstol=1e-4::Float64,reltol=1e-3::Float64,maxErr=Inf::Float64) where {PRTYPE,T,Z,D,CS,Sparsity}    
```

- back to top [Quantized System Solver](@ref)
