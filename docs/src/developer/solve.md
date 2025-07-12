# Solve

The solve function is the primary interface for solving ODE problems using various QSS (Quantized State Systems) algorithms. It does
this by dispatching on the problem and algorithm types to select the right solver.

Based on the QSS Algorithm provided, the solve function either selects a basic QSS integration method or, in the case of LiQSS, constructs additional data
structures needed for implicit integration. The method defaults to cycle detection (Val(2)), tolerances (abstol=1e-3 and reltol=1e-3), and a maximum number of iterations (maxiters=1e7). These parameters can be adjusted based on the problem's complexity and desired accuracy.
When a user doesn't provide a solver explicitly, the solve function defaults to using the modified second-order implicit algorith (mLiQSS2). 

## Helper Functions:

### createCommonData:
Sets up the common data required by the QSS solver. It initializes all necessary vectors (x, q, tx, tq,
nextStateTime, etc.) used in the integration process, and it pre-allocates a cache of Taylor series to avoid repeated memory allocation during integration (taylorOpsCache).

### getClosure:
Creates a closure for the Jacobian and the dependency matrices, allowing easy
extension. The closure is used to access the Jacobian matrix during the integration process.

### createLiqssData:
Initializes internal data structures required by LiQSS solvers for managing linear approximations.
This function has two methods, depending on how the Jacobian is specified during problem definition (via the `ODEProblem` function):

**Symbolic Jacobian (`:symbolic`):**

```julia
  createLiqssData(::Val{:symbolic}, ::Val{M}, ::Val{T}, ::Val{Order}) where {T, Order, M}
```
Constructs and returns an AexprLiQSS_data object, which stores symbolic expressions for the Jacobian entries.

**Approximated Jacobian (`:approximate`):**
```julia
createLiqssData(::Val{:approximate}, ::Val{M}, ::Val{T}, ::Val{Order}) where {T, Order, M}
```
Constructs and returns an AmanualLiQSS_data object, which stores coefficient data to be used for the approximation of the jacobian entries.

## Internals

```@docs
QuantizedSystemSolver.custom_Solve(prob::ODEProblemData{F,PRTYPE,T,D,Z,CS},al::QSSAlgorithm{Solver, Order},::Val{M},finalTime::Float64,saveat::Float64,initialTime::Float64,abstol::Float64,reltol::Float64,maxErr::Float64,maxiters::Int,verbose::Bool) where{F,PRTYPE,T,D,Z,CS,Solver,Order,M}    
```


```@docs 
QuantizedSystemSolver.createCommonData(prob::ODEProblemData{F,PRTYPE,T,D,Z,CS},::Val{Order},finalTime::Float64,saveat::Float64,initialTime::Float64,abstol::Float64,reltol::Float64,maxErr::Float64,maxiters::Int,verbose::Bool) where{F,PRTYPE,T,D,Z,CS,Order}
```

```@docs
QuantizedSystemSolver.createLiqssData(::Val{:symbolic},::Val{M},::Val{T},::Val{Order}) where{T,Order,M} 
```



```@docs
QuantizedSystemSolver.CommonQSS_Data{Z}
```

  

```@docs 
 QuantizedSystemSolver.LiQSS_Data{O,M}
```

```@docs 
 QuantizedSystemSolver.AexprLiQSS_data{O,M}
```

```@docs 
 QuantizedSystemSolver.AmanualLiQSS_data{O,M}
```