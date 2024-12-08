# Solve

The solve function is the primary interface for solving non-linear ODE problems using various QSS (Quantized State Systems) algorithms. It does
this by dispatching on the problem and algorithm types to select the right solver.

Based on the QSS Algorithm provided, the solve function either selects a basic QSS integration method or, in the case of LiQSS, constructs additional data
structures needed for implicit integration. The method defaults to sparse handling as false (Val(false)), tolerances (abstol=1e-4 and reltol=1e-3), and a maximum number of iterations (maxiters=1e7). These parameters can be adjusted based on the problem's complexity and desired accuracy.
When a user doesn't provide a solver explicitly, the solve function defaults to using the modified second-order implicit algorith (mLiQSS2). 

## Helper Functions:

*createCommonData:* Sets up the common data required by the QSS solver. It initializes all necessary vectors (x, q, tx, tq,
nextStateTime, etc.) used in the integration process, and it pre-allocates a cache of Taylor series to avoid repeated memory allocation during integration (taylorOpsCache).

*getClosure:* Creates a closure for the Jacobian and the dependency matrices, allowing in a flexible in a way that enables easy
extension. The closure is used to access the Jacobian matrix during the integration process.

*createLiqssData:* For LiQSS solvers, this function sets up additional data structures and auxiliary variables used for storing linear
approximation coefficients.


## Internals

```@docs
QuantizedSystemSolver.custom_Solve(prob::NLODEProblem{PRTYPE,T,Z,D,CS}, al::QSSAlgorithm{Solver, Order}, ::Val{Sparsity}, finalTime::Float64, saveat::Float64, initialTime::Float64, abstol::Float64, reltol::Float64, maxErr::Float64, maxiters::Int) where {PRTYPE,T,Z,D,CS,Solver,Order,Sparsity}     
```


```@docs
QuantizedSystemSolver.createCommonData(prob::NLODEProblem{PRTYPE,T,Z,D,CS}, ::Val{Order}, finalTime::Float64, saveat::Float64, initialTime::Float64, abstol::Float64, reltol::Float64, maxErr::Float64, maxiters::Int) where {PRTYPE,T,Z,D,CS,Order}
```

```@docs
QuantizedSystemSolver.createLiqssData(prob::NLODEProblem{PRTYPE,T,Z,D,CS},::Val{false},::Val{T},::Val{Order}) where{PRTYPE,T,Z,D,CS,Order}  
```



```@docs
QuantizedSystemSolver.CommonQSS_Data{Z}
```



```@docs 
 QuantizedSystemSolver.LiQSS_Data{O,Sparsity}
```
