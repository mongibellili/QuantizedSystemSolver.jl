# Developper Guide

While the package is optimized to be fast, extensibility is not compromised. It is divided into 3 entities that can be extended separately: Problem, Algorithm, and Solution. The package uses other packages such as MacroTools.jl for user-code parsing, SymEngine.jl for Jacobian computation, and a modified TaylorSeries.jl that uses caching to obtain free Taylor variables. The approximation through Taylor variables transforms any complicated equations to polynomials, which makes root finding cheaper.




### [Problem extension](./Problem.md#problem)

### [Algorithm extension](./algorithmDev.md#algorithm-extension)

### [Solution extension](./solutionDev.md#solution-extension)

### [Taylor0 extension](./Taylor0.md#taylor0)



