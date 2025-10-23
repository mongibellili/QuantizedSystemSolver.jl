# Problem Internals


### Problem definition




```@docs
QuantizedSystemSolver.ODEDiscProblem{JACMODE,T,D,Z,CS,F,JAC,CLS}
```



### Problem construction

!!! note "The examples of the continuous problem"

    All the examples that explain the problem construction functions use the following problem. The examples are reproducible as they are shown. To see the problem construction process in one step, add `print()` statements inside the functions of the package while solving this problem as shown in the tutorial section.
      ```julia
    du[1] = u[2]-2.0*u[1]*u[2]
    for k in 2:9  
    du[k]=u[k]*u[k-1];
    end 
    du[10]=u[1]-u[10]
    ```
 


```@docs 
QuantizedSystemSolver.odeProblemFunc(ir::ODEFunctionIR,::Val{T},::Val{D},::Val{0},initCond::Vector{Float64},discrVars::Union{Vector{EM}, Tuple{Vararg{EM}}},preProcessData::PreProcessData,jac_mode ::Symbol) where {T,D,EM} 

```

```@docs
QuantizedSystemSolver.odeProblemFunc(ir::ODEFunctionIR,::Val{T},::Val{D},::Val{Z},initCond::Vector{Float64},discrVars::Union{Vector{EM}, Tuple{Vararg{EM}}},preProcessData::PreProcessData,jac_mode ::Symbol) where {T,D,Z,EM}
```

### Problem construction helpers



```@docs
QuantizedSystemSolver.changeExprToFirstValue(ex::Expr) 
```
```@docs
QuantizedSystemSolver.symbolFromRef(el::Symbol,refEx::Union{Int64,Expr,Symbol}) 
```

```@docs
QuantizedSystemSolver.restoreRef(coefExpr,symDict)
```



```@docs
QuantizedSystemSolver.handleEvents(argI::Expr,eventequs::Vector{Expr},length_zcequs::Int64,evsArr::Vector{EventDependencyStruct})
```

```@docs
QuantizedSystemSolver.transformFSimplecase(ex::Union{Float64,Int64,Expr,Symbol})
```

```@docs
QuantizedSystemSolver.transformF(ex::Expr)
```



```@docs
QuantizedSystemSolver.createJacVect(jac:: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},::Val{T}) where {T}
```

```@docs
QuantizedSystemSolver.createSDVect(jac:: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},::Val{T}) where {T}
```

```@docs
QuantizedSystemSolver.createExactJacFun(otherCode::Expr,Exactjac:: Dict{Expr,Union{Float64,Int,Symbol,Expr}},funName::Symbol,f)
```

```@docs
QuantizedSystemSolver.createContEqFun(otherCode::Expr,equs::Dict{Union{Int,Expr},Union{Int,Symbol,Expr}},fname::Symbol,f)
```

```@docs
QuantizedSystemSolver.extractJacDep(b::Int,niter::Int,rhs::Expr,jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},exactJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},jac_mode ::Symbol,symDict::Dict{Symbol,Expr},dD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 
```



```@docs
QuantizedSystemSolver.extractZCJacDep(counter::Int,zcf::Expr,zcjac :: Vector{Vector{Int}},SZ ::Dict{Int,Set{Int}},dZ :: Dict{Int,Set{Int}})
```

```@docs
QuantizedSystemSolver.EventDependencyStruct
```

```@docs
QuantizedSystemSolver.createSZVect(SZ :: Dict{Int64, Set{Int64}},::Val{T}) where {T} 
```

```@docs
QuantizedSystemSolver.createdDVect(dD::Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}},::Val{D}) where {D}
```

```@docs
QuantizedSystemSolver.createDependencyToEventsDiscr(dD::Vector{Vector{Int}},dZ::Dict{Int64, Set{Int64}},eventDep::Vector{EventDependencyStruct}) 
```

```@docs
QuantizedSystemSolver.createDependencyToEventsCont(SD::Vector{Vector{Int}},sZ::Dict{Int64, Set{Int64}},eventDep::Vector{EventDependencyStruct}) 
```

```@docs
QuantizedSystemSolver.unionDependency(HZD1::Vector{Vector{Int}},HZD2::Vector{Vector{Int}})
```



## Index

```@index
Pages = ["problem.md"]
Order = [:type, :function]
```