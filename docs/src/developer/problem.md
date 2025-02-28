# Problem 

The NLodeProblem function is the entry point for defining a new problem to be solved by the QSS solver. It takes user-provided code, which includes system parameters, variables, equations, and event logic, and constructs a Problem object that encapsulates all the necessary information for the solver to simulate the system such as problem dimensions, dependencies, and equations. The function works by parsing the user code and extracting relevant data to populate the Problem object.

## Problem extension
Problem extension can be achieved easily via PRTYPE which is of type Val, or another subtype of this superclass can be created.
```@docs
 QuantizedSystemSolver.NLODEProblem{F,PRTYPE,T,D,Z,CS}
```
### What is needed with a new problem:
The more different the new problem from the `NLODEContProblem`, the more functions are needed to be extended. In general the following functions need to be extended.
  - The `NLodeProblemFunc` method to handle this problem.
  - The `integrate` method for this new type of problem.
  - The `custom_Solve` method if needed.


### Example
```julia
struct SmallODEProblem{CS}<: NLODEProblem{0,1,1,0,0,CS} 
  cacheSize::Val{CS}# CS= cache size 
  initConditions::Float64   
  eq::Function#function that holds the differential equation
end
```
This new problem type takes care of one differential equation. There is no need for the Jacobian nor for the dependencies. This needs an extension of the custom_Solve method that just removes the references to the `jac` and the `SD`. An extension of the integrate method is also needed since the implementation is a lot simpler than what is currently implemented.

## Further reading about the functions creating the problem
*NLODEDiscProblem{F,PRTYPE,T,D,Z,CS}:* This is the struct that holds all
the necessary data for a nonlinear ordinary differential equation (ODE) 
problem with discrete events. The structure includes various fields such
as initial conditions, discrete variables, Jacobians, event
dependencies, and other data related to how the problem is formulated.
This structure serves as the core data holder for the problem and will
be used in the solver. It is a parametric abstract type that has the
following parameters:

PRTYPE: The type of the problem (to distinguish between various types,
and allow future extension of the solver to handle new types).

T: The number of continuous variables (state variables).

Z: The number of zero-crossing functions, which are used to detect
events.

Y: The actual number of events.

CS: Cache size, which is used to store intermediate operations.

The use of abstract types in this context allows for flexibility and
extensibility in the solver. By defining these abstract types, the code
can be easily adapted to handle different types of problems, algorithms,
and solutions without needing to modify the core solver logic. This
design choice enhances the maintainability and scalability of the
solver, making it easier to add new features or support additional
problem types in the future.

**NLodeProblemFunc:** After an initial preparation performed by the The
NLodeProblem function, The function NLodeProblemFunc takes the resulting
expressions to continue constructing an instance of the NLODEDiscProblem
structure. It works in several key stages:

*Initialization:* The function begins by initializing vectors and
dictionaries that will hold equations (equs), Jacobian dependencies
(jac), zero-crossing functions (ZCjac), and event dependencies. These
serve to store the different types of equations and their relationships.

*Processing ODEs:* It loops through each of the ODE expressions provided
by the user. Depending on the type of expression (discrete variables,
differential equations, or loop constructs), it processes the right-hand
side (RHS) of the equation. For differential equations, it extracts
dependencies to build the Jacobian and transform the equations into a
more appropriate form for further use. Special cases are handled, such
as if the RHS is a number or a symbol.

*Handling Events:* The function also processes event-related constructs
(if conditions) that correspond to different points where the system
might undergo discrete changes. It process the RHS of the event
equations, transforms them into a suitable form, and builds the
necessary dependency structures. Specifically, it constructs how
discrete and continuous variables influence one another through the
events.

*Constructing the Function Code:* After processing all ODEs and events,
the function dynamically generates a Julia function code needed to store
the system of ODEs and events. This code is built into a function that
handles different cases (i.e., which equation to evaluate based on an
index of a state change or an event).

**Building Dependencies:** Several helper functions that build the
dependencies between variables, events. They build dependency vectors 
that track how discrete and continous variables influence the system.
This is used to know what variables to update and determine when
specific events should be checked. By tracking the relationships between
variables and events, the solver can determine the appropriate actions
to take at each time step. The dependencies are stored in the following
vectors:

\-$jac$: It determines which variables affect a derivative.

\-$ZCjac$: It determines which variables affect a zero-crossing
function.

\-$SD$: It determines which derivatives that are affected by a given
variable.

\-$SZ$: It determines which zero-crossing functions that are affected by
a given variable.

\-$HZ$: It tells which Zero-crossing functions influenced by a given
event.

\-$HD$: It tells which derivatives influenced by a given event.

Here's a quick summary and what each helper function is doing:

*extractJacDepNormal:* It Extracts the dependencies for normal
(non-loop) expressions. It updates the Jacobian matrix $jac$ and a
dictionary $dD$ for tracking dependencies of derivatives to discrete
variables.

*extractJacDepLoop:* Similar to extractJacDepNormal, but specifically
for loop expressions. It tracks dependencies across loop iterations.

*extractZCJacDepNormal:* It Extracts zero-crossing Jacobian dependencies
for discrete variables ($dZ$), and it updates $zcjac$, $SZ$.

*createDependencyToEventsDiscr:* It maps discrete dependencies (dD, dZ)
to specific events, it and constructs dependency matrices HZ and HD from
the discrete variables only.

*createDependencyToEventsCont:* Similar to
createDependencyToEventsDiscr, but for continuous dependencies (SD, sZ),
and it updates the matrices HZ and HD from the continuous variables
only.

*unionDependency:* Merges the two previous sets of dependencies
(continuous and discrete) into the final matrices HZ and HD.

## Helper packages
The 𝑝𝑜𝑠𝑡𝑤𝑎𝑙𝑘 function from the MacroTools.jl (Copyright (c) 2015: Mike Innes) package plugs parameters and helper functions directly into the equations, and traverses the right-hand side of differential equations and zero-crossing functions, facilitating the construction of the Jacobian matrix and identifying variable dependencies. It also transforms specific expressions like 𝑞[1] into 𝑞[1] [0] within events, and converts 𝑞[𝑖] to 𝑞𝑖, making the equations more tractable for differentiation and Jacobian construction. Additionally, the @𝑐𝑎𝑝𝑡𝑢𝑟𝑒 macro efficiently handles cases where differential equations are defined within a for loop. 

The diff(basi, symarg) function from the SymEngine.jl (Copyright (c) 2015-2017 Isuru Fernando) package is applied to perform symbolic differentiation, where basi is an expression and symarg is the symbol with respect to which the derivative is taken. This returns the partial derivative of the expression, making it particularly useful for deriving system Jacobians. 

@code_string macro from the CodeTracking.jl (Copyright (c) 2019 Tim Holy) is used to get the body expression of the function that holds the problem given by the user. 

@RuntimeGeneratedFunction from the RuntimeGeneratedFunctions.jl package (Copyright (c) 2020 Chris Rackauckas) is used to avoid world-age issues with the generated functions.

## Internals

### Problem definition

```@docs
QuantizedSystemSolver.NLODEContProblem{F,PRTYPE,T,D,Z,CS}
```
```@docs
QuantizedSystemSolver.NLODEContProblemSpan{PRTYPE,T,Z,Y,CS}
```

```@docs
QuantizedSystemSolver.NLODEDiscProblem{F,PRTYPE,T,D,Z,CS}
```


```@docs
QuantizedSystemSolver.NLODEDiscProblemSpan{PRTYPE,T,Z,D,CS}
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
QuantizedSystemSolver. NLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{D},::Val{0},initConditions::Vector{Float64},du::Symbol,tspan::Tuple{Float64, Float64},discrVars::Union{Vector{EM}, Tuple{Vararg{EM}}},prbName::Symbol) where {T,D,EM} 
```

```@docs
QuantizedSystemSolver.NLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{D},::Val{Z},initCond::Vector{Float64},du::Symbol,tspan::Tuple{Float64, Float64},discrVars::Union{Vector{EM}, Tuple{Vararg{EM}}},prbName::Symbol) where {T,D,Z,EM} 
```

### Problem construction helpers
```@docs
QuantizedSystemSolver.prepareInfo(odeExprs::Expr,stateVarName::Symbol,discrParamName::Symbol) 
```
```@docs
QuantizedSystemSolver.probHelper
```
```@docs
QuantizedSystemSolver.arrangeProb(x::Expr)
```



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
QuantizedSystemSolver.changeVarNames_params(ex::Expr,stateVarName::Symbol,discrParamName::Symbol,muteVar::Symbol,param::Dict{Symbol,Union{Float64,Int64,Expr,Symbol}})
```


```@docs
QuantizedSystemSolver.changeVarNames_params(element::Symbol,stateVarName::Symbol,discrParamName::Symbol,muteVar::Symbol,param::Dict{Symbol,Union{Float64,Int64,Expr,Symbol}})
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
QuantizedSystemSolver.extractJacDepNormal(varNum::Int,rhs::Union{Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}, exactJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},symDict::Dict{Symbol,Expr})
```
```@docs
QuantizedSystemSolver.extractJacDepLoop(b::Int,niter::Int,rhs::Union{Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}} ,exactJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},symDict::Dict{Symbol,Expr})
```
```@docs
QuantizedSystemSolver.createJacVect(jac:: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},::Val{T}) where {T}
```
```@docs
QuantizedSystemSolver.createSDVect(jac:: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},::Val{T}) where {T}
```
```@docs
QuantizedSystemSolver.createExactJacFun(otherCode::Expr,Exactjac:: Dict{Expr,Union{Float64,Int,Symbol,Expr}},funName::Symbol,f::F) where{F} 
```
```@docs
QuantizedSystemSolver.createContEqFun(otherCode::Expr,equs::Dict{Union{Int,Expr},Union{Int,Symbol,Expr}},fname::Symbol,f::F) where {F}
```
```@docs
QuantizedSystemSolver.extractJacDepNormalDiscrete(varNum::Int,rhs::Union{Symbol,Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},exactJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},symDict::Dict{Symbol,Expr},dD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 
```
```@docs
QuantizedSystemSolver.extractJacDepLoopDiscrete(b::Int,niter::Int,rhs::Union{Symbol,Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},exactJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},symDict::Dict{Symbol,Expr},dD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 
```

```@docs
QuantizedSystemSolver.extractZCJacDepNormal(counter::Int,zcf::Expr,zcjac :: Vector{Vector{Int}},SZ ::Dict{Int,Set{Int}},dZ :: Dict{Int,Set{Int}})
```

```@docs
QuantizedSystemSolver.EventDependencyStruct
```

```@docs
QuantizedSystemSolver.createSZVect(SZ :: Dict{Int64, Set{Int64}},::Val{T}) where{T} 
```

```@docs
QuantizedSystemSolver.createdDVect(dD::Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}},::Val{D}) where{D}
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