# Problem 

The ODEProblem function is the entry point for defining a new problem to be solved by the QSS solver. It takes user-provided code, which includes system parameters, variables, equations, and event logic, and constructs a Problem object that encapsulates all the necessary information for the solver to simulate the system such as problem dimensions, dependencies, and equations. The function works by parsing the user code and extracting relevant data to populate the Problem object.

## Problem extension 
Problem extension can be achieved easily via PRTYPE which is of type Val, or another subtype of this superclass can be created.
```@docs
 QuantizedSystemSolver.ODEProblemData{F,PRTYPE,T,D,Z,CS}
```
### What is needed with a new problem:
The more different the new problem from the `ODEContProblem`, the more functions are needed to be extended. In general the following functions need to be extended.
  - The `odeProblemFunc` method to handle this problem.
  - The `integrate` method for this new type of problem.
  - The `custom_Solve` method if needed.


### Example
```julia
struct SmallODEProblem{CS}<: ODEProblemData{0,1,1,0,0,CS} 
  cacheSize::Val{CS}# CS= cache size 
  initConditions::Float64   
  eq::Function#function that holds the differential equation
end
```
This new problem type takes care of one differential equation. There is no need for the Jacobian nor for the dependencies. This needs an extension of the custom_Solve method that just removes the references to the `jac` and the `SD`. An extension of the integrate method is also needed since the implementation is a lot simpler than what is currently implemented.

## Further reading about the functions creating the problem
*ODEDiscProblem{F,PRTYPE,T,D,Z,CS}:* This is the struct that holds all
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

**odeProblemFunc:** After an initial preparation performed by the The
NLodeProblem function, The function odeProblemFunc takes the resulting
expressions to continue constructing an instance of the ODEDiscProblem
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



*extractZCJacDep:* It Extracts zero-crossing Jacobian dependencies
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
