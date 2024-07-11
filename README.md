# QuantizedSystemSolver
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mongibellili.github.io/QuantizedSystemSolver.jl/dev/)
[![CI](https://github.com/mongibellili/QuantizedSystemSolver/actions/workflows/CI.yml/badge.svg)](https://github.com/mongibellili/QuantizedSystemSolver/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/mongibellili/QuantizedSystemSolver/branch/main/graph/badge.svg)](https://codecov.io/gh/mongibellili/QuantizedSystemSolver)

The growing intricacy of contemporary engineering systems, typically reduced to differential equations with events, poses a difficulty in digitally simulating them using traditional numerical integration techniques. The Quantized State System (QSS) and the Linearly Implicit Quantized State System (LIQSS) are different methods for tackling such problems.
The QuantizedSystemSolver aims to solve a set of Ordinary differential equations with a set of events. It implements the quantized state system methods: An approach that builds the solution by updating the system variables independently as opposed to classic integration methods that update all the system variables every step.
# Example: Buck circuit
[The Buck](https://en.wikipedia.org/wiki/Buck_converter) is a converter that decreases voltage and increases current with a greater power efficiency than linear regulators. After a mesh analysis we get the problem discribed below.
<img width="220" alt="buckcircuit" src="https://github.com/mongibellili/QuantizedSystemSolver/assets/59377156/c0bcfdbe-ed12-4bb0-8ad1-649ae72dfdd2">

## Problem
The NLodeProblem function takes the following user code:
```julia
    odeprob = NLodeProblem(quote
          name=(buck,)
          #parameters
          C = 1e-4; L = 1e-4; R = 10.0;U = 24.0; T = 1e-4; DC = 0.5; ROn = 1e-5;ROff = 1e5;
          #discrete and continous variables
          discrete = [1e5,1e-5,1e-4,0.0,0.0];u = [0.0,0.0]
          #rename for convenience
          rd=discrete[1];rs=discrete[2];nextT=discrete[3];lastT=discrete[4];diodeon=discrete[5]
          il=u[1] ;uc=u[2]
          #helper equations
          id=(il*rs-U)/(rd+rs) # diode's current
          #differential equations
          du[1] =(-id*rd-uc)/L
          du[2]=(il-uc/R)/C
          #events 
          if t-nextT>0.0 
            lastT=nextT;nextT=nextT+T;rs=ROn
          end
          if t-lastT-DC*T>0.0 
            rs=ROff
          end                          
          if diodeon*(id)+(1.0-diodeon)*(id*rd-0.6)>0
            rd=ROn;diodeon=1.0
          else
            rd=ROff;diodeon=0.0
          end     
    end)
```
The output is an object that subtypes the

```julia 
abstract type NLODEProblem{PRTYPE,T,Z,Y,CS} end
```

## Solve
The solve function takes the previous problem (NLODEProblem{PRTYPE,T,Z,Y,CS}) with a chosen algorithm (ALGORITHM{N,O}) and some simulation settings:

 ```julia
tspan = (0.0, 0.001)
 sol= solve(odeprob,nmliqss2(),tspan,abstol=1e-4,reltol=1e-3)
```
It outputs a solution of type Sol{T,O}.

## Query the solution

```julia
# The value of variable 2  at time 0.0005
sol(2,0.0005)
19.209921627620943
# The total number of steps to end the simulation
sol.totalSteps
498
# The number of simultaneous steps during the simulation
sol.simulStepCount
132
# The total number of events during the simulation
sol.evCount
53
# The actual data is stored in two vectors:
sol.savedTimes
sol.savedVars
```

## Plot the solution

```julia
# If the user wants to perform other tasks with the plot:
plot_Sol(sol)
```
```julia
# If the user wants to save the plot to a file:
save_Sol(sol)
```
![plot_buck_nmLiqss2_()_0 0001_ _ft_0 001_2024_7_2_16_43_54](https://github.com/mongibellili/QuantizedSystemSolver/assets/59377156/00bee649-d337-445a-9ddb-9669076d8ffa)

## Error Analysis

```julia
#Compute the error against an analytic solution
error=getError(sol::Sol{T,O},index::Int,f::Function)
#Compute the error against a reference solution
error=getAverageErrorByRefs(sol,solRef::Vector{Any})
```

