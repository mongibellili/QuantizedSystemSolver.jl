# QuantizedSystemSolver
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mongibellili.github.io/QuantizedSystemSolver/dev/)
[![CI](https://github.com/mongibellili/QuantizedSystemSolver/actions/workflows/CI.yml/badge.svg)](https://github.com/mongibellili/QuantizedSystemSolver/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/mongibellili/QuantizedSystemSolver/branch/main/graph/badge.svg)](https://codecov.io/gh/mongibellili/QuantizedSystemSolver)


The QuantizedSystemSolver implements the quantized state system methods: An approach that builds the solution by updating the system variables independently as opposed to classic integration methods that update all the system variables every step.
## Problem
the NLodeProblem function takes the user code as shown in the example below:
```
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
The output is an object that subtypes 
``` abstract type NLODEProblem{PRTYPE,T,Z,Y,CS} end ```

## solve
Input: NLODEProblem{PRTYPE,T,Z,Y,CS}
       ALGORITHM{N,O}
 ```
tspan = (0.0, 0.001)
 sol= solve(odeprob,nmliqss2(),tspan,abstol=1e-4,reltol=1e-3)
```
Output:Sol{T,O}

## Query the solution

```xp=sol(2,0.0005)```

## Plot the solution

```getPlot(sol)```
```save_Sol(sol)```

## Compute the error against a reference solution
``` error=getAverageErrorByRefs(solRef::Vector{Any},sol)```

