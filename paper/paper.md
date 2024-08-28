---
title: 'QuantizedSystemSolver: A discontinuous ODE system solver for Julia.'
tags:
  - Quantised State system
  - Discontinuities
  - Stiff Ordinary Differential Equations
authors:
 - name: Elmongi Elbellili
   orcid: 0000-0003-1230-5488
   affiliation: "1, 2"
 - name: Daan Huybrechs
   orcid: 0000-0002-0536-2647
   affiliation: 2
 - name: Ben Lauwens
   orcid: 0000-0003-0761-6265
   affiliation: 1
affiliations:
 - name: Royal Military Academy, Brussels, Belgium
   index: 1
 - name: Kuleuven, Leuven, Belgium
   index: 2
date: August 2024
bibliography: paper.bib
---

# Summary
The growing intricacy of contemporary engineering systems, typically reduced to stiff differential equations with events, poses a difficulty in digitally simulating them using traditional numerical integration techniques. The Quantized State System (QSS) and the Linearly Implicit Quantized State System (LIQSS) are different methods for tackling such problems. It is an approach that builds the solution by updating the system variables independently as opposed to classic integration methods that update all the system variables every step. [@QSS]. 

# Statement of need
Traditional solvers are challenged by frequent discontinuities where the state of the system abruptly changes
at specific points or intervals. They struggle to accurately capture the dynamics around discontinuities. They either undergo expensive iterations to pinpoint exact discontinuity instances or resort to interpolating their locations, resulting in unreliable outcomes. 
Written in the easy-to-learn Julia language [Julia programming language](https://julialang.org) [@julia], 
and taking advantage of its features such as multiple dispatch and metaprogramming, the QuantizedSystemSolver.jl is a solver that aims to efficiently solve a set of Ordinary differential Equations with a set of events via implementing the QSS and LIQSS methods. It is the first such tool to be published in the Julia ecosystem.

# Package description
While the package is optimized to be fast, extensibility is not compromised. It is divided into 3 entities that can be extended separately: Problem, Algorithm, and Solution. The rest of the code is to create these entities and glue them together as shown in the Figure. The API was designed to provide an easier way to handle events than the approach provided by the classic integration solvers. Inside an NLodeProblem function, the user may introduce any parameters, variables, equations, and events:

```
odeprob = NLodeProblem(quote 
          name
          parameters
          discrete and continous variables
          helper expressions
          differential equations
          if-statments for events)
```

The output of this function is an object of type Problem, and it is passed to the solve function along any other configuration arguments such as the algorithm type, the time span and the tolerance. The solve function dispatches on the given algorithm and start the numerical integration. 

```
tspan = (0.0, 0.001)
sol= solve(odeprob,nmliqss2(),tspan,abstol=1e-4,reltol=1e-3)
```

A the end, a solution object is produced that can be queried and plotted. 
```
sol(0.0005,idxs=2) 
sol.stats
plot(sol)
```
![alt text](diagram.png)


In addition, the package contains several shared helper functions used during the integration process by the algorithm such as the scheduler that organizes which variable of the system to update at any specific time of the simulation. 
The solver uses other packages such as  [`MacroTools.jl `]( https://github.com/FluxML/MacroTools.jl)[@MacroTools] for user-code parsing, [`SymEngine.jl `]( https://github.com/symengine/SymEngine.jl)[@SymEngine]  for Jacobian computation and dependencies extraction, and a modified [`TaylorSeries.jl`](https://github.com/JuliaDiff/TaylorSeries.jl/)[@TaylorSeries] that uses caching to obtain free Taylor variable operations as the current version of TaylorSeries creates a heap allocated object for every operation. The approximation through Taylor variables transforms any complicated equations to polynomials, which makes root finding cheaper, which the QSS methods relies heavily on it. 

In conclusion, the package offers extensive functionality to facilitate practical research tasks, all while being lightweight and well documented to be easily used by researchers and students to efficiently model various dynamical systems with discontinuities, as well as further study and improve the newly developed QSS methods. In Fact, Anyone that has a different problem type that can be handled by a specific QSS algorithm that outputs a different solution, they can define their types as subclasses of the abstract types (Problem, Algorithm, and Solution), and use the solver as usual.



# Acknowledgements
This research has received no external funding.

# References
