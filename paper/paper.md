---
title: 'QuantizedSystemSolver: discontinous ODE system solver for Julia.'
tags:
  - Quantised State system
  - discontinuities
  - Stiff Ordinary differential equations
authors:
 - name: Elmongi Elbellili
   orcid: 0000-0003-1230-5488
   affiliation: "1, 2",
 - name: Daan Huybrechs
   orcid: 0000-0002-0536-2647
   affiliation: 2,
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

The growing intricacy of contemporary engineering systems, typically reduced to differential equations with events, poses a difficulty in digitally simulating them using traditional numerical integration techniques. The Quantized State System (QSS) and the Linearly Implicit Quantized State System (LIQSS) are different methods for tackling such problems [@QSS]. The QuantizedSystemSolver is a package written in the Julia language [Julia programming language](https://julialang.org) [@julia]. that aims to solve a set of Ordinary differential equations with a set of events. It implements the quantized state system methods: An approach that builds the solution by updating the system variables independently as opposed to classic integration methods that update all the system variables every step.

While the package is optimized to be fast, extensibility is not compromised. It is divided into 3 entities that can be extended separately: Problem, Algorithm, and Solution. The API was designed to provide an easier way to handle events than the approach provided by the classic integration solvers. Inside an NLodeProblem function, the user may introduce any parameters, variables, equations, and events. The output of this function (object of type Problem) is sent to the solve function along any other configuration arguments such as the alroithm type, the tolerance, and time span. Finally, a solution is produced that can be queried and plotted. The package uses other packages such as  [`MacroTools.jl `]( https://github.com/FluxML/MacroTools.jl)[@MacroTools] for user-code parsing, [`SymEngine.jl `]( https://github.com/symengine/SymEngine.jl)[@SymEngine]  for Jacobian computation, and a modified [`TaylorSeries.jl`](https://github.com/JuliaDiff/TaylorSeries.jl/)[@TaylorSeries] that uses caching to obtain free Taylor variables. The approximation through Taylor variables transforms any complicated equations to polynomials, which makes root finding cheaper. Finally, It can be used by researchers and students to efficiently model various dynamical systems with discontinuities. 



# References
