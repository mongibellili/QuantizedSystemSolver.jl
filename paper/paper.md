---
title: 'QuantizedSystemSolver: A discontinuous ODE system solver in Julia.'
tags:
  - Quantised State System
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
 - name: KU Leuven, Leuven, Belgium
   index: 2
date: August 2024
bibliography: paper.bib
---

# Summary

Contemporary engineering systems, such as electrical circuits, mechanical systems with shocks, and chemical reactions with rapid kinetics, are often characterized by dynamics that can be modeled using stiff differential equations with events. Stiffness typically arises in these systems due to the presence of both rapidly changing and slowly changing components. This stiffness requires extremely small time steps to maintain stability when using traditional numerical integration techniques. Recently, quantization-based techniques have emerged as an effective alternative for handling such complex models. Methods like the Quantized State System (QSS) and the Linearly Implicit Quantized State System (LIQSS) offer promising results, particularly for large sparse stiff models. Unlike classic numerical integration methods, which update all system variables at each time step, the quantized approach updates individual system variables independently. Specifically, in quantized methods, each variable is updated only when its value changes by a predefined quantization level. Moreover, these methods are advantageous when dealing with discontinuous events. An event is a discontinuity where the state of the system abruptly changes at a specific point. Classic methods may struggle with events: They either undergo expensive iterations to pinpoint the exact discontinuity instance or resort to interpolating its location, resulting in unreliable outcomes. Therefore, this QSS strategy can significantly reduce computational effort and improve efficiency in large sparse stiff models with frequent discontinuities [@improveliqss].

# Statement of need

Traditional solvers are challenged by large sparse stiff models and systems with frequent discontinuities. The buck converter is a stiff system with frequent discontinuities that classic solvers from the DifferentialEquations.jl [@Rackauckas2017] are currently unable to handle properly.  
Written in the easy-to-learn Julia programming language [@julia] and inspired by the ``qss-solver`` written in C [@qssC], the QuantizedSystemSolver.jl package takes advantage of Julia features such as multiple dispatch and metaprogramming. The package shares the same interface as the DifferentialEquations.jl package and aims to efficiently solve a large set of stiff Ordinary Differential Equations (ODEs) with events by implementing the QSS and LIQSS methods. It is the first such tool to be published in the Julia ecosystem. 

# Quantization-based techniques

The general form of a problem composed of a set of ODEs and a set of events that QSS is able to solve is described in the following equations: 

$\dot X = f(X,P,t)$; if $zc(X,P,t)$ ; Set $x_i=H(X,P,t)$ and $p_j=L(X,P,t)$

where $X = [x_1,x_2...,x_n]^T$ is the state vector, $f:\mathbb{R}^n.\; {R}^m. \;{R}^+ \rightarrow \mathbb{R}^n$ is the derivative function, and $t$ is the independent variable. $P = [p_1,p_2...,p_m]^T$ is the vector of the system discrete variables. $n$ and $m$ are the number of state variables and discrete variables of the system respectively. $zc$ is an event condition, $H$ and $L$ are functions used in the effects of the event $zc$.

In QSS, besides the step size, the difference between $x_i(t_k)$ (the current value) and $x_i(t_{k+1})$ (the next value) is called the quantum $\Delta_i$. Depending on the type of the QSS method (explicit or implicit), a new variable $q_i$ is set to equal $x_i(t_k)$  or $x_i(t_{k+1})$ respectively. $q_i$ is called the quantized state of $x_i$, and it is used in updating the derivative function [@elbellili].  A general description of a QSS algorithm is given as follows:
![](alg.png)

# Package description

While the package is optimized to be fast, extensibility is not compromised. It is divided into three entities that can be extended separately: The ``problem``, the ``algorithm``, and the ``solution``. The rest of the code creates these entities and glue them together. The API was designed to match the differentialEquations.jl interface while providing an easier way to handle events. The problem is defined inside a function, in which the user may introduce any parameters, variables, equations, and events:
```julia
function func(du,u,p,t) 
  # parameters, helpers, differential eqs., if-statements for events; e.g.:
  du[1] = p[1] * u[1]
  if (t > 1) p[1] = 42.0 end
end
```
Then, this function is passed to an `ODEProblem` function along with the initial conditions, the time span, and any parameters or discrete variables. 
```julia
tspan = (initial_time,final_time)
u = [10.0] #initial conditions
p = [-1.0] # parameters and discrete variables
odeprob=ODEProblem(func,u,tspan,p)
```

The output of the previous `ODEProblem` function, which is a QSS problem, is passed to a ``solve`` function with other configuration arguments such as the algorithm type and the tolerance. The ``solve`` function dispatches on the given algorithm and starts the numerical integration.
```julia
sol = solve(odeprob,algorithm,abstol = ...,reltol = ...)    
```
At the end, a solution object is produced that can be queried, plotted, and analyzed for error.

```julia
sol(0.05,idxs = 1) # get the value of variable 1 at time 0.05
sol.stats          # get statistics about the simulation
plot(sol)          # plot the solution
```

The solver uses other packages such as  [`MacroTools.jl`]( https://github.com/FluxML/MacroTools.jl) [@MacroTools] for user-code parsing, [`SymEngine.jl`]( https://github.com/symengine/SymEngine.jl) [@SymEngine] for Jacobian computation and dependency extraction. It also uses and a modified [`TaylorSeries.jl`](https://github.com/JuliaDiff/TaylorSeries.jl/) [@TaylorSeries] that implements caching to obtain free Taylor variable operations, since the current version of TaylorSeries creates a heap allocated object for every operation. The approximation through Taylor variables transforms any complicated equations to polynomials, making root finding cheaper--a process that QSS methods rely on heavily. 

# The buck converter example

The Buck is a converter that decreases voltage and increases current with a greater power efficiency than linear regulators [@buck]. Its circuit is shown in Fig.1(a).

![The buck converter](buck.png)

The diode $D$ and the switch $S$ can be modeled as two variables resistors $RD$ and $RS$. A mesh and a nodal analysis give the relationship between the different components in the circuit as follows:

$i_d = \frac{RS.i_l-V1}{RS+RD}$; 
$\frac{du_c}{dt} = \frac{i_l-\frac{u_c}{R}}{C}$; 
$\frac{di_l}{dt} = \frac{-uc-i_d.RD}{L}$

The buck problem contains frequent discontinuities and can be solved by the QuantizedSystemSolver.jl package using the following code:

```julia
using QuantizedSystemSolver
function buck(du,u,p,t)
  #Constant parameters
  C = 1e-4; L = 1e-4; R = 10.0; V1 = 24.0; T = 1e-4; ROn = 1e-5; ROff = 1e5
  #Optional rename of the continuous and discrete variables
  RD = p[1]; RS = p[2]; nextT = p[3]; lastT = p[4]; il = u[1]; uc = u[2]
  #Equations
  id = (il*RS-V1)/(RD+RS) # diode's current
  du[1] = (-id*RD-uc)/L; du[2] = (il-uc/R)/C
  #Events
  if t-nextT > 0.0 # model when the switch is ON
    lastT = nextT; nextT = nextT+T; RS = ROn
  end
  if t-lastT-0.5*T > 0.0 # model when the switch is OFF
    RS = ROff
  end
  if id > 0 # model when the Diode is ON
    RD = ROn;
  else
    RD = ROff;
  end
end
# Initial conditions and time settings
p = [1e5,1e-5,1e-4,0.0]; u0 = [0.0,0.0]; tspan = (0.0,0.001)
# Define the problem
QSSproblem = ODEProblem(buck,u0,tspan,p)
# solve the problem
sol = solve(QSSproblem,nmliqss2(),abstol = 1e-3,reltol = 1e-2)
# Get the value of variable 2 at time 0.0005 
sol(0.0005,idxs = 2)
# plot the solution
plot(sol)
```
The solution plot is presented in Fig.1(b).


# Conclusion
The package provides robust functionality to efficiently solve stiff ODEs with events using the quantized state method. It is well-documented, making it accessible for researchers across various domains. Additionally, users can extend its capabilities to handle a variety of problems.

# Acknowledgements
This research has received no external funding.

# References
