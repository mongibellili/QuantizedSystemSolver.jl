# <img width="50" height="50" style="position:relative; top:5px" src="docs\src\logo.png"> QuantizedSystemSolver



| **Documentation** |**Build Status** | **Code coverage** |**Citation** |
|:------------ |------------|------------|------------|
| [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mongibellili.github.io/QuantizedSystemSolver.jl/dev/) | [![CI](https://github.com/mongibellili/QuantizedSystemSolver/actions/workflows/CI.yml/badge.svg)](https://github.com/mongibellili/QuantizedSystemSolver/actions/workflows/CI.yml)|[![Coverage](https://codecov.io/gh/mongibellili/QuantizedSystemSolver/branch/main/graph/badge.svg)](https://codecov.io/gh/mongibellili/QuantizedSystemSolver)|[![DOI](https://zenodo.org/badge/729081556.svg)](https://doi.org/10.5281/zenodo.14361142)|


Contemporary engineering systems, such as electrical circuits, mechanical systems with shocks, and chemical reactions with rapid kinetics, are often characterized by dynamics that can be modeled using stiff differential equations with events. Recently, quantization-based techniques have emerged as an effective alternative for handling such complex models. Methods like the Quantized State System (QSS) and the Linearly Implicit Quantized State System (LIQSS) offer promising results, particularly for large sparse stiff models. Unlike classic numerical integration methods, which update all system variables at each time step, the quantized approach updates individual system variables independently. Moreover, these methods are advantageous when dealing with discontinuous events, where traditional integrators may struggle with accuracy.  
The QuantizedSystemSolver aims to solve a set of Ordinary differential equations with a set of events. It implements the quantized state system methods.
#   <span style="color:red">Installation</span>
QuantizedSystemSolver.jl is a [registered package](http://pkg.julialang.org), and is
simply installed by the foloowing:

Run julia in the terminal, then enter ] to bring up Julia's package manager, and add the QuantizedSystemSolver.jl package:
```console
julia

julia> ]
 
(@v1.x) pkg> add QuantizedSystemSolver
```
The general form of a problem composed of a set of ODEs and a set of events that QSS is able to solve is described in the following: 

**System of $n$ ODEs:**

```math
\begin{align*}
  & \dot X = f(X,D,t) , 
\end{align*}
```

**System of $v$ events:**
```math
\begin{align*}
& if \; zc_v(x_i...,d_p...,t) \; i \in [1,n]  \;  ; \; p  \in [1,m] \\
& \qquad x_i = H(x_i...,d_p...,t) \\
& \qquad \qquad...\\
& \qquad d_p = L(x_i...,d_p...,t)  \\
& \qquad \qquad...\\
\end{align*}
```
where $X = [x_1,x_2...,x_n]^T$ is the state vector, $f:\mathbb{R}^n \rightarrow \mathbb{R}^n$ is the derivative function, and $t$ is the independent variable. $D = [d_1,d_2...,d_m]^T$ is the vector of the system discrete variables. $n$ and $m$ are the number of state variables and discrete variables of the system respectively. $v$ is the number of events and $zc$ is an event condition, $H$ and $L$ are functions used in the effects of the event $zc$.

For new users, take a look at the [Tutorial](https://mongibellili.github.io/QuantizedSystemSolver.jl/dev/guide/userTutorial/) section. If you see something wrong,
please open an [issue](https://github.com/mongibellili/QuantizedSystemSolver.jl/issues)

For developpers, take a look at the [Developer Guide](https://mongibellili.github.io/QuantizedSystemSolver.jl/dev/developer/devIntro/) section. Then, if you have an idea,
do a [pull request](https://github.com/mongibellili/QuantizedSystemSolver.jl/pulls)!





