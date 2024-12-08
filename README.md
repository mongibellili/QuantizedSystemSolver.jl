# <img width="50" height="50" style="position:relative; top:15px" src="docs\src\logo.png"> QuantizedSystemSolver
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mongibellili.github.io/QuantizedSystemSolver.jl/dev/)
[![CI](https://github.com/mongibellili/QuantizedSystemSolver/actions/workflows/CI.yml/badge.svg)](https://github.com/mongibellili/QuantizedSystemSolver/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/mongibellili/QuantizedSystemSolver/branch/main/graph/badge.svg)](https://codecov.io/gh/mongibellili/QuantizedSystemSolver)
<a style="float: right" href="https://mongibellili.github.io/QuantizedSystemSolver.jl/dev/">Documentation</a>

The growing intricacy of contemporary engineering systems, typically reduced to differential equations with events, poses a difficulty in digitally simulating them using traditional numerical integration techniques. The Quantized State System (QSS) and the Linearly Implicit Quantized State System (LIQSS) are different methods for tackling such problems.
The QuantizedSystemSolver aims to solve a set of Ordinary differential equations with a set of events. It implements the quantized state system methods: An approach that builds the solution by updating the system variables independently as opposed to classic integration methods that update all the system variables every step.

#   <span style="color:red">Installation</span>
QuantizedSystemSolver.jl is a [registered package](http://pkg.julialang.org), and is
simply installed by the foloowing:

Run julia in the terminal, then enter ] to bring up Julia's package manager, and add the QuantizedSystemSolver.jl package:
```console
julia

julia> ]
 
(@v1.x) pkg> add QuantizedSystemSolver
```

In order to solve any problem using the quantizedSystemSolver, we have to construct it in the follwing form:

$\dot X=f(X,P,t) $

$ if \; zc_v(x_i...,p_d...,t) \; i \in [1,n]  \;  ; \; d  \in [1,m] $

$ \qquad x_i=H_v(x_i...,p_d...,t) $

$ \qquad p_d=L_v(x_i...,p_d...,t)  $

$ \qquad \qquad...$

where $X=[x_1,x_2...,x_n]^T$ and $P=[p_1,p_2...,p_m]^T$ are the state variables and discrete variables of the system respectively. $v$ is the number of events and $zc$ is an event condition, $H$ and $L$ are functions used in the effects of the event $zc$.


For new users, take a look at the [Tutorial](@ref) section. If you see something wrong,
please open an [issue](https://github.com/mongibellili/QuantizedSystemSolver.jl/issues)

For developpers, take a look at the [Developer Guide](@ref) section. The, if you have an idea,
do a [pull request](https://github.com/mongibellili/QuantizedSystemSolver.jl/pulls)!





