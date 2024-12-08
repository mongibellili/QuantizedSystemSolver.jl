# QuantizedSystemSolver.jl

A [Julia](http://julialang.org) package for solving systems of ODEs with events.

---
### Statement of need

Traditional solvers are challenged by frequent discontinuities where the state of the system abruptly changes
at specific points or intervals. They struggle to accurately capture the dynamics around discontinuities especially in large sparse and stiff systems. They either undergo expensive iterations to pinpoint exact discontinuity instances or resort to interpolating their locations, resulting in unreliable outcomes. 
Written in the easy-to-learn [Julia programming language](https://julialang.org), 
and taking advantage of its features such as multiple dispatch and metaprogramming, the QuantizedSystemSolver.jl is a solver that aims to efficiently solve a set of Ordinary differential Equations with a set of events via implementing the QSS and LIQSS methods. It is the first such tool to be published in the Julia ecosystem.

### Authors

- Mongi Bellili, Belgian Royal Military Academy and Ku Leuven.



### License

QuantizedSystemSolver is licensed under the MIT license; see
[LICENSE](https://github.com/mongibellili/QuantizedSystemSolver.jl/blob/main/LICENSE) for
the full license text.

### Installation

QuantizedSystemSolver.jl is a [registered package](http://pkg.julialang.org), and is
simply installed by running

```console
julia

julia> ]

(@v1.x) pkg> add QuantizedSystemSolver
```





For new users, take a look at the [Tutorial](@ref) section. If you see something wrong,
please open an [issue](https://github.com/mongibellili/QuantizedSystemSolver.jl/issues)

For developpers, take a look at the [Developer Guide](@ref) section. The, if you have an idea,
do a [pull request](https://github.com/mongibellili/QuantizedSystemSolver.jl/pulls)!


### Related packages

- [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl): Differential Equations solver using classic integration methods.
- [QuantizedStateSystems.jl](https://github.com/hurak/QuantizedStateSystems.jl): julia qss-solver
- [qss-solver](https://github.com/CIFASIS/qss-solver): C qss-solver









