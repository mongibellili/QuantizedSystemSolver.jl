#= using Symbolics

 @variables t x y
 D = Differential(t)

z = t + t^2
@show  D(z) # symbolic representation of derivative(t + t^2, t)
 =#
 using QuantizedSystemSolver
 
 @show testfunc(1.0)
 