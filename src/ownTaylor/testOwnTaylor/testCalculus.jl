using qss
using StaticArrays
using BenchmarkTools
t2=Taylor0([2.0,2.0,3.0],2)
t1=Taylor0([1.0,1.1],1)
#= cache1=Taylor0([0.0,0.0,0.0],2)
cache2=Taylor0([0.0,0.0,0.0],2) =#
@show t1