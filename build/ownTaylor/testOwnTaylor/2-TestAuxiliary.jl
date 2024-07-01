using qss
using BenchmarkTools
t2=Taylor0([2.0,2.0,3.0],2)
t1=Taylor0([1.0,1.1],1)
cache1=Taylor0([0.0,0.0,0.0],2)
cache2=Taylor0([0.0,0.0,0.0],2)

function test0(t2::Taylor0{Float64},t1::Taylor0{Float64})
    qss.resize_coeffs1!(t2.coeffs,1)
    #@show 
    t2.coeffs[1]
   # @show 
    qss.getcoeff(t2,1)
   # @show 
    qss.getindex(t2,1)
   # @show 
    qss.setindex!(t2,2.22,1)
    return nothing
end

function test1(t2::Taylor0{Float64},t1::Taylor0{Float64})
     qss.iterate(t1,0)
      qss.iterate(t1,1)
end

function test2(t2::Taylor0{Float64},t1::Taylor0{Float64},cache1::Taylor0{Float64},cache2::Taylor0{Float64})
    qss.fixorder(t2,t1)
end

#= @show test0(t2,t1)
@show t2 
@btime test0(t2,t1)=#

#@btime test1(t2,t1)
#test1(t2,t1)
#=  
t1,t2=test2(t1,t2,cache1,cache2)
@show t1
@show t2 =#
#@btime test2(t2,t1,cache1,cache2)


function test4(t2::Taylor0{Float64},t1::Taylor0{Float64},cache1::Taylor0{Float64},cache2::Taylor0{Float64})
      qss.copyto!(cache1,t2)
end

#= @show test4(t2,t1,cache1,cache2)
@show cache1 =#
@btime test4(t2,t1,cache1,cache2)