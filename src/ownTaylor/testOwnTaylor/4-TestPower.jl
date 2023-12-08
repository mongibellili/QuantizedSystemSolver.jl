using qss
using BenchmarkTools
t2=Taylor0([0.5,2.0,3.0],2)
t1=Taylor0([1.0,1.1],1)
x=Taylor0([0.0,0.0,0.0],2)
cache1=Taylor0([0.0,0.0,0.0],2)
cache2=Taylor0([0.0,0.0,0.0],2)
cache3=Taylor0([0.0,0.0,0.0],2)
function testsquare0(t2::Taylor0{Float64},cache1::Taylor0{Float64})
    qss.square(t2,cache1)
end
function testsquare1(t2::Taylor0{Float64},t1::Taylor0{Float64})
     t2^2
end
function testsqrt0(t2::Taylor0{Float64},cache1::Taylor0{Float64})
    qss.sqrt(t2,cache1)
end
function testsqrt1(t2::Taylor0{Float64},t1::Taylor0{Float64})
     sqrt(t2)
end

function testpower0(t2::Taylor0{Float64},cache1::Taylor0{Float64},cache2::Taylor0{Float64})
    qss.powerT(t2,3,cache1)
end
function testpower1(t2::Taylor0{Float64},cache1::Taylor0{Float64},cache2::Taylor0{Float64})
    ^(t2,3)
end
#################real power###############
function testrealpower0(t2::Taylor0{Float64},cache1::Taylor0{Float64})
    qss.powerT(t2,3.22,cache1)
end
function testrealpower1(t2::Taylor0{Float64},cache1::Taylor0{Float64})
    ^(t2,3.22)
end


#= 
@show testsquare0(t2,cache1)
@show testsquare1(t2,cache1) =#

#@btime testsquare0(t2,cache1)#21.193 ns (0 allocations: 0 bytes)
#@btime testsquare1(t2,cache1)# 102.524 ns (2 allocations: 160 bytes)


#= @show testsqrt0(t2,cache1)
@show testsqrt1(t2,cache1) =#
#@btime testsqrt0(t2,cache1)#20.749 ns (0 allocations: 0 bytes)
#@btime testsqrt1(t2,cache1)#99.014 ns (2 allocations: 160 bytes)

#= @show testpower0(t2,cache1,cache2)
@show testpower1(t2,cache1,cache2) =#
#@btime testpower0(t2,cache1,cache2)#60.724 ns (0 allocations: 0 bytes)
#@btime testpower1(t2,cache1,cache2)#132.472 ns (2 allocations: 160 bytes)


#= @show testrealpower0(t2,cache1)
@show testrealpower1(t2,cache1) =#
#@btime testrealpower0(t2,cache1)#142.59 ns (0 allocations: 0 bytes)
#@btime testrealpower1(t2,cache1)#216.271 ns (2 allocations: 160 bytes)




