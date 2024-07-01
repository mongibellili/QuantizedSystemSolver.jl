using qssgenerated
using BenchmarkTools
t2=Taylor0([2.0,2.0,3.0],2)
t1=Taylor0([1.0,1.0,0.0],2)
cache1=Taylor0([0.0,0.0,0.0],2)
cache2=Taylor0([0.0,0.0,0.0],2)

function test0(t2::Taylor0{Float64})
    one(t2)
end
function test0(t2::Taylor0{Float64},cache::Taylor0{Float64})
    qssgenerated.one(t2,cache)
end
function test1(t2::Taylor0{Float64},t1::Taylor0{Float64})
    res=t2/t1/t1
end

function test2(t2::Taylor0{Float64},t1::Taylor0{Float64},cache1::Taylor0{Float64},cache2::Taylor0{Float64})
    res2=qssgenerated.divT(qss.divT(t2,t1,cache1),t1,cache2)
end

#= @show test0(t2)
@show test0(t2,cache1) =#
#@btime test0(t2)
#@btime test0(t2,cache1)

#= @show test1(t2,t1)
@show test2(t2,t1,cache1,cache2) =#
#@btime test1(t2,t1)
#@btime test2(t2,t1,cache1,cache2)
