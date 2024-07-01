using qss
using BenchmarkTools
t2=Taylor0([0.5,2.0,3.0],2)
t1=Taylor0([1.0,1.1],1)
x=Taylor0([0.0,0.0,0.0],2)
cache1=Taylor0([0.0,0.0,0.0],2)
cache2=Taylor0([0.0,0.0,0.0],2)
cache3=Taylor0([0.0,0.0,0.0],2)
function testexp0(t2::Taylor0{Float64},cache1::Taylor0{Float64})
    qss.exp(t2,cache1)
end
function testexp1(t2::Taylor0{Float64},t1::Taylor0{Float64})
     exp(t2)
end
function testlog0(t2::Taylor0{Float64},cache1::Taylor0{Float64})
    qss.log(t2,cache1)
end
function testlog1(t2::Taylor0{Float64},t1::Taylor0{Float64})
     log(t2)
end
function testsincos0(t2::Taylor0{Float64},cache1::Taylor0{Float64},cache2::Taylor0{Float64})
    qss.sincos(t2,cache1,cache2)
end
function testsincos1(t2::Taylor0{Float64},cache1::Taylor0{Float64},cache2::Taylor0{Float64})
    sincos(t2)
end
function testtan0(t2::Taylor0{Float64},cache1::Taylor0{Float64},cache2::Taylor0{Float64})
    qss.tan(t2,cache1,cache2)
end
function testtan1(t2::Taylor0{Float64},cache1::Taylor0{Float64},cache2::Taylor0{Float64})
    tan(t2)
end
function testasin0(t2::Taylor0{Float64},cache1::Taylor0{Float64},cache2::Taylor0{Float64},cache3::Taylor0{Float64})
    qss.asin(t2,cache1,cache2,cache3)
end
function testasin1(t2::Taylor0{Float64},cache1::Taylor0{Float64},cache2::Taylor0{Float64})
   asin(t2)
end
function testacos0(t2::Taylor0{Float64},cache1::Taylor0{Float64},cache2::Taylor0{Float64},cache3::Taylor0{Float64})
    qss.acos(t2,cache1,cache2,cache3)
end
function testacos1(t2::Taylor0{Float64},cache1::Taylor0{Float64},cache2::Taylor0{Float64})
   acos(t2)
end

function testabs0(t2::Taylor0{Float64},cache1::Taylor0{Float64})
    qss.abs(t2,cache1)
 end
 function testabs1(t2::Taylor0{Float64},cache1::Taylor0{Float64})
    abs(t2)
 end


#= @show testlog0(t2,cache1)
@show testlog1(t2,cache1) =#
#@btime testlog0(t2,cache1)
#@btime testlog1(t2,cache1)

#= @show testsincos0(t2,cache1,cache2)
@show testsincos1(t2,cache1,cache2) =#
#@btime testsincos0(t2,cache1,cache2)
#@btime testsincos1(t2,cache1,cache2)

#= @show testtan0(t2,cache1,cache2)
@show testtan1(t2,cache1,cache2) =#
#@btime testtan0(t2,cache1,cache2)
#@btime testtan1(t2,cache1,cache2)

#= @show testasin0(t2,cache1,cache2,cache3)
@show testasin1(t2,cache1,cache2) =#
#@btime testasin0(t2,cache1,cache2,cache3)
#@btime testasin1(t2,cache1,cache2)
#= 
@show testacos0(t2,cache1,cache2,cache3)
@show testacos1(t2,cache1,cache2) =#
#@btime testasin0(t2,cache1,cache2,cache3) #(0 allocations: 0 bytes)
#@btime testasin1(t2,cache1,cache2) #(15 allocations: 1.08 KiB)


#= @show testabs0(t2,cache1)
@show testabs1(t2,cache1) =#
@btime testabs0(t2,cache1)
#@btime testabs1(t2,cache1)

