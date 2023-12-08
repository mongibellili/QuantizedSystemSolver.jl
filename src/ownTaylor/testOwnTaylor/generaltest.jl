using TaylorSeries
using BenchmarkTools
using qssgenerated
av = vec(3(rand(1,3) .+ 1.5)) #changed + for log and sqrt only
bv = vec(3(rand(1,3) .+ 1.5))
cv = vec(3(rand(1,3) .+ 1.5))
dv = vec(3(rand(1,3) .+ 1.5))
a1 =  Taylor1(av,2)
b1 =  Taylor1(bv,2)
c1 =  Taylor1(cv,2)
d1 =  Taylor1(dv,2)
a0 =  Taylor0(av,2)
b0 =  Taylor0(bv,2)
c0 =  Taylor0(cv,2)
d0 =  Taylor0(dv,2)
e=3rand() + 1.5
f=3rand() + 1.5
cache=[Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2)]

function clearCache(cache::Vector{Taylor0{Float64}}) #where {T<:Number}
    for i=1:10
      cache[i].coeffs.=0.0
    end
  end
function same(sum1::Taylor1{Float64},sum2::Taylor0{Float64}) 
    if sum1[0]==Inf && sum2[0]==Inf
        return true
    end
    for i=0:length(sum1.coeffs)-1
      if abs(sum1[i]-sum2[i])>1e-9
      #if (sum1[i]!=sum2[i])
        return false
      end
    end
    return true
end
function same(sum1::Float64,sum2::Taylor0{Float64}) 
    if sum1==Inf && sum2[0]==Inf
        return true
    end
    return abs(sum1-sum2[0])<1e-9
    #return (sum1==sum2[0])
  end
macro changeAST(ex)
  Base.remove_linenums!(ex)
  # dump(ex; maxdepth=18)
  qssgenerated.twoInOne(ex)
# dump( ex.args[1])
#@show ex.args[1]
#=  # return  
return nothing =#
esc(ex.args[1])# return 
end

function investigationTest(a1::Taylor1{Float64},b1::Taylor1{Float64},c1::Taylor1{Float64},d1::Taylor1{Float64},a0::Taylor0{Float64},b0::Taylor0{Float64},c0::Taylor0{Float64},d0::Taylor0{Float64},cache::Vector{Taylor0{Float64}},e::Float64,f::Float64)
    clearCache(cache)
    @show @changeAST((((b0 - f) - a0) * d0, 1))
    @show ((b1 - f) + a1) * d1
end

function outerInvestigation()
    for _=1:20
        av = vec(3(rand(1,3) .- 1.5))
        bv = vec(3(rand(1,3) .- 1.5))
        cv = vec(3(rand(1,3) .- 1.5))
        dv = vec(3(rand(1,3) .- 1.5))
        a1 =  Taylor1(av,2)
        b1 =  Taylor1(bv,2)
        c1 =  Taylor1(cv,2)
        d1 =  Taylor1(dv,2)
        a0 =  Taylor0(av,2)
        b0 =  Taylor0(bv,2)
        c0 =  Taylor0(cv,2)
        d0 =  Taylor0(dv,2)
        e=3rand() - 1.5
        f=3rand() - 1.5
        investigationTest(a1,b1,c1,d1,a0,b0,c0,d0,cache,e,f)
    end
end
#outerInvestigation()

function boundaryTest(a1::Taylor1{Float64},b1::Taylor1{Float64},c1::Taylor1{Float64},d1::Taylor1{Float64},a0::Taylor0{Float64},b0::Taylor0{Float64},c0::Taylor0{Float64},d0::Taylor0{Float64},cache::Vector{Taylor0{Float64}},e::Float64,f::Float64)

    mulTT(1.2, 2.1, a0, b0, 3.2, c0,d0, cache[1], cache[2])
   #@show 1.2*2.1*a1*b1*3.2*c1*d1
   return nothing
end
@btime boundaryTest(a1,b1,c1,d1,a0,b0,c0,d0,cache,e,f)
#= function stressTest(a1::Taylor1{Float64},b1::Taylor1{Float64},c1::Taylor1{Float64},d1::Taylor1{Float64},a0::Taylor0{Float64},b0::Taylor0{Float64},c0::Taylor0{Float64},d0::Taylor0{Float64},cache::Vector{Taylor0{Float64}},e::Float64,f::Float64)
    for op1 in (:+,:-,:*,:/)
    for op2 in (:+,:-,:*,:/)
    for op3 in (:+,:-,:*,:/)
        for (term11,term01) in((:a1,:a0),(:e,:e))
            for (term12,term02) in((:a1,:a0),(:b1,:b0))
                for (term13,term03) in((:c1,:c0),(:f,:f))
                    for (term14,term04) in((:d1,:d0),(:a1,:a0))
                    @eval begin
                        clearCache(cache)
                        @assert same($op3($op1($op2($term12,$term13),$term11),$term14), @changeAST ($op3($op1($op2($term02,$term03),$term01),$term04),1))
                       #=  if !same($op3($op1($op2($term12,$term13),$term11),$term14), @changeAST ($op3($op1($op2($term02,$term03),$term01),$term04),1))
                            @show cache
                            @show $op3($op1($op2($term12,$term13),$term11),$term14)
                            @show @changeAST ($op3($op1($op2($term02,$term03),$term01),$term04),1)
                        end =#
                    end 
                    end
                end
            end
        end
    end
    end
    end
end
=#
#stressTest(a1,b1,c1,d1,a0,b0,c0,d0,cache,e,f) #successful
function stressTest2(a1::Taylor1{Float64},b1::Taylor1{Float64},c1::Taylor1{Float64},d1::Taylor1{Float64},a0::Taylor0{Float64},b0::Taylor0{Float64},c0::Taylor0{Float64},d0::Taylor0{Float64},cache::Vector{Taylor0{Float64}},e::Float64,f::Float64)
    for op1 in (:abs,:exp,:cos,:sin)
    for op2 in (:abs,:exp,:cos,:sin)
    for op3 in (:abs,:exp,:cos,:sin)
        for (term11,term01) in((:a1,:a0),(:e,:e))
            for (term12,term02) in((:a1,:a0),(:b1,:b0))
                        @eval begin
                            clearCache(cache)
                            @assert same($op3($op1($op2($term12))), @changeAST ($op3($op1($op2($term02))),1))
                            clearCache(cache)
                            @assert same($op3($op1($op2($term11))), @changeAST ($op3($op1($op2($term01))),1))
                    end                               
            end
        end
    end
    end
    end
end
#stressTest2(a1,b1,c1,d1,a0,b0,c0,d0,cache,e,f)#successful

#= function stressTest3(a1::Taylor1{Float64},b1::Taylor1{Float64},c1::Taylor1{Float64},d1::Taylor1{Float64},a0::Taylor0{Float64},b0::Taylor0{Float64},c0::Taylor0{Float64},d0::Taylor0{Float64},cache::Vector{Taylor0{Float64}},e::Float64,f::Float64)
    for op1 in (:sqrt,:log)
    for op2 in (:sqrt,:log)

        for (term11,term01) in((:a1,:a0),(:e,:e))
            for (term12,term02) in((:a1,:a0),(:b1,:b0))
                    @eval begin
                    clearCache(cache)
                    @assert same($op1($op2($term12)), @changeAST ($op1($op2($term02)),1))
                    clearCache(cache)
                    @assert same($op1($op2($term11)), @changeAST ($op1($op2($term01)),1))
                    end                               
            end
        end

    end
    end
end =#
#stressTest3(a1,b1,c1,d1,a0,b0,c0,d0,cache,e,f)#successful