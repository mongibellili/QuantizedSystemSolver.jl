using TaylorSeries
using qssgenerated
import Base.:/
av = vec(3(rand(1,3) .- 1.5))
bv = vec(3(rand(1,3) .- 1.5))
cv = vec(3(rand(1,3) .- 1.5))
dv = vec(3(rand(1,3) .- 1.5))
ev = vec(3(rand(1,3) .- 1.5))
fv = vec(3(rand(1,3) .- 1.5))
gv = vec(3(rand(1,3) .- 1.5))
hv = vec(3(rand(1,3) .- 1.5))
a1 =  Taylor1(av,2)
b1 =  Taylor1(bv,2)
c1 =  Taylor1(cv,2)
d1 =  Taylor1(dv,2)
e1 =  Taylor1(ev,2)
f1 =  Taylor1(fv,2)
g1 =  Taylor1(gv,2)
h1 =  Taylor1(hv,2)
a0 =  Taylor0(av,2)
b0 =  Taylor0(bv,2)
c0 =  Taylor0(cv,2)
d0 =  Taylor0(dv,2)
e0 =  Taylor0(ev,2)
f0 =  Taylor0(fv,2)
g0 =  Taylor0(gv,2)
h0 =  Taylor0(hv,2)
i=3rand() - 1.5
j=3rand() - 1.5
k=3rand() - 1.5
cache=[Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2),Taylor0([0.0,0.0,0.0],2)]

function clearCache(cache::Vector{Taylor0{Float64}}) #where {T<:Number}
    for i=1:12
      cache[i].coeffs.=0.0
    end
  end
function same(sum1::Taylor1{Float64},sum2::Taylor0{Float64}) 
    isclose=true
    for i=0:length(sum1.coeffs)-1
      if abs(sum1[i]-sum2[i])>1e-9
      #if (sum1[i]!=sum2[i])
        isclose=false
        break
      end
    end
    return isclose
end
function same(sum1::Float64,sum2::Taylor0{Float64}) 
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

function /(a::T, b::Taylor1{T}) where {T<:Number}
    c=^(b, -1.0)
    c2 = Taylor1( b.order )
    @__dot__ c2.coeffs=c.coeffs*a ### broadcast dimension mismatch------------------div can have 2 caches then...:(
  
    return c2
  end
function stressTest(a1::Taylor1{Float64},b1::Taylor1{Float64},c1::Taylor1{Float64},d1::Taylor1{Float64},e1::Taylor1{Float64},f1::Taylor1{Float64},g1::Taylor1{Float64},h1::Taylor1{Float64},a0::Taylor0{Float64},b0::Taylor0{Float64},c0::Taylor0{Float64},d0::Taylor0{Float64},e0::Taylor0{Float64},f0::Taylor0{Float64},g0::Taylor0{Float64},h0::Taylor0{Float64},cache::Vector{Taylor0{Float64}},i::Float64,j::Float64,k::Float64)
    for op1 in (:+,:-,:*,:/)
    for op2 in (:+,:-,:*,:/)
    for op3 in (:+,:-,:*,:/)
        for (term11,term01) in((:a1,:a0),(:j,:j))
            for (term12,term02) in((:a1,:a0),(:b1,:b0))
                for (term13,term03) in((:c1,:c0),(:i,:i))
                    for (term14,term04) in((:d1,:d0),(:a1,:a0))
                        for (term15,term05) in((:a1,:a0),(:j,:j))
                            for (term16,term06) in((:a1,:a0),(:b1,:b0))
                                for (term17,term07) in((:c1,:c0),(:i,:i))
                                    for (term18,term08) in((:d1,:d0),(:a1,:a0))
                                            @eval begin
                                                clearCache(cache)
                                                @assert same($op1($op3($op2($op1($op3($op2($op1($term11,$term12),$term13),$term14),$term15),$term16),$term17),$term18),@changeAST ($op1($op3($op2($op1($op3($op2($op1($term01,$term02),$term03),$term04),$term05),$term06),$term07),$term08),1))
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
    end
    end
    end
end

stressTest(a1,b1,c1,d1,e1,f1,g1,h1,a0,b0,c0,d0,e0,f0,g0,h0,cache,i,j,k)#successful

#= function stressTest2(a1::Taylor1{Float64},b1::Taylor1{Float64},c1::Taylor1{Float64},d1::Taylor1{Float64},a0::Taylor0{Float64},b0::Taylor0{Float64},c0::Taylor0{Float64},d0::Taylor0{Float64},cache::Vector{Taylor0{Float64}},e::Float64,f::Float64)
    for op1 in (:abs,:exp,:cos,:sin)
    for op2 in (:abs,:exp,:cos,:sin)
    for op3 in (:abs,:exp,:cos,:sin)
        for (term11,term01) in((:a1,:a0),(:e,:e))
            for (term12,term02) in((:a1,:a0),(:b1,:b0))
                    @eval begin
                    clearCache(cache)
                    @assert same($op3($op1($op2($term12))), @changeAST ($op3($op1($op2($term02))),1))
                    @assert same($op3($op1($op2($term11))), @changeAST ($op3($op1($op2($term01))),1))
                    end                               
            end
        end
    end
    end
    end
end
#stressTest2(a1,b1,c1,d1,a0,b0,c0,d0,cache,e,f)

function stressTest3(a1::Taylor1{Float64},b1::Taylor1{Float64},c1::Taylor1{Float64},d1::Taylor1{Float64},a0::Taylor0{Float64},b0::Taylor0{Float64},c0::Taylor0{Float64},d0::Taylor0{Float64},cache::Vector{Taylor0{Float64}},e::Float64,f::Float64)
    for op1 in (:abs,:exp,:cos,:sin)
    for op2 in (:abs,:exp,:cos,:sin)
    for op3 in (:abs,:exp,:cos,:sin)
        for (term11,term01) in((:a1,:a0),(:e,:e))
            for (term12,term02) in((:a1,:a0),(:b1,:b0))
                    @eval begin
                    clearCache(cache)
                    @assert same($op3($op1($op2($term12))), @changeAST ($op3($op1($op2($term02))),1))
                    @assert same($op3($op1($op2($term11))), @changeAST ($op3($op1($op2($term01))),1))
                    end                               
            end
        end
    end
    end
    end
end
#stressTest3(a1,b1,c1,d1,a0,b0,c0,d0,cache,e,f) =#