using TaylorSeries
using BenchmarkTools
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
function normalTest(a1::Taylor1{Float64},b1::Taylor1{Float64},c1::Taylor1{Float64},d1::Taylor1{Float64},e1::Taylor1{Float64},f1::Taylor1{Float64},g1::Taylor1{Float64},h1::Taylor1{Float64},a0::Taylor0{Float64},b0::Taylor0{Float64},c0::Taylor0{Float64},d0::Taylor0{Float64},e0::Taylor0{Float64},f0::Taylor0{Float64},g0::Taylor0{Float64},h0::Taylor0{Float64},cache::Vector{Taylor0{Float64}},i::Float64,j::Float64,k::Float64)
    +(-(+(+(*(*(-(a1,b1),j),c1),d1),i),e1),f1)
end
function freeTest(a1::Taylor1{Float64},b1::Taylor1{Float64},c1::Taylor1{Float64},d1::Taylor1{Float64},e1::Taylor1{Float64},f1::Taylor1{Float64},g1::Taylor1{Float64},h1::Taylor1{Float64},a0::Taylor0{Float64},b0::Taylor0{Float64},c0::Taylor0{Float64},d0::Taylor0{Float64},e0::Taylor0{Float64},f0::Taylor0{Float64},g0::Taylor0{Float64},h0::Taylor0{Float64},cache::Vector{Taylor0{Float64}},i::Float64,j::Float64,k::Float64)
  clearCache(cache)
  @changeAST (+(-(+(+(*(*(-(a0,b0),j),c0),d0),i),e0),f0),1)
end
@btime normalTest(a1,b1,c1,d1,e1,f1,g1,h1,a0,b0,c0,d0,e0,f0,g0,h0,cache,i,j,k)
@btime freeTest(a1,b1,c1,d1,e1,f1,g1,h1,a0,b0,c0,d0,e0,f0,g0,h0,cache,i,j,k)
