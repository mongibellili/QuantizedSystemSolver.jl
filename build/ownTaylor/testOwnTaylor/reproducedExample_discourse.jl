
using BenchmarkTools
import Base:   eachindex, getindex, setindex!
struct MyType{T<:Number} 
    coeffs :: Array{T,1}
    size :: Int
end
getindex(a::MyType, n::Int) = a.coeffs[n+1]
setindex!(a::MyType{T}, x::T, n::Int) where {T<:Number} = a.coeffs[n+1] = x
@inline firstindex(a::MyType) = 0
@inline lastindex(a::MyType) = a.size
@inline eachindex(a::MyType) = firstindex(a):lastindex(a)
function mulTfoldl(res::P, bs...) where {P<:MyType}
    l = length(bs)
    i =  2;  l == i && return res 
    i =  3;  l == i && return              mulTT(res, bs[1],bs[end-1],bs[end])
    i =  4; l == i && return           mulTT(mulTT(res, bs[1],bs[end-1],bs[end]),bs[2],bs[end-1],bs[end])
    i =  5;  l == i && return        mulTT(mulTT(mulTT(res, bs[1],bs[end-1],bs[end]),bs[2],bs[end-1],bs[end]),bs[3],bs[end-1],bs[end])
    i =  6;  l == i && return      mulTT(mulTT(mulTT(mulTT(res, bs[1],bs[end-1],bs[end]),bs[2],bs[end-1],bs[end]),bs[3],bs[end-1],bs[end]),bs[4],bs[end-1],bs[end])
end
function mulTT(a::P, b::R, c::Q, xs...)where {P,Q,R <:Union{MyType,Number}}
    if length(xs)>1
      mulTfoldl( mulTT( mulTT(a,b,xs[end-1],xs[end]) , c,xs[end-1] ,xs[end] ), xs...)
    end
end
function mulTT(a::MyType{T}, b::T,cache1::MyType{T},cache2::MyType{T}) where {T<:Number}
  fill!(cache2.coeffs, b)
  @__dot__ cache1.coeffs = a.coeffs * cache2.coeffs  ##fixed broadcast dimension mismatch
  return cache1
end
mulTT(a::T,b::MyType{T}, cache1::MyType{T},cache2::MyType{T}) where {T<:Number} = mulTT(b , a,cache1,cache2)
function mulTT(a::MyType{T}, b::MyType{T},cache1::MyType{T},cache2::MyType{T}) where {T<:Number}
 for k in eachindex(a)
      @inbounds cache2[k] = a[0] * b[k]
      @inbounds for i = 1:k
        cache2[k] += a[i] * b[k-i]
      end
  end
  @__dot__ cache1.coeffs = cache2.coeffs
  return cache1
end
function mulTT(a::T, b::T,cache1::MyType{T},cache2::MyType{T}) where {T<:Number}
  cache1[0]=a*b  
  return cache1
end

av = vec(3(rand(1,3) .- 1.5)) 
bv = vec(3(rand(1,3) .- 1.5))
cv = vec(3(rand(1,3) .- 1.5))
dv = vec(3(rand(1,3) .- 1.5))

a0 =  MyType(av,2)
b0 =  MyType(bv,2)
c0 =  MyType(cv,2)
d0 =  MyType(dv,2)
e=3rand() - 1.5
f=3rand() - 1.5
cache=[MyType([0.0,0.0,0.0],2),MyType([0.0,0.0,0.0],2)]
function boundaryTest(a0::MyType{Float64},b0::MyType{Float64},c0::MyType{Float64},d0::MyType{Float64},cache::Vector{MyType{Float64}},e::Float64,f::Float64)
    @show mulTT(e, f, a0, b0, e, c0,d0, cache[1], cache[2])
end
boundaryTest(a0,b0,c0,d0,cache,e,f)