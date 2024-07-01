# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

#= ## Differentiating ##
"""
    differentiate(a)

Return the `Taylor0` polynomial of the differential of `a::Taylor0`.
The order of the result is `a.order-1`.

The function `derivative` is an exact synonym of `differentiate`.
""" =#
function differentiate(a::Taylor0)
    res = Taylor0(zero(a[0]), a.order(a)-1)
    for ord in eachindex(res)
        differentiate!(res, a, ord)
    end
    return res
end

#= """
    derivative

An exact synonym of [`differentiate`](@ref).
""" =#
const derivative = differentiate

#= """
    differentiate!(res, a) --> nothing

In-place version of `differentiate`. Compute the `Taylor0` polynomial of the
differential of `a::Taylor0` and return it as `res` (order of `res` remains
unchanged).
""" =#
#= function differentiate!(p::Taylor0, a::Taylor0)
    for k in eachindex(a)
        #differentiate!(res, a, ord)
        @inbounds p[k] = (k+1)*a[k+1]
    end
    nothing
end =#

function differentiate!(p::Taylor0, a::Taylor0, k::Int)
    if k < a.order
        @inbounds p[k] = (k+1)*a[k+1]
    end
    return nothing
end
function differentiate!(::Val{O},cache::Taylor0, a::Taylor0)where {O}
    for k=0:O-1
       # differentiate!(res, a, ord)
       # if k < a.order
          @inbounds cache[k] = (k+1)*a[k+1]
      #end
    end
    cache[O]=0.0  #cache Letter O not zero
    nothing
end
function differentiate(a::Taylor0, n::Int) 
    if n > a.order
        return Taylor0(T, 0)
    elseif n==0
        return a
    else
        res = differentiate(a)
        for i = 2:n
            differentiate!(Val(a.order),res, res)
        end
        return Taylor0(view(res.coeffs, 1:a.order-n+1))
    end
end
function ndifferentiate!(cache::Taylor0,a::Taylor0, n::Int) 
    if n > a.order
        cache.coeffs.=0.0
    elseif n==0
        cache.coeffs.=a.coeffs
    else
        differentiate!(Val(a.order),cache,a)
        for i = 2:n
            differentiate!(Val(a.order),cache, cache)
        end
        #return Taylor0(view(res.coeffs, 1:a.order-n+1))
    end
end



#= """
    differentiate(n, a)

Return the value of the `n`-th differentiate of the polynomial `a`.
""" =#
function differentiate(n::Int, a::Taylor0) 
    @assert a.order ≥ n ≥ 0
    factorial( widen(n) ) * a[n] 
end

## Integrating ##
#= """
    integrate(a, [x])

Return the integral of `a::Taylor0`. The constant of integration
(0-th order coefficient) is set to `x`, which is zero if ommitted.
""" =#
function integrate(a::Taylor0, x::S) where {S<:Number}
    order = get_order(a)
    aa = a[0]/1 + zero(x)
    R = typeof(aa)
    coeffs = Array{typeof(aa)}(undef, order+1)
    fill!(coeffs, zero(aa))
    @inbounds for i = 1:order
        coeffs[i+1] = a[i-1] / i
    end
    @inbounds coeffs[1] = convert(R, x)
    return Taylor0(coeffs, a.order)
end
integrate(a::Taylor0)  = integrate(a, zero(a[0]))

function integrate!(p::Taylor0, a::Taylor0, x::S) where {S<:Number}
    p.coeffs[1]=x
    @inbounds for i = 1:a.order
        p[i] = a[i-1] / i
    end
    return nothing
end
function integrate!( a::Taylor0, x::S) where {S<:Number}
    
    @inbounds for i in a.order-1:-1:0
        a[i+1] = a[i] / (i+1)
    end
    a.coeffs[1]=x
    return nothing
end


