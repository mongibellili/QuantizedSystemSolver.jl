# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
# Modified by elmongi elbellili, 2024; cache added to avoid allocations
# FunctionsT

"""
    exp(a::Taylor0, c::Taylor0)

Calculates the exponential of `a` and stores the result in `c`.

# Arguments
- `a::Taylor0`: The input Taylor series.
- `c::Taylor0`: The cache to store the result.

# Returns
- `c::Taylor0`: The result of `exp(a)`.
"""
function exp(a::Taylor0, c::Taylor0)
    #order = a.order
    #aux = exp(constant_term(a))
    for k in eachindex(a)
        exp!(c, a, k) 
    end
    return c
end

"""
    exp(a::T, c::Taylor0) where {T<:Number}

Calculates the exponential of `a` and stores the result in `c`.

# Arguments
- `a::T`: The input number.
- `c::Taylor0`: The cache to store the result.

# Returns
- `c::Taylor0`: The result of `exp(a)`.
"""
function exp(a::T, c::Taylor0) where {T<:Number}
    c[0] = exp(a)
    return c
end

"""
    log(a::Taylor0, c::Taylor0)

Compute the logarithm of a `Taylor0` object `a` and stores the result in `Taylor0` object `c`.

# Arguments
- `a::Taylor0`: The `Taylor0` object for which the logarithm is to be computed.
- `c::Taylor0`: to store the result

# Returns
- object c .

"""
function log(a::Taylor0, c::Taylor0)
    iszero(constant_term(a)) && throw(DomainError(a,
        """The 0-th order coefficient must be non-zero in order to expand `log` around 0."""))
    for k in eachindex(a)
        log!(c, a, k)
    end
    return c
end
function log(a::T, c::Taylor0) where {T<:Number}
    iszero(a) && throw(DomainError(a,
        """ log(0) undefined."""))
    c[0] = log(a)
    return c
end

function sincos(a::Taylor0, s::Taylor0, c::Taylor0)
    for k in eachindex(a)
        sincos!(s, c, a, k)
    end
    return s, c
end
sin(a::Taylor0, s::Taylor0, c::Taylor0) = sincos(a, s, c)[1]
cos(a::Taylor0, c::Taylor0, s::Taylor0) = sincos(a, s, c)[2]
"""
    sin(a::T, s::Taylor0, c::Taylor0) where {T<:Number}

Calculates the sine of the number `a` and stores the result in the Taylor series `s`.

# Arguments
- `a::T`: The input number.
- `s::Taylor0`: The Taylor series to store the result.
- `c::Taylor0`: An auxiliary Taylor series (not used in this function).

# Returns
- `s::Taylor0`: The result of `sin(a)`.
"""
function sin(a::T, s::Taylor0, c::Taylor0) where {T<:Number}
    s[0] = sin(a)
    return s
end
"""
    cos(a::T, s::Taylor0, c::Taylor0) where {T<:Number}

Calculates the cosine of the number `a` and stores the result in the Taylor series `s`.

# Arguments
- `a::T`: The input number.
- `s::Taylor0`: The Taylor series to store the result.
- `c::Taylor0`: An auxiliary Taylor series (not used in this function).

# Returns
- `s::Taylor0`: The result of `cos(a)`.
"""
function cos(a::T, s::Taylor0, c::Taylor0) where {T<:Number}
    s[0] = cos(a)
    return s
end
"""
    tan(a::Taylor0, c::Taylor0, c2::Taylor0)

Calculates the tangent of the Taylor series `a` and stores the result in `c`.

# Arguments
- `a::Taylor0`: The input Taylor series.
- `c::Taylor0`: The Taylor series to store the result.
- `c2::Taylor0`: An auxiliary Taylor series used in the calculation.

# Returns
- `c::Taylor0`: The result of `tan(a)`.
"""
function tan(a::Taylor0, c::Taylor0, c2::Taylor0)
    for k in eachindex(a)
        tan!(c, a, c2, k)
    end
    return c
end
function tan(a::T, c::Taylor0, c2::Taylor0) where {T<:Number}
    c[0] = tan(a)
    return c
end

"""
    asin(a::Taylor0, c::Taylor0, r::Taylor0, cache3::Taylor0)

Calculates the arcsine of `a` and stores the result in `c`.

# Arguments
- `a::Taylor0`: The input Taylor series.
- `c::Taylor0`: The cache to store the result.
- `r::Taylor0`: An auxiliary Taylor series.
- `cache3::Taylor0`: Another auxiliary Taylor series.

# Returns
- `c::Taylor0`: The result of `acos(a)`.
"""
function asin(a::Taylor0, c::Taylor0, r::Taylor0, cache3::Taylor0)
    a0 = constant_term(a)
    a0^2 == one(a0) && throw(DomainError(a,
        """Series expansion of asin(x) diverges at x = ±1."""))
    for k in eachindex(a)
        asin!(c, a, r, cache3, k)
    end
    return c
end

"""
    acos(a::Taylor0, c::Taylor0, r::Taylor0, cache3::Taylor0)

Calculates the arccosine of `a` and stores the result in `c`.

# Arguments
- `a::Taylor0`: The input Taylor series.
- `c::Taylor0`: The cache to store the result.
- `r::Taylor0`: An auxiliary Taylor series.
- `cache3::Taylor0`: Another auxiliary Taylor series.

# Returns
- `c::Taylor0`: The result of `acos(a)`.
"""
function acos(a::Taylor0, c::Taylor0, r::Taylor0, cache3::Taylor0)
    a0 = constant_term(a)
    a0^2 == one(a0) && throw(DomainError(a,
        """Series expansion of asin(x) diverges at x = ±1."""))
    for k in eachindex(a)
        acos!(c, a, r, cache3, k)
    end
    return c
end

"""
    atan(a::Taylor0, c::Taylor0, r::Taylor0)

Calculates the arctangent of `a` and stores the result in `c`.

# Arguments
- `a::Taylor0`: The input Taylor series.
- `c::Taylor0`: The cache to store the result.
- `r::Taylor0`: An auxiliary Taylor series.

# Returns
- `c::Taylor0`: The result of `atan(a)`.
"""
function atan(a::Taylor0, c::Taylor0, r::Taylor0)
    for k in eachindex(a)
        atan!(c, a, r, k)
    end
    return c
end

function atan2(a::Taylor0, c::Taylor0, r::Taylor0)
    r[0]=Base.atan2(a[0],c[0])
    return r
end

"""
    asin!(c::Taylor0, a::Taylor0, r::Taylor0, cache3::Taylor0, k::Int)

Compute the arcsine of a Taylor series `a` and store the result in `c`. 
The computation uses intermediate results stored in `r` and `cache3`. 
The parameter `k` specifies the order of the Taylor series expansion. 
This is a copy of the function `asin!` as written in the package
TaylorSeries.jl but with an added cache to avoid allocation.

# Arguments
- `c::Taylor0`: The Taylor series where the result will be stored.
- `a::Taylor0`: The input Taylor series for which the arcsine is computed.
- `r::Taylor0`: An intermediate Taylor series used in the computation.
- `cache3::Taylor0`: Another intermediate Taylor series used in the computation to avoid one allocation of function square inside.
- `k::Int`: The order of the Taylor series expansion.

# Returns
- nothing: The result of the arcsine computation stored in `c`.

"""
@inline function asin!(c::Taylor0, a::Taylor0, r::Taylor0, cache3::Taylor0, k::Int)
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = asin(a0)
        @inbounds r[0] = sqrt(1 - a0^2)
        return nothing
    end
    @inbounds c[k] = (k - 1) * r[1] * c[k-1]
    @inbounds for i in 2:k-1
        c[k] += (k - i) * r[i] * c[k-i]
    end
    cache3 = square(a, cache3)
    @__dot__ cache3.coeffs = (-)(cache3.coeffs)
    cache3[0] = 1 + cache3[0]  #1-square(a,cache3)=1-a^2
    sqrt!(r, cache3, k)
    r0=constant_term(r)
    if r0==0.0
        r0=1e-9
    end
    @inbounds c[k] = (a[k] - c[k] / k) / r0
    return nothing
end
"""
        acos!(c::Taylor0, a::Taylor0, r::Taylor0, cache3::Taylor0, k::Int)

Compute the arccosine of a Taylor series `a` and store the result in `c`. 
The computation uses intermediate results stored in `r` and `cache3`. 
The parameter `k` specifies the order of the Taylor series expansion.
This is a copy of the function `acos!` as written in the package
TaylorSeries.jl but with an added cache to avoid one allocation.

# Arguments
- `c::Taylor0`: The Taylor series where the result will be stored.
- `a::Taylor0`: The input Taylor series for which the arcsine is computed.
- `r::Taylor0`: An intermediate Taylor series used in the computation.
- `cache3::Taylor0`: Another intermediate Taylor series used in the computation to avoid allocation of function square inside.
- `k::Int`: The order of the Taylor series expansion.

# Returns
- nothing: The result of the arcsine computation stored in `c`.

"""
@inline function acos!(c::Taylor0, a::Taylor0, r::Taylor0, cache3::Taylor0, k::Int)
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = acos(a0)
        @inbounds r[0] = sqrt(1 - a0^2)
        return nothing
    end
    @inbounds c[k] = (k - 1) * r[1] * c[k-1]
    @inbounds for i in 2:k-1
        c[k] += (k - i) * r[i] * c[k-i]
    end
    cache3 = square(a, cache3)
    @__dot__ cache3.coeffs = (-)(cache3.coeffs)
    cache3[0] = 1 + cache3[0]  #1-square(a,cache3)=1-a^2
    sqrt!(r, cache3, k)
    r0=constant_term(r)
    if r0==0.0
        r0=1e-9
    end
    @inbounds c[k] = -(a[k] + c[k] / k) / r0
    return nothing
end

"""
    abs(a::Taylor0, cache1::Taylor0)

Compute the absolute value of a `Taylor0` object `a` and store the result in `cache1`.

# Arguments
- `a::Taylor0`: The input `Taylor0` object whose absolute value is to be computed.
- `cache1::Taylor0`: A `Taylor0` object to store the result of the absolute value computation.

# Returns
- `cache1::Taylor0`: The `Taylor0` object containing the absolute value of `a`.

"""
function abs(a::Taylor0, cache1::Taylor0)
    if constant_term(a) > 0
        @__dot__ cache1.coeffs = (a.coeffs)
        return cache1
    elseif constant_term(a) < 0
        @__dot__ cache1.coeffs = (-)(a.coeffs)
        return cache1
    else
        cache1.coeffs .= Inf # no need to throw error, Inf is fine...for my solver i deal with it by guarding against small steps
        cache1[0] = 0.0
        return cache1
    end
end
function abs(a::T, cache1::Taylor0) where {T<:Number}
    cache1[0] = abs(a)
    return cache1
end


#powers

powerT(a::Taylor0, n::Integer, cache1::Taylor0) = powerT(a, float(n), cache1)
## Real power ##
function powerT(a::Taylor0, r::S, cache1::Taylor0) where {S<:Real}
    if iszero(r) #later test r=0.0
        cache1[0] = 1.0
        return cache1
    end
    r == 1.0 && return a
    r == 2.0 && return square(a, cache1)
    r == 1 / 2 && return sqrt(a, cache1)
    l0 = findfirst(a)
    lnull = trunc(Int, r * l0)
    if (a.order - lnull < 0) || (lnull > a.order)
        return cache1  #empty
    end
    c_order = l0 == 0 ? a.order : min(a.order, trunc(Int, r * a.order))
    for k = 0:c_order
        pow!(cache1, a, r, k)
    end
    return cache1
end

"""
    powerT(a::T, r::S, cache1::Taylor0) where {S<:Real,T<:Number}

Raises `a` to the power `r` and stores the result in `cache1`.

# Arguments
- `a::T`: The base, a number.
- `r::S`: The exponent, a real number.
- `cache1::Taylor0`: The cache to store the result.

# Returns
- `cache1::Taylor0`: The result of `a^r`.
"""
function powerT(a::T, r::S, cache1::Taylor0) where {S<:Real,T<:Number}
    cache1[0] = a^r
    return cache1
end

"""
    square(a::Taylor0, cache1::Taylor0)

Calculates the square of `a` and stores the result in `cache1`.

# Arguments
- `a::Taylor0`: The input Taylor series.
- `cache1::Taylor0`: The cache to store the result.

# Returns
- `cache1::Taylor0`: The result of `a^2`.
"""
function square(a::Taylor0, cache1::Taylor0)
    cache1[0] = constant_term(a)^2
    for k in 1:a.order
        sqr!(cache1, a, k)
    end
    return cache1
end

## Square root ##
"""
    sqrt(a::Taylor0, cache1::Taylor0)

Calculates the square root of `a` and stores the result in `cache1`.

# Arguments
- `a::Taylor0`: The input Taylor series.
- `cache1::Taylor0`: The cache to store the result.

# Returns
- `cache1::Taylor0`: The result of `sqrt(a)`.
"""
function sqrt(a::Taylor0, cache1::Taylor0)
   
    l0nz = findfirst(a)
    if l0nz < 0
        return cache1
    end
    cache1[0] = sqrt(a[0])
    for k = 1:a.order
        sqrt!(cache1, a, k, 0)
    end
    #@show cache1
    return cache1
end

"""
    sqrt(a::T, cache1::Taylor0) where {T<:Number}

Calculates the square root of `a` and stores the result in `cache1`.

# Arguments
- `a::T`: The input number.
- `cache1::Taylor0`: The cache to store the result.

# Returns
- `cache1::Taylor0`: The result of `sqrt(a)`.
"""
function sqrt(a::T, cache1::Taylor0) where {T<:Number}
    cache1[0] = sqrt(a)
    return cache1
end


#= function LinearAlgebra.norm(v::AbstractVector{Taylor0})
    return norm([x.coeffs[1] for x in v])
end
 =#
#= function norm(v::AbstractVector{Taylor0})
    return norm([x.coeffs[1] for x in v])
end =#

#= function norm(v::AbstractVector{Taylor0})
    return sqrt(sum(x.coeffs[1]^2 for x in v)) 
end =#

function rad2deg(x::Taylor0) 
    return x[0] * 180 / π
end