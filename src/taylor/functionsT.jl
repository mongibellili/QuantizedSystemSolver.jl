# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#
# Functions

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
    order = a.order
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
    @inbounds c[k] = (a[k] - c[k] / k) / constant_term(r)
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
    @inbounds c[k] = -(a[k] + c[k] / k) / constant_term(r)
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
