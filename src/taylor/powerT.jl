# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#
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