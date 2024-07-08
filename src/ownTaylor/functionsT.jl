# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#
# Functions
function exp(a::Taylor0, c::Taylor0)
    order = a.order
    #aux = exp(constant_term(a))
    for k in eachindex(a)
        exp!(c, a, k)
    end
    return c
end
function exp(a::T, c::Taylor0) where {T<:Number}
    c[0] = exp(a)
    return c
end
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
function sin(a::T, s::Taylor0, c::Taylor0) where {T<:Number}
    s[0] = sin(a)
    return s
end
function cos(a::T, s::Taylor0, c::Taylor0) where {T<:Number}
    s[0] = cos(a)
    return s
end
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
#####################################constant case not implemented yet
function asin(a::Taylor0, c::Taylor0, r::Taylor0, cache3::Taylor0)
    a0 = constant_term(a)
    a0^2 == one(a0) && throw(DomainError(a,
        """Series expansion of asin(x) diverges at x = ±1."""))
    for k in eachindex(a)
        asin!(c, a, r, cache3, k)
    end
    return c
end
function acos(a::Taylor0, c::Taylor0, r::Taylor0, cache3::Taylor0)
    a0 = constant_term(a)
    a0^2 == one(a0) && throw(DomainError(a,
        """Series expansion of asin(x) diverges at x = ±1."""))
    for k in eachindex(a)
        acos!(c, a, r, cache3, k)
    end
    return c
end
function atan(a::Taylor0, c::Taylor0, r::Taylor0)
    for k in eachindex(a)
        atan!(c, a, r, k)
    end
    return c
end
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
