# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#
# Functions
function exp(a::Taylor0)
    order = a.order
    aux = exp(constant_term(a))
    aa = one(aux) * a
    c = Taylor0(aux, order)
    for k in eachindex(a)
        exp!(c, aa, k)
    end
    return c
end
function log(a::Taylor0)
    iszero(constant_term(a)) && throw(DomainError(a,
        """The 0-th order coefficient must be non-zero in order to expand `log` around 0."""))
    order = a.order
    aux = log(constant_term(a))
    aa = one(aux) * a
    c = Taylor0(aux, order)
    for k in eachindex(a)
        log!(c, aa, k)
    end
    return c
end
sin(a::Taylor0) = sincos(a)[1]
cos(a::Taylor0) = sincos(a)[2]
function sincos(a::Taylor0)
    order = a.order
    aux = sin(constant_term(a))
    aa = one(aux) * a
    s = Taylor0(aux, order)
    c = Taylor0(cos(constant_term(a)), order)
    for k in eachindex(a)
        sincos!(s, c, aa, k)
    end
    return s, c
end
function tan(a::Taylor0)
    order = a.order
    aux = tan(constant_term(a))
    aa = one(aux) * a
    c = Taylor0(aux, order)
    c2 = Taylor0(aux^2, order)
    for k in eachindex(a)
        tan!(c, aa, c2, k)
    end
    return c
end
function asin(a::Taylor0)
    a0 = constant_term(a)
    a0^2 == one(a0) && throw(DomainError(a,
        """Series expansion of asin(x) diverges at x = ±1."""))
    order = a.order
    aux = asin(a0)
    aa = one(aux) * a
    c = Taylor0(aux, order)
    r = Taylor0(sqrt(1 - a0^2), order)
    for k in eachindex(a)
        asin!(c, aa, r, k)
    end
    return c
end
function acos(a::Taylor0)
    a0 = constant_term(a)
    a0^2 == one(a0) && throw(DomainError(a,
        """Series expansion of acos(x) diverges at x = ±1."""))
    order = a.order
    aux = acos(a0)
    aa = one(aux) * a
    c = Taylor0(aux, order)
    r = Taylor0(sqrt(1 - a0^2), order)
    for k in eachindex(a)
        acos!(c, aa, r, k)
    end
    return c
end
function atan(a::Taylor0)
    order = a.order
    a0 = constant_term(a)
    aux = atan(a0)
    aa = one(aux) * a
    c = Taylor0(aux, order)
    r = Taylor0(one(aux) + a0^2, order)
    iszero(constant_term(r)) && throw(DomainError(a,
        """Series expansion of atan(x) diverges at x = ±im."""))
    for k in eachindex(a)
        atan!(c, aa, r, k)
    end
    return c
end
sinh(a::Taylor0) = sinhcosh(a)[1]
cosh(a::Taylor0) = sinhcosh(a)[2]
function sinhcosh(a::Taylor0)
    order = a.order
    aux = sinh(constant_term(a))
    aa = one(aux) * a
    s = Taylor0(aux, order)
    c = Taylor0(cosh(constant_term(a)), order)
    for k in eachindex(a)
        sinhcosh!(s, c, aa, k)
    end
    return s, c
end
function tanh(a::Taylor0)
    order = a.order
    aux = tanh(constant_term(a))
    aa = one(aux) * a
    c = Taylor0(aux, order)
    c2 = Taylor0(aux^2, order)
    for k in eachindex(a)
        tanh!(c, aa, c2, k)
    end
    return c
end
function asinh(a::Taylor0)
    order = a.order
    a0 = constant_term(a)
    aux = asinh(a0)
    aa = one(aux) * a
    c = Taylor0(aux, order)
    r = Taylor0(sqrt(a0^2 + 1), order)
    iszero(constant_term(r)) && throw(DomainError(a,
        """Series expansion of asinh(x) diverges at x = ±im."""))
    for k in eachindex(a)
        asinh!(c, aa, r, k)
    end
    return c
end
function acosh(a::Taylor0)
    a0 = constant_term(a)
    a0^2 == one(a0) && throw(DomainError(a,
        """Series expansion of acosh(x) diverges at x = ±1."""))
    order = a.order
    aux = acosh(a0)
    aa = one(aux) * a
    c = Taylor0(aux, order)
    r = Taylor0(sqrt(a0^2 - 1), order)
    for k in eachindex(a)
        acosh!(c, aa, r, k)
    end
    return c
end
function atanh(a::Taylor0)
    order = a.order
    a0 = constant_term(a)
    aux = atanh(a0)
    aa = one(aux) * a
    c = Taylor0(aux, order)
    r = Taylor0(one(aux) - a0^2, order)
    iszero(constant_term(r)) && throw(DomainError(a,
        """Series expansion of atanh(x) diverges at x = ±1."""))
    for k in eachindex(a)
        atanh!(c, aa, r, k)
    end
    return c
end



@inline function exp!(c::Taylor0, a::Taylor0, k::Int)
    if k == 0
        @inbounds c[0] = exp(constant_term(a))  #this was already computed before!!no need to do anything if k==0
        return nothing
    end
    @inbounds c[k] = k * a[k] * c[0]
    @inbounds for i = 1:k-1
        c[k] += (k - i) * a[k-i] * c[i]
    end
    @inbounds c[k] = c[k] / k
    return nothing
end

@inline function log!(c::Taylor0, a::Taylor0, k::Int)
    if k == 0
        @inbounds c[0] = log(constant_term(a))
        return nothing
    elseif k == 1
        a0 = constant_term(a)
        if a0==0.0
            a0 = 1e-9
        end
        @inbounds c[1] = a[1] / a0
        return nothing
    end
    @inbounds c[k] = (k - 1) * a[1] * c[k-1]
    @inbounds for i = 2:k-1
        c[k] += (k - i) * a[i] * c[k-i]
    end
    a0 = constant_term(a)
    if a0==0.0
        a0 = 1e-9
    end
    @inbounds c[k] = (a[k] - c[k] / k) / a0
    return nothing
end

@inline function sincos!(s::Taylor0, c::Taylor0, a::Taylor0, k::Int)
    if k == 0
        a0 = constant_term(a)
        @inbounds s[0], c[0] = sincos(a0)
        return nothing
    end
    x = a[1]
    @inbounds s[k] = x * c[k-1]
    @inbounds c[k] = -x * s[k-1]
    @inbounds for i = 2:k
        x = i * a[i]
        s[k] += x * c[k-i]
        c[k] -= x * s[k-i]
    end
    @inbounds s[k] = s[k] / k
    @inbounds c[k] = c[k] / k
    return nothing
end

@inline function tan!(c::Taylor0, a::Taylor0, c2::Taylor0, k::Int)
    if k == 0
        @inbounds aux = tan(constant_term(a))
        @inbounds c[0] = aux
        @inbounds c2[0] = aux^2
        return nothing
    end
    @inbounds c[k] = k * a[k] * c2[0]
    @inbounds for i = 1:k-1
        c[k] += (k - i) * a[k-i] * c2[i]
    end
    @inbounds c[k] = a[k] + c[k] / k
    sqr!(c2, c, k)
    return nothing
end

@inline function asin!(c::Taylor0, a::Taylor0, r::Taylor0, k::Int)
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
    sqrt!(r, 1 - a^2, k)
    r0=constant_term(r)
    if r0==0.0
        r0=1e-9
    end
    @inbounds c[k] = (a[k] - c[k] / k) / r0
    return nothing
end

@inline function acos!(c::Taylor0, a::Taylor0, r::Taylor0, k::Int)
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
    sqrt!(r, 1 - a^2, k)
    r0=constant_term(r)
    if r0==0.0
        r0=1e-9
    end
    @inbounds c[k] = (a[k] + c[k] / k) / r0
    return nothing
end
   

@inline function atan!(c::Taylor0, a::Taylor0, r::Taylor0, k::Int)
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = atan(a0)
        @inbounds r[0] = 1 + a0^2
        return nothing
    end
    @inbounds c[k] = (k - 1) * r[1] * c[k-1]
    @inbounds for i in 2:k-1
        c[k] += (k - i) * r[i] * c[k-i]
    end
    @inbounds sqr!(r, a, k)
    r0=constant_term(r)
    if r0==0.0
        r0=1e-9
    end
    @inbounds c[k] = (a[k] + c[k] / k) / r0
    return nothing
end
 

@inline function sinhcosh!(s::Taylor0, c::Taylor0, a::Taylor0, k::Int)
    if k == 0
        @inbounds s[0] = sinh(constant_term(a))
        @inbounds c[0] = cosh(constant_term(a))
        return nothing
    end
    x = a[1]
    @inbounds s[k] = x * c[k-1]
    @inbounds c[k] = x * s[k-1]
    @inbounds for i = 2:k
        x = i * a[i]
        s[k] += x * c[k-i]
        c[k] += x * s[k-i]
    end
    s[k] = s[k] / k
    c[k] = c[k] / k
    return nothing
end

@inline function tanh!(c::Taylor0, a::Taylor0, c2::Taylor0, k::Int)
    if k == 0
        @inbounds aux = tanh(constant_term(a))
        @inbounds c[0] = aux
        @inbounds c2[0] = aux^2
        return nothing
    end
    @inbounds c[k] = k * a[k] * c2[0]
    @inbounds for i = 1:k-1
        c[k] += (k - i) * a[k-i] * c2[i]
    end
    @inbounds c[k] = a[k] - c[k] / k
    sqr!(c2, c, k)
    return nothing
end

@inline function asinh!(c::Taylor0, a::Taylor0, r::Taylor0, k::Int)
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = asinh(a0)
        @inbounds r[0] = sqrt(a0^2 + 1)
        return nothing
    end
    @inbounds c[k] = (k - 1) * r[1] * c[k-1]
    @inbounds for i in 2:k-1
        c[k] += (k - i) * r[i] * c[k-i]
    end
    sqrt!(r, a^2 + 1, k)   # a^...allocates will need a third cache
    r0=constant_term(r)
    if r0==0.0
        r0=1e-9
    end
    @inbounds c[k] = (a[k] - c[k] / k) / r0
    return nothing
end

@inline function acosh!(c::Taylor0, a::Taylor0, r::Taylor0, k::Int)
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = acosh(a0)
        @inbounds r[0] = sqrt(a0^2 - 1)
        return nothing
    end
    @inbounds c[k] = (k - 1) * r[1] * c[k-1]
    @inbounds for i in 2:k-1
        c[k] += (k - i) * r[i] * c[k-i]
    end
    sqrt!(r, a^2 - 1, k)
    r0=constant_term(r)
    if r0==0.0
        r0=1e-9
    end
    @inbounds c[k] = (a[k] - c[k] / k) / r0
    return nothing
end

@inline function atanh!(c::Taylor0, a::Taylor0, r::Taylor0, k::Int)
    if k == 0
        a0 = constant_term(a)
        @inbounds c[0] = atanh(a0)
        @inbounds r[0] = 1 - a0^2
        return nothing
    end
    @inbounds c[k] = (k - 1) * r[1] * c[k-1]
    @inbounds for i in 2:k-1
        c[k] += (k - i) * r[i] * c[k-i]
    end
    @inbounds sqr!(r, a, k)
    r0=constant_term(r)
    if r0==0.0
        r0=1e-9
    end
    @inbounds c[k] = (a[k] + c[k] / k) / r0
    return nothing
end
function abs(a::Taylor0)
    if constant_term(a) > 0
        return a
    elseif constant_term(a) < 0
        return -a
    else 
        throw(DomainError(a,
            """The 0th order Taylor0 coefficient must be non-zero
            (abs(x) is not differentiable at x=0)."""))
    end
end
function Derivative(y, x)
    #partial derivative of abs(x)/dx...needed for symbolic differentiation in the diff function inside extractdependency
    if x > 0
        return 1
    elseif x < 0
        return -1
    else
        return 0
    end
end







