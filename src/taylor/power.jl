# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#
function ^(a::Taylor0, n::Integer)
    n == 0 && return one(a)
    n == 1 && return (a)
    n == 2 && return square(a)
    n < 0 && throw(DomainError("taylor^n & n<0 !!"))
    return power_by_squaring(a, n)
end

^(a::Taylor0, b::Taylor0) = exp(b * log(a))

function power_by_squaring(x::Taylor0, p::Integer)
    p == 1 && return (x)
    p == 0 && return one(x)
    p == 2 && return square(x)
    t = trailing_zeros(p) + 1
    p >>= t
    while (t -= 1) > 0
        x = square(x)
    end
    y = x
    while p > 0
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) ≥ 0
            x = square(x)
        end
        y *= x
    end
    return y
end
## Real power ##
function ^(a::Taylor0, r::S) where {S<:Real}
    a0 = constant_term(a)
    aux = one(a0)^r
    iszero(r) && return Taylor0(1.0, a.order)
    aa = aux * a
    r == 1 && return a
    r == 2 && return square(a)
    r == 1 / 2 && return sqrt(a)
    l0 = findfirst(a)
    lnull = trunc(Int, r * l0)
    if (a.order - lnull < 0) || (lnull > a.order)
        return Taylor0(0.0, a.order)
    end
    c_order = l0 == 0 ? a.order : min(a.order, trunc(Int, r * a.order))
    c = Taylor0(0.0, c_order)
    for k = 0:c_order
        pow!(c, aa, r, k)
    end
    return c
end

@inline function pow!(c::Taylor0, a::Taylor0, r::S, k::Int) where {S<:Real}
    l0 = findfirst(a)
    if l0 < 0
        c[k] = zero(a[0])
        return nothing
    end
    lnull = trunc(Int, r * l0)
    kprime = k - lnull
    if (kprime < 0) || (lnull > a.order)
        @inbounds c[k] = zero(a[0])
        return nothing
    end
    if isinteger(r) && r > 0 && (k > r * findlast(a))
        @inbounds c[k] = zero(a[0])
        return nothing
    end
    if k == lnull
        @inbounds c[k] = (a[l0])^r
        return nothing
    end
    # The recursion formula
    if l0 + kprime ≤ a.order
        @inbounds c[k] = r * kprime * c[lnull] * a[l0+kprime]
    else
        @inbounds c[k] = zero(a[0])
    end
    for i = 1:k-lnull-1
        ((i + lnull) > a.order || (l0 + kprime - i > a.order)) && continue
        aux = r * (kprime - i) - i
        @inbounds c[k] += aux * c[i+lnull] * a[l0+kprime-i]
    end
    @inbounds c[k] = c[k] / (kprime * a[l0])
    return nothing
end
function square(a::Taylor0)
    c = Taylor0(constant_term(a)^2, a.order)
    for k in 1:a.order
        sqr!(c, a, k)
    end
    return c
end

@inline function sqr!(c::Taylor0, a::Taylor0, k::Int)
    if k == 0
        @inbounds c[0] = constant_term(a)^2
        return nothing
    end
    kodd = k % 2
    kend = div(k - 2 + kodd, 2)
    @inbounds for i = 0:kend
        c[k] += a[i] * a[k-i]
    end
    @inbounds c[k] = 2 * c[k]
    kodd == 1 && return nothing
    @inbounds c[k] += a[div(k, 2)]^2
    return nothing
end
## Square root ##
function sqrt(a::Taylor0)
    l0nz = findfirst(a)
    aux = zero(a[0])
    if l0nz < 0
        return Taylor0(aux, a.order)
    end
    c = Taylor0(sqrt(a[0]), a.order)
    aa = one(aux) * a
    for k = 1:a.order
        sqrt!(c, aa, k, 0)
    end
    return c
end

@inline function sqrt!(c::Taylor0, a::Taylor0, k::Int, k0::Int=0)
    kodd = (k - k0) % 2
    kend = div(k - k0 - 2 + kodd, 2)
    imax = min(k0 + kend, a.order)
    imin = max(k0 + 1, k + k0 - a.order)
    if imin ≤ imax
        c[k] = c[imin] * c[k+k0-imin]
    end
    @inbounds for i = imin+1:imax
        c[k] += c[i] * c[k+k0-i]
    end
    if k + k0 ≤ a.order
        @inbounds aux = a[k+k0] - 2 * c[k]
    else
        @inbounds aux = -2 * c[k]
    end
    if kodd == 0
        @inbounds aux = aux - (c[kend+k0+1])^2
    end
    @inbounds c[k] = aux / (2 * c[k0])
    return nothing
end






