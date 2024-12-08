# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
function ==(a::Taylor0, b::Taylor0)
    return a.coeffs == b.coeffs
end
zero(a::Taylor0) = Taylor0(zero.(a.coeffs))
function zero(a::Taylor0, cache::Taylor0)
    cache.coeffs .= 0.0
    return cache
end
function one(a::Taylor0, cache::Taylor0)
    cache.coeffs .= 0.0
    cache[0] = 1.0
    return cache
end
function one(a::Taylor0)
    b = zero(a)
    b[0] = one(b[0])
    return b
end
function (+)(a::Taylor0, b::Taylor0)
    v = similar(a.coeffs)
    @__dot__ v = (+)(a.coeffs, b.coeffs)
    return Taylor0(v, a.order)
end
function (+)(a::Taylor0)
    v = similar(a.coeffs)
    @__dot__ v = (+)(a.coeffs)
    return Taylor0(v, a.order)
end
function (+)(a::Taylor0, b::T) where {T<:Number}
    coeffs = copy(a.coeffs)
    @inbounds coeffs[1] = (+)(a[0], b)
    return Taylor0(coeffs, a.order)
end
(+)(b::T, a::Taylor0) where {T<:Number} = (+)(a, b)
function (-)(a::Taylor0, b::Taylor0)
    v = similar(a.coeffs)
    @__dot__ v = (-)(a.coeffs, b.coeffs)
    return Taylor0(v, a.order)
end
function (-)(a::Taylor0)
    v = similar(a.coeffs)
    @__dot__ v = (-)(a.coeffs)
    return Taylor0(v, a.order)
end
function (-)(a::Taylor0, b::T) where {T<:Number}
    coeffs = copy(a.coeffs)
    @inbounds coeffs[1] = (-)(a[0], b)
    return Taylor0(coeffs, a.order)
end
function (-)(b::T, a::Taylor0) where {T<:Number}
    coeffs = similar(a.coeffs)
    @__dot__ coeffs = (-)(a.coeffs)
    @inbounds coeffs[1] = (-)(b, a[0])
    return Taylor0(coeffs, a.order)
end
function *(a::Taylor0, b::Taylor0)
    c = Taylor0(zero(a[0]), a.order)
    for ord in eachindex(c)
        mul!(c, a, b, ord)
    end
    return c
end
function (*)(a::T, b::Taylor0) where {T<:Number}
    v = Array{T}(undef, length(b.coeffs))
    @__dot__ v = a * b.coeffs
    return Taylor0(v, b.order)
end
*(b::Taylor0, a::T) where {T<:Number} = a * b

@inline function mul!(c::Taylor0, a::Taylor0, b::Taylor0, k::Int)
    @inbounds c[k] = a[0] * b[k]
    @inbounds for i = 1:k
        c[k] += a[i] * b[k-i]
    end
    return nothing
end
function /(a::Taylor0, b::T) where {T<:Number}
    @inbounds aux = a.coeffs[1] / b
    v = Array{typeof(aux)}(undef, length(a.coeffs))
    @__dot__ v = a.coeffs / b
    return Taylor0(v, a.order)
end
function /(a::Taylor0, b::Taylor0)
    iszero(a) && !iszero(b) && return zero(a)
    ordfact, cdivfact = divfactorization(a, b)
    c = Taylor0(cdivfact, a.order - ordfact)
    for ord in eachindex(c)
        div!(c, a, b, ord)
    end
    return c
end
@inline function divfactorization(a1::Taylor0, b1::Taylor0)
    a1nz = findfirst(a1)
    b1nz = findfirst(b1)
    a1nz = a1nz ≥ 0 ? a1nz : a1.order
    b1nz = b1nz ≥ 0 ? b1nz : a1.order
    ordfact = min(a1nz, b1nz)
    cdivfact = a1[ordfact] / b1[ordfact]
    iszero(b1[ordfact]) && throw(ArgumentError(
        """Division does not define a Taylor0 polynomial;
        order k=$(ordfact) => coeff[$(ordfact)]=$(cdivfact)."""))
    return ordfact, cdivfact
end

@inline function div!(c::Taylor0, a::Taylor0, b::Taylor0, k::Int)
    ordfact, cdivfact = divfactorization(a, b)
    if k == 0
        @inbounds c[0] = cdivfact
        return nothing
    end
    imin = max(0, k + ordfact - b.order)
    @inbounds c[k] = c[imin] * b[k+ordfact-imin]
    @inbounds for i = imin+1:k-1
        c[k] += c[i] * b[k+ordfact-i]
    end
    if k + ordfact ≤ b.order
        @inbounds c[k] = (a[k+ordfact] - c[k]) / b[ordfact]
    else
        @inbounds c[k] = -c[k] / b[ordfact]
    end
    return nothing
end


