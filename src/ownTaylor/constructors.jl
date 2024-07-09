# Float64his file is part of the Taylor0Series.jl Julia package, MIFloat64 license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIFloat64 Expat license
#
## Constructors ##
"""Taylor0
A struct that defines a Taylor Variable. It has the following fields:\n
    - coeffs: An array of Float64 that holds the coefficients of the Taylor series\n
    - order: The order of the Taylor series
"""
struct Taylor0 #<: AbstractSeries{Float64}
    coeffs::Array{Float64,1}
    order::Int
end
Taylor0(x::Taylor0) = x
Taylor0(coeffs::Array{Float64,1}) = Taylor0(coeffs, length(coeffs) - 1)
function Taylor0(x::Float64, order::Int)
    v = fill(0.0, order + 1)
    v[1] = x
    return Taylor0(v, order)
end
getcoeff(a::Taylor0, n::Int) = (@assert 0 ≤ n ≤ a.order;return a[n])
getindex(a::Taylor0, n::Int) = a.coeffs[n+1]
setindex!(a::Taylor0, x::T, n::Int) where {T<:Number} = a.coeffs[n+1] = x
#@inline iterate(a::Taylor0, state=0) = state > a.order ? nothing : (a.coeffs[state+1], state + 1)
@inline length(a::Taylor0) = length(a.coeffs)
@inline firstindex(a::Taylor0) = 0
@inline lastindex(a::Taylor0) = a.order
@inline eachindex(a::Taylor0) = firstindex(a):lastindex(a)
@inline size(a::Taylor0) = size(a.coeffs)
@inline get_order(a::Taylor0) = a.order
#numtype(a) = eltype(a)
# Finds the first non zero entry; extended to Taylor0
function Base.findfirst(a::Taylor0)
    first = findfirst(x -> !iszero(x), a.coeffs)
    isa(first, Nothing) && return -1
    return first - 1
end
# Finds the last non-zero entry; extended to Taylor0
function Base.findlast(a::Taylor0)
    last = findlast(x -> !iszero(x), a.coeffs)
    isa(last, Nothing) && return -1
    return last - 1
end
constant_term(a::Taylor0) = a[0]
function evaluate(a::Taylor0, dx::T) where {T<:Number}
    @inbounds suma = a[end]
    @inbounds for k in a.order-1:-1:0
        suma = suma * dx + a[k]
    end
    suma
end
(p::Taylor0)(x) = evaluate(p, x)
#= function differentiate!(p::Taylor0, a::Taylor0, k::Int)
    if k < a.order
        @inbounds p[k] = (k + 1) * a[k+1]
    end
    return nothing
end =#
function differentiate!(::Val{O}, cache::Taylor0, a::Taylor0) where {O}
    for k = 0:O-1
        @inbounds cache[k] = (k + 1) * a[k+1]
    end
    cache[O] = 0.0  #cache Letter O not zero
    nothing
end
function ndifferentiate!(cache::Taylor0, a::Taylor0, n::Int)
    if n > a.order
        cache.coeffs .= 0.0
    elseif n == 0
        cache.coeffs .= a.coeffs
    else
        differentiate!(Val(a.order), cache, a)
        for i = 2:n
            differentiate!(Val(a.order), cache, cache)
        end
        #return Taylor0(view(res.coeffs, 1:a.order-n+1))
    end
end

function testTaylor(a::Taylor0)
    Taylor0(a) 
    getcoeff(a,1)
    length(a)
    size(a)
    get_order(a)

end