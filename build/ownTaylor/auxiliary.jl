# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Auxiliary function ##

#= """
    resize_coeffs1!{T<Number}(coeffs::Array{T,1}, order::Int)

If the length of `coeffs` is smaller than `order+1`, it resizes
`coeffs` appropriately filling it with zeros.
""" =#
#= function resize_coeffs1!(coeffs::Array{T,1}, order::Int) where {T<:Number}
    lencoef = length(coeffs)
    resize!(coeffs, order+1)
    #println("resize-coeffs")
    if order > lencoef-1
        z = zero(coeffs[1])
        @simd for ord in lencoef+1:order+1
            @inbounds coeffs[ord] = z
        end
    end
    return nothing
end =#

#= """
    resize_coeffsHP!{T<Number}(coeffs::Array{T,1}, order::Int)

If the length of `coeffs` is smaller than the number of coefficients
correspondinf to `order` (given by `size_table[order+1]`), it resizes
`coeffs` appropriately filling it with zeros.
""" =#
#= function resize_coeffsHP!(coeffs::Array{T,1}, order::Int) where {T<:Number}
    lencoef = length( coeffs )
    @inbounds num_coeffs = size_table[order+1]
    @assert order ≤ get_order() && lencoef ≤ num_coeffs
    num_coeffs == lencoef && return nothing
    resize!(coeffs, num_coeffs)
    z = zero(coeffs[1])
    @simd for ord in lencoef+1:num_coeffs
        @inbounds coeffs[ord] = z
    end
    return nothing
end =#

## Minimum order of an HomogeneousPolynomial compatible with the vector's length
#= function orderH(coeffs::Array{T,1}) where {T<:Number}
    ord = 0
    ll = length(coeffs)
    for i = 1:get_order()+1
        @inbounds num_coeffs = size_table[i]
        ll ≤ num_coeffs && break
        ord += 1
    end
    return ord
end =#





## getcoeff ##
#= """
    getcoeff(a, n)

Return the coefficient of order `n::Int` of a `a::Taylor0` polynomial.
""" =#
getcoeff(a::Taylor0, n::Int) = (@assert 0 ≤ n ≤ a.order; return a[n])

getindex(a::Taylor0, n::Int) = a.coeffs[n+1]
getindex(a::Taylor0, u::UnitRange{Int}) = view(a.coeffs, u .+ 1 )
getindex(a::Taylor0, c::Colon) = view(a.coeffs, c)
getindex(a::Taylor0, u::StepRange{Int,Int})  =
    view(a.coeffs, u[:] .+ 1)

setindex!(a::Taylor0, x::T, n::Int) where {T<:Number} = a.coeffs[n+1] = x
setindex!(a::Taylor0, x::T, u::UnitRange{Int}) where {T<:Number} =
    a.coeffs[u .+ 1] .= x
function setindex!(a::Taylor0, x::Array{T,1}, u::UnitRange{Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a.coeffs[u[ind]+1] = x[ind]
    end
end
setindex!(a::Taylor0, x::T, c::Colon) where {T<:Number} = a.coeffs[c] .= x
setindex!(a::Taylor0, x::Array{T,1}, c::Colon) where {T<:Number} = a.coeffs[c] .= x
setindex!(a::Taylor0, x::T, u::StepRange{Int,Int}) where {T<:Number} =
    a.coeffs[u[:] .+ 1] .= x
function setindex!(a::Taylor0, x::Array{T,1}, u::StepRange{Int,Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a.coeffs[u[ind]+1] = x[ind]
    end
end


#= function setorder(a::Taylor0, n::Int) #setfield!: immutable struct of type Taylor0 cannot be changed
    a.order=n
end =#


## eltype, length, get_order, etc ##

       
@inline iterate(a::Taylor0, state=0) = state > a.order ? nothing : (a.coeffs[state+1], state+1)
# Base.iterate(rS::Iterators.Reverse{Taylor0}, state=rS.itr.order) = state < 0 ? nothing : (a.coeffs[state], state-1)
@inline length(a::Taylor0) = length(a.coeffs)
@inline firstindex(a::Taylor0) = 0
@inline lastindex(a::Taylor0) = a.order

@inline eachindex(a::Taylor0) = firstindex(a):lastindex(a)
#@inline numtype(::Taylor0{S}) where {S<:Number} = S
@inline size(a::Taylor0) = size(a.coeffs)
@inline get_order(a::Taylor0) = a.order

@inline axes(a::Taylor0) = ()

numtype(a) = eltype(a)

#= @doc doc"""
    numtype(a::AbstractSeries)

Returns the type of the elements of the coefficients of `a`.
""" numtype =#

# Dumb methods included to properly export normalize_taylor (if IntervalArithmetic is loaded)
#@inline normalize_taylor(a::AbstractSeries) = a


## fixorder ##

#= @inline function fixorder(a::Taylor0, b::Taylor0)
    a.order == b.order && return a, b
    minorder, maxorder = minmax(a.order, b.order)
    if minorder ≤ 0
        minorder = maxorder
    end
    return Taylor0(copy(a.coeffs), minorder), Taylor0(copy(b.coeffs), minorder)
end =#

# can i go to the higher order when fixing orders??
@inline function fixorder(a::Taylor0, b::Taylor0)
   # a.order == b.order && return a, b
    if a.order < b.order
        resize_coeffs1!(a.coeffs,b.order)
       return Taylor0(a.coeffs, b.order),b
    elseif a.order > b.order
        resize_coeffs1!(b.coeffs,a.order)
        return a, Taylor0(b.coeffs, a.order)
    end
    
end

# Finds the first non zero entry; extended to Taylor0
function Base.findfirst(a::Taylor0) 
    first = findfirst(x->!iszero(x), a.coeffs)
    isa(first, Nothing) && return -1
    return first-1
end
# Finds the last non-zero entry; extended to Taylor0
function Base.findlast(a::Taylor0) 
    last = findlast(x->!iszero(x), a.coeffs)
    isa(last, Nothing) && return -1
    return last-1
end


## copyto! ##
# Inspired from base/abstractarray.jl, line 665
 function copyto!(dst::Taylor0, src::Taylor0) 
        length(dst) < length(src) && throw(ArgumentError(string("Destination has fewer elements than required; no copy performed")))
        destiter = eachindex(dst)
        y = iterate(destiter)
        for x in src
            dst[y[1]] = x
            y = iterate(destiter, y[2])
        end
        return dst
end


#= 
"""
    constant_term(a)

Return the constant value (zero order coefficient) for `Taylor0`
and `TaylorN`. The fallback behavior is to return `a` itself if
`a::Number`, or `a[1]` when `a::Vector`.
""" =#
constant_term(a::Taylor0) = a[0]

constant_term(a::Vector{T}) where {T<:Number} = constant_term.(a)

constant_term(a::Number) = a

#= """
    linear_polynomial(a)

Returns the linear part of `a` as a polynomial (`Taylor0` or `TaylorN`),
*without* the constant term. The fallback behavior is to return `a` itself.
""" =#
#= linear_polynomial(a::Taylor0) = Taylor0([zero(a[1]), a[1]], a.order)

linear_polynomial(a::Vector{T}) where {T<:Number} = linear_polynomial.(a)

linear_polynomial(a::Number) = a

#= """
    nonlinear_polynomial(a)

Returns the nonlinear part of `a`. The fallback behavior is to return `zero(a)`.
""" =#
nonlinear_polynomial(a::AbstractSeries) = a - constant_term(a) - linear_polynomial(a)

nonlinear_polynomial(a::Vector{T}) where {T<:Number} = nonlinear_polynomial.(a)

nonlinear_polynomial(a::Number) = zero(a) =#
