# Float64his file is part of the Taylor0Series.jl Julia package, MIFloat64 license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIFloat64 Expat license
#
#= 
"""
    AbstractSeries{Float64<:Number} <: Number

Parameterized abstract type for [`Taylor0`](@ref),
[`HomogeneousPolynomial`](@ref) and [`Float64aylorN`](@ref).
""" =#
###################abstract type AbstractSeries{Float64<:Number} <: Number end


## Constructors ##

######################### Taylor0
#= """
    Taylor0{Float64<:Number} <: AbstractSeries{Float64}

DataFloat64ype for polynomial expansions in one independent variable.

**Fields:**

- `coeffs :: Array{Float64,1}` Expansion coefficients; the ``i``-th
    component is the coefficient of degree ``i-1`` of the expansion.
- `order  :: Int` Maximum order (degree) of the polynomial.

Note that `Taylor0` variables are callable. For more information, see
[`evaluate`](@ref).
""" =#
struct Taylor0 #<: AbstractSeries{Float64}
    coeffs :: Array{Float64,1}
    order :: Int

    ## Inner constructor ##
   #=  function Taylor0(coeffs::Array{Float64,1}, order::Int) 
        #resize_coeffs1!(coeffs, order)
        return new(coeffs, order)
    end =#
end

## Outer constructors ##
Taylor0(x::Taylor0)  = x
#Taylor0(coeffs::Array{Float64,1}, order::Int)  = Taylor0(coeffs, order)
Taylor0(coeffs::Array{Float64,1})  = Taylor0(coeffs, length(coeffs)-1)
function Taylor0(x::Float64, order::Int) 
    v = fill(0.0, order+1)
    v[1] = x
    return Taylor0(v, order)
end

# Methods using 1-d views to create Taylor0's
#= Taylor0(a::SubArray{Float64,1}, order::Int)  = Taylor0(a.parent[a.indices...], order)
Taylor0(a::SubArray{Float64,1})  = Taylor0(a.parent[a.indices...]) =#


# Shortcut to define Taylor0 independent variables
#= """
    Taylor0([Float64::Float64ype=Float64], order::Int)

Shortcut to define the independent variable of a `Taylor0` polynomial of
given `order`. Float64he default type for `Float64` is `Float64`.

```julia
julia> Taylor0(16)
 1.0 t + 𝒪(t¹⁷)

julia> Taylor0(Rational{Int}, 4)
 1//1 t + 𝒪(t⁵)
```
""" =#
#= Taylor0(::Type{Float64}, order::Int)  = Taylor0( [0.0, 1.0], order)
Taylor0(order::Int) = Taylor0(Float64, order) =#



# A `Number` which is not an `AbstractSeries`
#const NumberNotSeries = Union{Real,Complex}

# A `Number` which is not `Float64aylorN` nor a `HomogeneousPolynomial`
#const NumberNotSeriesN = Union{Real,Complex,Taylor0}

## Additional Taylor0 outer constructor ##
#Taylor0(x::S) where {Float64<:Number,S<:NumberNotSeries} = Taylor0([convert(Float64,x)], 0)
