# This file is part of the Taylor0Series.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Conversion
convert(::Type{Taylor0}, a::Taylor0)  = Taylor0(convert(Array{Float64,1}, a.coeffs), a.order)

convert(::Type{Taylor0}, a::Taylor0)  = a

#convert(::Type{Taylor0{Rational{T}}}, a::Taylor0{S}) where {T<:Integer,S<:AbstractFloat} = Taylor0(rationalize.(a[:]), a.order)

convert(::Type{Taylor0}, b::Array{Float64,1})  =
    Taylor0(b, length(b)-1)

convert(::Type{Taylor0}, b::Array{S,1}) where {S<:Number} =
    Taylor0(convert(Array{Float64,1},b), length(b)-1)

convert(::Type{Taylor0}, b::S)  where {T<:Number,S<:Number} =
    Taylor0([convert(T,b)], 0)

convert(::Type{Taylor0}, b::T)   = Taylor0([b], 0)

convert(::Type{Taylor0}, a::T)  = Taylor0(a, 0)








# Promotion
#= promote_rule(::Type{Taylor0}, ::Type{Taylor0})  = Taylor0

promote_rule(::Type{Taylor0}, ::Type{Taylor0{S}}) where {T<:Number,S<:Number} =
    Taylor0{promote_type(T,S)} =#

promote_rule(::Type{Taylor0}, ::Type{Array{Float64,1}})  = Taylor0

#= promote_rule(::Type{Taylor0}, ::Type{Array{S,1}}) where {T<:Number,S<:Number} =
    Taylor0{promote_type(T,S)} =#


#promote_rule(::Type{Taylor0}, ::Type{T}) where {T<:Number} = Taylor0

#= promote_rule(::Type{Taylor0}, ::Type{S}) where {T<:Number,S<:Number} =
    Taylor0{promote_type(T,S)} =#

#= promote_rule(::Type{Taylor0{Taylor0}}, ::Type{Taylor0}) where {T<:Number} =
    Taylor0{Taylor0} =#


# Order may matter
#= promote_rule(::Type{S}, ::Type{T}) where {S<:NumberNotSeries,T<:AbstractSeries} =
    promote_rule(T,S)
# disambiguation with Base
promote_rule(::Type{Bool}, ::Type{T}) where {T<:AbstractSeries} = promote_rule(T, Bool)

promote_rule(::Type{S}, ::Type{T}) where
    {S<:AbstractIrrational,T<:AbstractSeries} = promote_rule(T,S)

 =#



# Nested Taylor0's
#= function promote(a::Taylor0{Taylor0}, b::Taylor0) where {T<:NumberNotSeriesN}
    order_a = get_order(a)
    order_b = get_order(b)
    zb = zero(b)
    new_bcoeffs = similar(a.coeffs)
    new_bcoeffs[1] = b
    @inbounds for ind in 2:order_a+1
        new_bcoeffs[ind] = zb
    end
    return a, Taylor0(b, order_a)
end
promote(b::Taylor0, a::Taylor0{Taylor0}) where {T<:NumberNotSeriesN} =
    reverse(promote(a, b))
 =#