# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Evaluating ##
#= """
    evaluate(a, [dx])

Evaluate a `Taylor0` polynomial using Horner's rule (hand coded). If `dx` is
ommitted, its value is considered as zero. Note that the syntax `a(dx)` is
equivalent to `evaluate(a,dx)`, and `a()` is equivalent to `evaluate(a)`.
""" =#
function evaluate(a::Taylor0, dx::T) where {T<:Number}
    @inbounds suma = a[end]
    @inbounds for k in a.order-1:-1:0
        suma = suma*dx + a[k]
    end
    suma
end
#= function evaluate(a::Taylor0, dx::S) where {S<:Number}
    suma = a[end]*one(dx)
    @inbounds for k in a.order-1:-1:0
        suma = suma*dx + a[k]
    end
    suma
end =#
evaluate(a::Taylor0)  = a[0]

#= """
    evaluate(x, δt)

Evaluates each element of `x::AbstractArray{Taylor0}`,
representing the dependent variables of an ODE, at *time* δt. Note that the
syntax `x(δt)` is equivalent to `evaluate(x, δt)`, and `x()`
is equivalent to `evaluate(x)`.
""" =#
#= evaluate(x::AbstractArray{Taylor0}, δt::S) where
    {T<:Number,S<:Number} = evaluate.(x, δt)
evaluate(a::AbstractArray{Taylor0}) where {T<:Number} = evaluate.(a, zero(T)) =#

#= """
    evaluate!(x, δt, x0)

Evaluates each element of `x::AbstractArray{Taylor0}`,
representing the Taylor expansion for the dependent variables
of an ODE at *time* `δt`. It updates the vector `x0` with the
computed values.
""" =#
#= function evaluate!(x::AbstractArray{Taylor0}, δt::S,
        x0::AbstractArray{T}) where {T<:Number,S<:Number}

    @inbounds for i in eachindex(x, x0)
        x0[i] = evaluate( x[i], δt )
    end
    nothing
end =#

#= """
    evaluate(a, x)

Substitute `x::Taylor0` as independent variable in a `a::Taylor0` polynomial.
Note that the syntax `a(x)` is equivalent to `evaluate(a, x)`.
""" =#
#= evaluate(a::Taylor0, x::Taylor0{S}) where {T<:Number,S<:Number} =
    evaluate(promote(a,x)...) =#

function evaluate(a::Taylor0, x::Taylor0) 
    if a.order != x.order
        a, x = fixorder(a, x)#allocates...delete later
    end
    @inbounds suma = a[end]*one(x)
    @inbounds for k = a.order-1:-1:0
        suma = suma*x + a[k]
    end
    suma
end
function evaluateX(a::Taylor0, x::Taylor0)  #not tested
    if a.order != x.order
        a, x = fixorder(a, x)
    end
    @inbounds suma = a[end]*one(x)
    @inbounds for k = a.order-1:-1:0
        suma = suma*x + a[k]
    end
    #println("suma",suma)
    @inbounds for k = 0:a.order-1
        a[k]=suma[k]
    end
    return nothing
end
#= function evaluate(a::Taylor0{Taylor0}, x::Taylor0) where {T<:Number}
    @inbounds suma = a[end]*one(x)
    @inbounds for k = a.order-1:-1:0
        suma = suma*x + a[k]
    end
    suma
end
function evaluate(a::Taylor0, x::Taylor0{Taylor0}) where {T<:Number}
    @inbounds suma = a[end]*one(x)
    @inbounds for k = a.order-1:-1:0
        suma = suma*x + a[k]
    end
    suma
end =#

evaluate(p::Taylor0, x::Array{S}) where {S<:Number} =
    evaluate.([p], x)

#function-like behavior for Taylor0
(p::Taylor0)(x) = evaluate(p, x)
(p::Taylor0)() = evaluate(p)

#function-like behavior for Vector{Taylor0}
#= if VERSION >= v"1.3"
    (p::AbstractArray{Taylor0})(x) where {T<:Number} = evaluate.(p, x)
    (p::AbstractArray{Taylor0})() where {T<:Number} = evaluate.(p)
else
    (p::Array{Taylor0})(x) where {T<:Number} = evaluate.(p, x)
    (p::SubArray{Taylor0})(x) where {T<:Number} = evaluate.(p, x)
    (p::Array{Taylor0})() where {T<:Number} = evaluate.(p)
    (p::SubArray{Taylor0})() where {T<:Number} = evaluate.(p)
end =#

