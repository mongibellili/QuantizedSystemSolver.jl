# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#



    ## real, imag, conj and ctranspose ##
  
 #=    real(a::Taylor0) = Taylor0(real.(a.coeffs), a.order)
    imag(a::Taylor0) = Taylor0(imag.(a.coeffs), a.order)
    conj(a::Taylor0) = Taylor0(conj.(a.coeffs), a.order)
     adjoint(a::Taylor0) = conj(a)

    ## isinf and isnan ##
    isinf(a::Taylor0) = any( isinf.(a.coeffs) )

    isnan(a::Taylor0) = any( isnan.(a.coeffs) ) =#



## Division functions: rem and mod ##
#= for op in (:mod, :rem)
   
        @eval begin
            function ($op)(a::Taylor0, x::T) where {T<:Real}
                coeffs = copy(a.coeffs)
                @inbounds coeffs[1] = ($op)(constant_term(a), x)
                return Taylor0(coeffs, a.order)
            end

            function ($op)(a::Taylor0, x::S) where {T<:Real,S<:Real}
                R = promote_type(T, S)
                a = convert(Taylor0{R}, a)
                return ($op)(a, convert(R,x))
            end
        end
   

    
end =#


## mod2pi and abs ##

    
#= function mod2pi(a::Taylor0) where {T<:Real}
    coeffs = copy(a.coeffs)
    @inbounds coeffs[1] = mod2pi( constant_term(a) )
    return Taylor0( coeffs, a.order)
end =#

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
function abs(a::Taylor0,cache1::Taylor0)
    if constant_term(a) > 0
        return a
    elseif constant_term(a) < 0
        @__dot__ cache1.coeffs = (-)(a.coeffs)
        return cache1
    else
        cache1.coeffs .=Inf # no need to throw error, Inf is fine...for my solver i deal with it by guarding against small steps
        cache1[0]=0.0
        return cache1
        #= throw(DomainError(a, 
        """The 0th order Taylor0 coefficient must be non-zero
        (abs(x) is not differentiable at x=0).""")) =#
    end
end
function abs(a::T,cache1::Taylor0) where {T<:Number}
        cache1[0]=abs(a)
        return cache1
        #= throw(DomainError(a, 
        """The 0th order Taylor0 coefficient must be non-zero
        (abs(x) is not differentiable at x=0).""")) =#

end
#abs2(a::Taylor0) = real(a)^2 + imag(a)^2
#abs(x::Taylor0) where {T<:Complex} = sqrt(abs2(x))
#abs(x::Taylor0{Taylor0}) where {T<:Complex} = sqrt(abs2(x))

#= @doc doc"""
    abs(a)

For a `Real` type returns `a` if `constant_term(a) > 0` and `-a` if `constant_term(a) < 0` for
`a <:Union{Taylor0,TaylorN}`.
For a `Complex` type, such as `Taylor0{ComplexF64}`, returns `sqrt(real(a)^2 + imag(a)^2)`. 

Notice that `typeof(abs(a)) <: AbstractSeries` and 
that for a `Complex` argument a `Real` type is returned (e.g. `typeof(abs(a::Taylor0{ComplexF64})) == Taylor0{Float64}`).

""" abs


#norm
@doc doc"""
    norm(x::AbstractSeries, p::Real)

Returns the p-norm of an `x::AbstractSeries`, defined by

```math
\begin{equation*}
\left\Vert x \right\Vert_p =  \left( \sum_k | x_k |^p \right)^{\frac{1}{p}},
\end{equation*}
```
which returns a non-negative number.

""" norm =#

#= norm(x::AbstractSeries, p::Real=2) = norm( norm.(x.coeffs, p), p)
#norm for Taylor vectors
norm(v::Vector{T}, p::Real=2) where {T<:AbstractSeries} = norm( norm.(v, p), p)

# rtoldefault

    rtoldefault(::Type{Taylor0}) where {T<:Number} = rtoldefault(T)
    rtoldefault(::Taylor0) where {T<:Number} = rtoldefault(T) =#


#= # isfinite
"""
    isfinite(x::AbstractSeries) -> Bool

Test whether the coefficients of the polynomial `x` are finite.
"""
isfinite(x::AbstractSeries) = !isnan(x) && !isinf(x)

# isapprox; modified from Julia's Base.isapprox
"""
    isapprox(x::AbstractSeries, y::AbstractSeries; rtol::Real=sqrt(eps), atol::Real=0, nans::Bool=false)

Inexact equality comparison between polynomials: returns `true` if
`norm(x-y,1) <= atol + rtol*max(norm(x,1), norm(y,1))`, where `x` and `y` are
polynomials. For more details, see [`Base.isapprox`](@ref).
""" =#
#= function isapprox(x::T, y::S; rtol::Real=rtoldefault(x,y,0), atol::Real=0.0,
        nans::Bool=false) where {T<:AbstractSeries,S<:AbstractSeries}

    x == y || (isfinite(x) && isfinite(y) &&
        norm(x-y,1) <= atol + rtol*max(norm(x,1), norm(y,1))) ||
        (nans && isnan(x) && isnan(y))
end
#isapprox for vectors of Taylors
function isapprox(x::Vector{T}, y::Vector{S}; rtol::Real=rtoldefault(T,S,0), atol::Real=0.0,
        nans::Bool=false) where {T<:AbstractSeries,S<:AbstractSeries}

    x == y || norm(x-y,1) <= atol + rtol*max(norm(x,1), norm(y,1)) ||
        (nans && isnan(x) && isnan(y))
end

#taylor_expand function for Taylor0
#= """
    taylor_expand(f, x0; order)

Computes the Taylor expansion of the function `f` around the point `x0`.

If `x0` is a scalar, a `Taylor0` expansion will be returned. If `x0` is a vector,
a `TaylorN` expansion will be computed. If the dimension of x0 (`length(x0)`)
is different from the variables set for `TaylorN` (`get_numvars()`), an
`AssertionError` will be thrown.
""" =#
function taylor_expand(f::F; order::Int=15) where {F}
   a = Taylor0(order)
   return f(a)
end

function taylor_expand(f::F, x0::T; order::Int=15) where {F,T<:Number}
   a = Taylor0([x0, one(T)], order)
   return f(a)
end



#update! function for Taylor0
#= """
    update!(a, x0)

Takes `a <: Union{Taylo1,TaylorN}` and expands it around the coordinate `x0`.
""" =#
function update!(a::Taylor0, x0::T) where {T<:Number}
    a.coeffs .= evaluate(a, Taylor0([x0, one(x0)], a.order) ).coeffs
    nothing
end




    deg2rad(z::Taylor0) where {T<:AbstractFloat} = z * (convert(T, pi) / 180)
     deg2rad(z::Taylor0) where {T<:Real} = z * (convert(float(T), pi) / 180)

    rad2deg(z::Taylor0) where {T<:AbstractFloat} = z * (180 / convert(T, pi))
     rad2deg(z::Taylor0) where {T<:Real} = z * (180 / convert(float(T), pi))


# Internal mutating deg2rad!, rad2deg! functions

     @inline function deg2rad!(v::Taylor0, a::Taylor0, k::Int) where {T<:AbstractFloat}
        @inbounds v[k] = a[k] * (convert(T, pi) / 180)
        return nothing
    end
     @inline function deg2rad!(v::Taylor0{S}, a::Taylor0, k::Int) where {S<:AbstractFloat,T<:Real}
        @inbounds v[k] = a[k] * (convert(float(T), pi) / 180)
        return nothing
    end
    @inline function rad2deg!(v::Taylor0, a::Taylor0, k::Int) where {T<:AbstractFloat}
        @inbounds v[k] = a[k] * (180 / convert(T, pi))
        return nothing
    end
    @inline function rad2deg!(v::Taylor0{S}, a::Taylor0, k::Int) where {S<:AbstractFloat,T<:Real}
        @inbounds v[k] = a[k] * (180 / convert(float(T), pi))
        return nothing
    end

 =#