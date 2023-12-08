# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#



#= The following method computes `a^float(n)` (except for cases like
Taylor0{Interval{T}}^n, where `power_by_squaring` is used), to
use internally `pow!`.
=#
#^(a::Taylor0, n::Integer) = a^float(n)





    function ^(a::Taylor0, n::Integer) 
        n == 0 && return one(a)
        n == 1 && return copy(a)
        n == 2 && return square(a)
        n < 0 && throw(DomainError("taylor^n & n<0 !!"))
        return power_by_squaring(a, n)
    end

 #=   function ^(a::Taylor0{Rational{T}}, n::Integer) where {T<:Integer}
        n == 0 && return one(a)
        n == 1 && return copy(a)
        n == 2 && return square(a)
        n < 0 && return inv( a^(-n) )
        return power_by_squaring(a, n)
    end =#

    ^(a::Taylor0, x::Rational) = a^(x.num/x.den)

     ^(a::Taylor0, b::Taylor0) = exp( b*log(a) )

     ^(a::Taylor0, x::T) where {T<:Complex} = exp( x*log(a) )



# power_by_squaring; slightly modified from base/intfuncs.jl
# Licensed under MIT "Expat"

     function power_by_squaring(x::Taylor0, p::Integer)
        p == 1 && return copy(x)
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
    # println()
    #a0 = constant_term(a)
    # @show(a, a0)
    #aux = one(a0)^r
    # @show(aux)

    iszero(r) && return Taylor0(1.0, a.order)
   # aa = aux*a
  #  @show(aa)
    r == 1 && return a
    r == 2 && return square(a)
    r == 1/2 && return sqrt(a)

    l0 = findfirst(a)
   # @show l0
    lnull = trunc(Int, r*l0 )
   # @show lnull
    if (a.order-lnull < 0) || (lnull > a.order)
       # @show a.order
        return Taylor0( 0.0, a.order)
    end
    c_order = l0 == 0 ? a.order : min(a.order, trunc(Int,r*a.order))
    # @show(c_order)
    c = Taylor0(0.0, c_order)
    for k = 0:c_order
        pow!(c, aa, r, k)
    end
    # println()
    return c
end



# Homogeneous coefficients for real power
#= @doc doc"""
    pow!(c, a, r::Real, k::Int)

Update the `k`-th expansion coefficient `c[k]` of `c = a^r`, for
both `c` and `a` either `Taylor0` or `TaylorN`.

The coefficients are given by

```math
c_k = \frac{1}{k a_0} \sum_{j=0}^{k-1} \big(r(k-j) -j\big)a_{k-j} c_j.
```

For `Taylor0` polynomials, a similar formula is implemented which
exploits `k_0`, the order of the first non-zero coefficient of `a`.

"""  =#

@inline function pow!(c::Taylor0, a::Taylor0, r::S, k::Int) where {S<:Real}

    if r == 0
        return one!(c, a, k)
    elseif r == 1
        return identity!(c, a, k)
    elseif r == 2
        return sqr!(c, a, k)
    elseif r == 0.5
        return sqrt!(c, a, k)
    end

    # First non-zero coefficient
    l0 = findfirst(a)
    if l0 < 0
        c[k] = zero(a[0])
        return nothing
    end

    # The first non-zero coefficient of the result; must be integer
   #=  !isinteger(r*l0) && throw(DomainError(a,
        """The 0th order Taylor0 coefficient must be non-zero
        to raise the Taylor0 polynomial to a non-integer exponent.""")) =#
    lnull = trunc(Int, r*l0 )
    kprime = k-lnull
    if (kprime < 0) || (lnull > a.order)
        @inbounds c[k] = zero(a[0])
        return nothing
    end

    # Relevant for positive integer r, to avoid round-off errors
    if isinteger(r) && r > 0 && (k > r*findlast(a))
        @inbounds c[k] = zero(a[0])
        return nothing
    end

    if k == lnull
        @inbounds c[k] = (a[l0])^r
        return nothing
    end

    # The recursion formula
    if l0+kprime ≤ a.order
        @inbounds c[k] = r * kprime * c[lnull] * a[l0+kprime]
    else
        @inbounds c[k] = zero(a[0])
    end
    for i = 1:k-lnull-1
        ((i+lnull) > a.order || (l0+kprime-i > a.order)) && continue
        aux = r*(kprime-i) - i
        @inbounds c[k] += aux * c[i+lnull] * a[l0+kprime-i]
    end
    @inbounds c[k] = c[k] / (kprime * a[l0])
    return nothing
end




#= 
## Square ##
"""
    square(a::AbstractSeries) --> typeof(a)

Return `a^2`; see [`TaylorSeries.sqr!`](@ref).
""" 
 =#

   function square(a::Taylor0)
        c = Taylor0( constant_term(a)^2, a.order)
        for k in 1:a.order
            sqr!(c, a, k)
        end
        return c
    end




#= # Homogeneous coefficients for square
@doc doc"""
    sqr!(c, a, k::Int) --> nothing

Update the `k-th` expansion coefficient `c[k]` of `c = a^2`, for
both `c` and `a` either `Taylor0` or `TaylorN`.

The coefficients are given by

```math
\begin{eqnarray*}
c_k & = & 2 \sum_{j=0}^{(k-1)/2} a_{k-j} a_j,
    \text{ if k is odd,} \\
c_k & = & 2 \sum_{j=0}^{(k-2)/2} a_{k-j} a_j + (a_{k/2})^2,
    \text{ if k is even. }
\end{eqnarray*}
```

""" sqr! =#


   
        @inline function sqr!(c::Taylor0, a::Taylor0, k::Int) 
            if k == 0
                @inbounds c[0] = constant_term(a)^2
                return nothing
            end

            kodd = k%2
            kend = div(k - 2 + kodd, 2)
            @inbounds for i = 0:kend
               # if Taylor0 == Taylor0
                    c[k] += a[i] * a[k-i]
               #=  else
                    mul!(c[k], a[i], a[k-i])
                end =#
            end
            @inbounds c[k] = 2 * c[k]
            kodd == 1 && return nothing

            #if Taylor0 == Taylor0
                @inbounds c[k] += a[div(k,2)]^2
            #= else
                sqr!(c[k], a[div(k,2)])
            end =#

            return nothing
        end
  






## Square root ##
function sqrt(a::Taylor0)

    # First non-zero coefficient
    l0nz = findfirst(a)
   # println("l0nz= ",l0nz ) 
    #aux = zero(sqrt( constant_term(a) ))
    aux=zero(a[0])
   # @show aux
#=     println(" constant_term(a)= ",constant_term(a))
    println("sqrt( constant_term(a)= ",sqrt(constant_term(a)))
    println("aux= ",aux) =#
    if l0nz < 0
     #   println("l0nz < 0 ",Taylor0(aux, a.order))
        return Taylor0(aux, a.order)
    #elseif l0nz%2 == 1 # l0nz must be pair
        #= throw(DomainError(a,
        """First non-vanishing Taylor0 coefficient must correspond
        to an **even power** in order to expand `sqrt` around 0.""")) =#
    end

    # The last l0nz coefficients are set to zero.
   # lnull = l0nz >> 1 # integer division by 2
  
    #println("lnull= ",lnull)
   
   # c_order = l0nz == 0 ? a.order : a.order >> 1
  
   # println("a_order= ",a.order)
   # println("c_order= ",c_order)
    c = Taylor0( sqrt( a[0] ), a.order )
   # println("c= ",c)
   
   # println("c= ",c)
    aa = one(aux) * a
    #println("aa= ",aa) #a int aa 
    #println("aa=a ",aa===a)
    for k = 1:a.order
        sqrt!(c, aa, k, 0)
       # println("after sqrt!c= ",c)
    end

    return c
end



# Homogeneous coefficients for the square-root
#= @doc doc"""
    sqrt!(c, a, k::Int, k0::Int=0)

Compute the `k-th` expansion coefficient `c[k]` of `c = sqrt(a)`
for both`c` and `a` either `Taylor0` or `TaylorN`.

The coefficients are given by

```math
\begin{eqnarray*}
c_k &=& \frac{1}{2 c_0} \big( a_k - 2 \sum_{j=1}^{(k-1)/2} c_{k-j}c_j\big),
    \text{ if k is odd,} \\
c_k &=& \frac{1}{2 c_0} \big( a_k - 2 \sum_{j=1}^{(k-2)/2} c_{k-j}c_j
    - (c_{k/2})^2\big), \text{ if k is even.}
\end{eqnarray*}
```

For `Taylor0` polynomials, `k0` is the order of the first non-zero
coefficient, which must be even.

""" sqrt! =#

@inline function sqrt!(c::Taylor0, a::Taylor0, k::Int, k0::Int=0) 
    
    #println("k= ",k)
   # println("k0= ",k0)
    if k == k0
        @inbounds c[k] = sqrt(a[2*k0])
       println("i think this is a dead branch")
        return nothing
    end

    kodd = (k - k0)%2
   # println("kodd= ",kodd)
    kend = div(k - k0 - 2 + kodd, 2)
   # println("kend= ",kend)
    imax = min(k0+kend, a.order)
   # println("imax= ",imax)
    imin = max(k0+1, k+k0-a.order)
   # println("imin= ",imin)
    #imin ≤ imax && ( @inbounds  println("c in loop= ",c) c[k] = c[imin] * c[k+k0-imin] )
    if imin ≤ imax 
       # println("c in imin imax= ",c) 
         c[k] = c[imin] * c[k+k0-imin] 
    end
   # println("c= ",c)
    @inbounds for i = imin+1:imax
        c[k] += c[i] * c[k+k0-i]
     #   println("c in loop= ",c)
    end
    if k+k0 ≤ a.order
        @inbounds aux = a[k+k0] - 2*c[k]
      #  println("aux if <order= ",aux)
    else
        @inbounds aux = - 2*c[k]
      #  println("aux else= ",aux)
    end
    if kodd == 0
        @inbounds aux = aux - (c[kend+k0+1])^2
      #  println("aux if kodd==0  = ",aux)
    end
    @inbounds c[k] = aux / (2*c[k0])
    #println("final c inside sqrt!= ",c)
    return nothing
end


