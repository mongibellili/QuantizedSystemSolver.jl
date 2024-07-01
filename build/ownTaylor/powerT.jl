# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#



powerT(a::Taylor0, n::Integer,cache1::Taylor0) = powerT(a, float(n),cache1)


## Real power ##
function powerT(a::Taylor0, r::S,cache1::Taylor0) where {S<:Real}
    if iszero(r) #later test r=0.0
        cache1[0]=1.0
         return cache1
    end
    
  #  @show(aa)
    r == 1 && return a
    r == 2 && return square(a,cache1)
    r == 1/2 && return sqrt(a,cache1)

    l0 = findfirst(a)
   # @show l0
    lnull = trunc(Int, r*l0 )
   # @show lnull
    if (a.order-lnull < 0) || (lnull > a.order)
       # @show a.order
        return cache1  #empty
    end
    c_order = l0 == 0 ? a.order : min(a.order, trunc(Int,r*a.order))
    # @show(c_order)
    #c = Taylor0(zero(aux), c_order)
    for k = 0:c_order
        pow!(cache1, a, r, k)
    end
    # println()
    return cache1
end
function powerT(a::T, r::S,cache1::Taylor0) where {S<:Real,T<:Number}
    cache1[0]=a^r
    return cache1
end



function square(a::Taylor0,cache1::Taylor0) #internal helper function ...no need to take care of a==number
        #c = Taylor0( constant_term(a)^2, a.order)
        cache1[0]=constant_term(a)^2
        for k in 1:a.order
            sqr!(cache1, a, k)
        end
        return cache1
 end




## Square root ##
function sqrt(a::Taylor0,cache1::Taylor0)
    l0nz = findfirst(a)
    if l0nz < 0
        cache1
    end
    cache1[0]=sqrt(a[0])
    for k = 1:a.order
        sqrt!(cache1, a, k, 0)
    end
    return cache1
end

function sqrt(a::T,cache1::Taylor0)where {T<:Number}
    cache1[0]=sqrt(a)
    return cache1
end