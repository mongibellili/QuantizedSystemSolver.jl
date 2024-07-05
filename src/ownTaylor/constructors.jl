# Float64his file is part of the Taylor0Series.jl Julia package, MIFloat64 license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIFloat64 Expat license
#




## Constructors ##


struct Taylor0 #<: AbstractSeries{Float64}
    coeffs :: Array{Float64,1}
    order :: Int
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

getcoeff(a::Taylor0, n::Int) = (@assert 0 ≤ n ≤ a.order; return a[n])

getindex(a::Taylor0, n::Int) = a.coeffs[n+1]
#getindex(a::Taylor0, u::UnitRange{Int}) = view(a.coeffs, u .+ 1 )
#getindex(a::Taylor0, c::Colon) = view(a.coeffs, c)
#getindex(a::Taylor0, u::StepRange{Int,Int})  =view(a.coeffs, u[:] .+ 1)

setindex!(a::Taylor0, x::T, n::Int) where {T<:Number} = a.coeffs[n+1] = x
#setindex!(a::Taylor0, x::T, u::UnitRange{Int}) where {T<:Number} = a.coeffs[u .+ 1] .= x
#= function setindex!(a::Taylor0, x::Array{T,1}, u::UnitRange{Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a.coeffs[u[ind]+1] = x[ind]
    end
end =#
#= setindex!(a::Taylor0, x::T, c::Colon) where {T<:Number} = a.coeffs[c] .= x
setindex!(a::Taylor0, x::Array{T,1}, c::Colon) where {T<:Number} = a.coeffs[c] .= x
setindex!(a::Taylor0, x::T, u::StepRange{Int,Int}) where {T<:Number} =
    a.coeffs[u[:] .+ 1] .= x
function setindex!(a::Taylor0, x::Array{T,1}, u::StepRange{Int,Int}) where {T<:Number}
    @assert length(u) == length(x)
    for ind in eachindex(x)
        a.coeffs[u[ind]+1] = x[ind]
    end
end =#



       
@inline iterate(a::Taylor0, state=0) = state > a.order ? nothing : (a.coeffs[state+1], state+1)
# Base.iterate(rS::Iterators.Reverse{Taylor0}, state=rS.itr.order) = state < 0 ? nothing : (a.coeffs[state], state-1)
@inline length(a::Taylor0) = length(a.coeffs)
@inline firstindex(a::Taylor0) = 0
@inline lastindex(a::Taylor0) = a.order

@inline eachindex(a::Taylor0) = firstindex(a):lastindex(a)
#@inline numtype(::Taylor0{S}) where {S<:Number} = S
@inline size(a::Taylor0) = size(a.coeffs)
@inline get_order(a::Taylor0) = a.order

#@inline axes(a::Taylor0) = ()

numtype(a) = eltype(a)



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



constant_term(a::Taylor0) = a[0]

constant_term(a::Vector{T}) where {T<:Number} = constant_term.(a)

constant_term(a::Number) = a

evaluate(p::Taylor0, x::Array{S}) where {S<:Number} =
    evaluate.([p], x)

#function-like behavior for Taylor0
function evaluate(a::Taylor0, dx::T) where {T<:Number}
    @inbounds suma = a[end]
    @inbounds for k in a.order-1:-1:0
        suma = suma*dx + a[k]
    end
    suma
end
evaluate(a::Taylor0)  = a[0]

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

(p::Taylor0)(x) = evaluate(p, x)
(p::Taylor0)() = evaluate(p)


function differentiate(a::Taylor0)
    res = Taylor0(zero(a[0]), a.order(a)-1)
    for ord in eachindex(res)
        differentiate!(res, a, ord)
    end
    return res
end

function differentiate!(p::Taylor0, a::Taylor0, k::Int)
    if k < a.order
        @inbounds p[k] = (k+1)*a[k+1]
    end
    return nothing
end
function differentiate!(::Val{O},cache::Taylor0, a::Taylor0)where {O}
    for k=0:O-1
       # differentiate!(res, a, ord)
       # if k < a.order
          @inbounds cache[k] = (k+1)*a[k+1]
      #end
    end
    cache[O]=0.0  #cache Letter O not zero
    nothing
end
function ndifferentiate!(cache::Taylor0,a::Taylor0, n::Int) 
    if n > a.order
        cache.coeffs.=0.0
    elseif n==0
        cache.coeffs.=a.coeffs
    else
        differentiate!(Val(a.order),cache,a)
        for i = 2:n
            differentiate!(Val(a.order),cache, cache)
        end
        #return Taylor0(view(res.coeffs, 1:a.order-n+1))
    end
end