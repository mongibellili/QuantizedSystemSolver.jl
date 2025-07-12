(addsub)(a, b, c) = (-)((+)(a, b), c)
(subsub)(a, b, c) = (-)((-)(a, b), c)
(subadd)(a, b, c) = (+)((-)(a, b), c)# 
function addTfoldl(op, res, bs...)# op is only addt so remove later
  l = length(bs)
  l == 1 && return res
  i = 2
  l == i && return op(res, bs[1], bs[end])
  i = 3
  l == i && return op(res, bs[1], bs[2], bs[end])
  i = 4
  l == i && return op(op(res, bs[1], bs[2], bs[end]), bs[3], bs[end])
  i = 5
  l == i && return op(op(res, bs[1], bs[2], bs[end]), bs[3], bs[4], bs[end])
  i = 6
  l == i && return op(op(op(res, bs[1], bs[2], bs[end]), bs[3], bs[4], bs[end]), bs[5], bs[end])
  i = 7
  l == i && return op(op(op(res, bs[1], bs[2], bs[end]), bs[3], bs[4], bs[end]), bs[5], bs[6], bs[end])
  i = 8
  l == i && return op(op(op(op(res, bs[1], bs[2], bs[end]), bs[3], bs[4], bs[end]), bs[5], bs[6], bs[end]), bs[7], bs[end])
  i = 9
  l == i && return op(op(op(op(res, bs[1], bs[2], bs[end]), bs[3], bs[4], bs[end]), bs[5], bs[6], bs[end]), bs[7], bs[8], bs[end])
  #last i creates y and passes it to the next loop: note that from this point, op will allocate
  #  warn users . (inserting -0 will avoid allocation!)
  for i in i:l-1 # -1 cuz last one is cache...but these lines are never reached cuz macro does not create mulT for args >7 in the first place
    y = op(y, bs[i], bs[end])
  end
  return y
end
function addT(a, b, c, d, xs...)
  if length(xs) != 0# if ==0 case below where d==cache     
    addTfoldl(addT, addT(addT(a, b, c, xs[end]), d, xs[end]), xs...)
  end
end
function mulTfoldl(res::P, bs...) where {P<:Taylor0}
  l = length(bs)
  i = 2
  l == i && return res
  i = 3
  l == i && return mulTT(res, bs[1], bs[end-1], bs[end])
  i = 4
  l == i && return mulTT(mulTT(res, bs[1], bs[end-1], bs[end]), bs[2], bs[end-1], bs[end])
  i = 5
  l == i && return mulTT(mulTT(mulTT(res, bs[1], bs[end-1], bs[end]), bs[2], bs[end-1], bs[end]), bs[3], bs[end-1], bs[end])
  i = 6
  l == i && return mulTT(mulTT(mulTT(mulTT(res, bs[1], bs[end-1], bs[end]), bs[2], bs[end-1], bs[end]), bs[3], bs[end-1], bs[end]), bs[4], bs[end-1], bs[end])
end
function mulTT(a::P, b::R, c::Q, xs...) where {P,Q,R<:Union{Taylor0,Number}}
  if length(xs) > 1
    mulTfoldl(mulTT(mulTT(a, b, xs[end-1], xs[end]), c, xs[end-1], xs[end]), xs...)
  end
end

"""
    createT(a::Taylor0, cache::Taylor0)

puts the coefficients of `a` into cache and returns cache.


```jldoctest
using QuantizedSystemSolver
a = Taylor0([1.0, 1.0])
cache = Taylor0([0.0, 0.0])
createT(a, cache)
cache[0]

# output

1.0
```
"""
function createT(a::Taylor0, cache::Taylor0)
  @__dot__ cache.coeffs = a.coeffs
  return cache
end
"""
    createT(a::T,cache::Taylor0) where {T<:Number}
   
    creates a Taylor0 from a constant. In case of order 2, cache=[a,0,0]
"""
function createT(a::T, cache::Taylor0) where {T<:Number} # requires cache1 clean
  cache[0] = a
  return cache
end
"""
    addsub(a::Taylor0, b::Taylor0,c::Taylor0,cache::Taylor0) 
   cache=a+b-c
"""
function addsub(a::Taylor0, b::Taylor0, c::Taylor0, cache::Taylor0)
  @__dot__ cache.coeffs = (-)((+)(a.coeffs, b.coeffs), c.coeffs)
  return cache
end
"""
    addsub(a::Taylor0, b::Taylor0, c::T, cache::Taylor0) where {T<:Number}

Performs an addition and subtraction operation on Taylor series objects.

# Arguments
- `a::Taylor0`: The first Taylor series object.
- `b::Taylor0`: The second Taylor series object.
- `c::T`: A scalar value to be subtracted from the constant term of `b`.
- `cache::Taylor0`: A Taylor series object used as a cache to store the result.

# Returns
- `cache::Taylor0`: The result of the operation, stored in the `cache` object.

# Description
This function performs the following operations:
1. Copies the coefficients of `b` into `cache`.
2. Subtracts the scalar `c` from the constant term of `cache`.
3. Adds the coefficients of `a` to `cache`.

The result is stored in the `cache` object and returned.
"""
function addsub(a::Taylor0, b::Taylor0, c::T, cache::Taylor0) where {T<:Number}
  cache.coeffs .= b.coeffs
  cache[0] = cache[0] - c
  @__dot__ cache.coeffs = (+)(cache.coeffs, a.coeffs)
  return cache
end
"""
    addsub(a::T, b::Taylor0,c::Taylor0,cache::Taylor0) where {T<:Number} 
#  Example:
  Order2 case: cache=[a+b[0]-c[0],b[1]-c[1],b[2]-c[2]]
"""
function addsub(a::T, b::Taylor0, c::Taylor0, cache::Taylor0) where {T<:Number}
  cache.coeffs .= b.coeffs
  cache[0] = a + cache[0]
  @__dot__ cache.coeffs = (-)(cache.coeffs, c.coeffs)
  return cache
end
addsub(a::Taylor0, b::T, c::Taylor0, cache::Taylor0) where {T<:Number} = addsub(b, a, c, cache)#addsub(b::T, a::Taylor0 ,c::Taylor0,cache::Taylor0) where {T<:Number}
function addsub(a::T, b::T, c::Taylor0, cache::Taylor0) where {T<:Number}
  @__dot__ cache.coeffs = (-)(c.coeffs)
  cache[0] = a + b + cache[0]
  return cache
end
function addsub(c::Taylor0, a::T, b::T, cache::Taylor0) where {T<:Number}
  cache.coeffs .= c.coeffs
  cache[0] = a + cache[0] - b
  return cache
end
addsub(a::T, c::Taylor0, b::T, cache::Taylor0) where {T<:Number} = addsub(c, a, b, cache)
function addsub(c::T, a::T, b::T, cache::Taylor0) where {T<:Number}#requires cache clean
  cache[0] = a + c - b
  return cache
end
###########"negate###########
"""
    negateT(a::Taylor0,cache::Taylor0)
cache=-a
"""
function negateT(a::Taylor0, cache::Taylor0) #where {T<:Number}
  @__dot__ cache.coeffs = (-)(a.coeffs)
  return cache
end
function negateT(a::T, cache::Taylor0) where {T<:Number} # requires cache1 clean
  cache[0] = -a
  return cache
end
#################################################subsub########################################################"
"""
    subsub(a::Taylor0, b::Taylor0,c::Taylor0,cache::Taylor0) 
cache=a-b-c
"""
function subsub(a::Taylor0, b::Taylor0, c::Taylor0, cache::Taylor0)
  @__dot__ cache.coeffs = subsub(a.coeffs, b.coeffs, c.coeffs)
  return cache
end
function subsub(a::Taylor0, b::Taylor0, c::T, cache::Taylor0) where {T<:Number}
  cache.coeffs .= b.coeffs
  cache[0] = cache[0] + c     #(a-b-c=a-(b+c))
  @__dot__ cache.coeffs = (-)(a.coeffs, cache.coeffs)
  return cache
end
subsub(a::Taylor0, b::T, c::Taylor0, cache::Taylor0) where {T<:Number} = subsub(a, c, b, cache)
function subsub(a::T, b::Taylor0, c::Taylor0, cache::Taylor0) where {T<:Number}
  @__dot__ cache.coeffs = (-)(b.coeffs)
  cache[0] = a + cache[0]
  @__dot__ cache.coeffs = (-)(cache.coeffs, c.coeffs)
  return cache
end
function subsub(a::Taylor0, b::T, c::T, cache::Taylor0) where {T<:Number}
  cache.coeffs .= a.coeffs
  cache[0] = cache[0] - b - c
  return cache
end
function subsub(a::T, b::Taylor0, c::T, cache::Taylor0) where {T<:Number}
  @__dot__ cache.coeffs = (-)(b.coeffs)
  cache[0] = cache[0] + a - c
  return cache
end
subsub(b::T, c::T, a::Taylor0, cache::Taylor0) where {T<:Number} = subsub(b, a, c, cache)
function subsub(a::T, b::T, c::T, cache::Taylor0) where {T<:Number}#require emptying the cache
  cache[0] = a - b - c
  return cache
end
#################################subadd####################################""
function subadd(a::P, b::Q, c::R, cache::Taylor0) where {P,Q,R<:Union{Taylor0,Number}}
  addsub(a, c, b, cache)
end
##################################""subT################################
"""
    subT(a::Taylor0, b::Taylor0,cache::Taylor0)
cache=a-b
"""
function subT(a::Taylor0, b::Taylor0, cache::Taylor0)# where {T<:Number}
  @__dot__ cache.coeffs = (-)(a.coeffs, b.coeffs)
  return cache
end
function subT(a::Taylor0, b::T, cache::Taylor0) where {T<:Number}
  @__dot__ cache.coeffs = a.coeffs
  cache[0] = cache[0] - b
  return cache
end
function subT(a::T, b::Taylor0, cache::Taylor0) where {T<:Number}
  @__dot__ cache.coeffs = (-)(b.coeffs)
  cache[0] = cache[0] + a
  return cache
end
function subT(a::T, b::T, cache::Taylor0) where {T<:Number} ##require emptying the cache because cache[1] and higher contains stuff
  cache[0] = a - b
  return cache
end
"""
    addT(a::Taylor0, b::Taylor0,cache::Taylor0) 
cache=a+b
"""
function addT(a::Taylor0, b::Taylor0, cache::Taylor0)
  @__dot__ cache.coeffs = (+)(a.coeffs, b.coeffs)
  return cache
end
function addT(a::Taylor0, b::T, cache::Taylor0) where {T<:Number}
  cache.coeffs .= a.coeffs
  cache[0] = cache[0] + b
  return cache
end
addT(b::T, a::Taylor0, cache::Taylor0) where {T<:Number} = addT(a, b, cache)
function addT(a::T, b::T, cache::Taylor0) where {T<:Number}#require emptying the cache
  #cache.coeffs.=0.0 # 
  cache[0] = a + b
  return cache
end
#add Three vars a b c
function addT(a::Taylor0, b::Taylor0, c::Taylor0, cache::Taylor0)
  @__dot__ cache.coeffs = (+)(a.coeffs, b.coeffs, c.coeffs)
  return cache
end
function addT(a::Taylor0, b::Taylor0, c::T, cache::Taylor0) where {T<:Number}
  @__dot__ cache.coeffs = (+)(a.coeffs, b.coeffs)
  cache[0] = cache[0] + c
  return cache
end
addT(a::Taylor0, c::T, b::Taylor0, cache::Taylor0) where {T<:Number} = addT(a, b, c, cache)
addT(c::T, a::Taylor0, b::Taylor0, cache::Taylor0) where {T<:Number} = addT(a, b, c, cache)
function addT(a::Taylor0, b::T, c::T, cache::Taylor0) where {T<:Number}
  cache.coeffs .= a.coeffs
  cache[0] = cache[0] + c + b
  return cache
end
addT(b::T, a::Taylor0, c::T, cache::Taylor0) where {T<:Number} = addT(a, b, c, cache)
addT(b::T, c::T, a::Taylor0, cache::Taylor0) where {T<:Number} = addT(a, b, c, cache)
function addT(a::T, b::T, c::T, cache::Taylor0) where {T<:Number}#require clean cache
  cache[0] = a + b + c
  return cache
end
#mul uses one cache when the original mul has only two elemnts a * b
function mulT(a::Taylor0, b::T, cache1::Taylor0) where {T<:Number}
  fill!(cache1.coeffs, b)
  @__dot__ cache1.coeffs = a.coeffs * cache1.coeffs  
  return cache1
end
mulT(a::T, b::Taylor0, cache1::Taylor0) where {T<:Number} = mulT(b, a, cache1)
"""
    mulT(a::Taylor0, b::Taylor0,cache1::Taylor0) 
cache1=a*b
"""
function mulT(a::Taylor0, b::Taylor0, cache1::Taylor0)
  for k in eachindex(a)
    @inbounds cache1[k] = a[0] * b[k]
    @inbounds for i = 1:k
      cache1[k] += a[i] * b[k-i]
    end
  end
  return cache1
end
function mulT(a::T, b::T, cache1::Taylor0) where {T<:Number}
  cache1[0] = a * b
  return cache1
end
#mul uses two caches when the op has many terms
function mulTT(a::Taylor0, b::T, cache1::Taylor0, cache2::Taylor0) where {T<:Number}#in middle ops: a is the returned cache1 so do not fudge a or cache1 before the end
  fill!(cache2.coeffs, b)
  @__dot__ cache1.coeffs = a.coeffs * cache2.coeffs  
  return cache1
end
mulTT(a::T, b::Taylor0, cache1::Taylor0, cache2::Taylor0) where {T<:Number} = mulTT(b, a, cache1, cache2)
function mulTT(a::Taylor0, b::Taylor0, cache1::Taylor0, cache2::Taylor0) #in middle ops: a is the returned cache1 so do not fudge a or cache1 before the end
  for k in eachindex(a)
    @inbounds cache2[k] = a[0] * b[k]
    @inbounds for i = 1:k
      cache2[k] += a[i] * b[k-i]
    end
  end
  @__dot__ cache1.coeffs = cache2.coeffs
  return cache1
end
function mulTT(a::T, b::T, cache1::Taylor0, cache2::Taylor0) where {T<:Number}
  #cache1 needs to be clean: a clear here will alloc...this is the only "mul" that wants a clean cache
  #sometimes the cache can be dirty at only first position: it will not throw an assert error!!!
  # cache[1].coeffs.=0.0 #it causes two allocs for mul(cst,cst,a,b)
  cache1[0] = a * b
  return cache1
end

"""
    muladdT(a::P,b::Q,c::R,cache1::Taylor0) where {P,Q,R <:Union{Taylor0,Number}}
cache1=a*b+c
"""
function muladdT(a::P, b::Q, c::R, cache1::Taylor0) where {P,Q,R<:Union{Taylor0,Number}}
  addT(mulT(a, b, cache1), c, cache1) 
end
"""
    mulsub(a::P,b::Q,c::R,cache1::Taylor0) where {P,Q,R <:Union{Taylor0,Number}}
cache1=a*b-c
"""
function mulsub(a::P, b::Q, c::R, cache1::Taylor0) where {P,Q,R<:Union{Taylor0,Number}}
  subT(mulT(a, b, cache1), c, cache1)
end
function divT(a::T, b::T, cache1::Taylor0) where {T<:Number}
  cache1[0] = a / b
  return cache1
end
function divT(a::Taylor0, b::T, cache1::Taylor0) where {T<:Number}
  fill!(cache1.coeffs, b)
  @__dot__ cache1.coeffs = a.coeffs / cache1.coeffs
  return cache1
end
function divT(a::T, b::Taylor0, cache1::Taylor0) where {T<:Number}
  powerT(b, -1.0, cache1)
  @__dot__ cache1.coeffs = cache1.coeffs * a 
  return cache1
end

"""
    divT(a::Taylor0, b::Taylor0,cache1::Taylor0) 
cache1=a/b
"""
function divT(a::Taylor0, b::Taylor0, cache1::Taylor0)
  iszero(a.coeffs) && !iszero(b.coeffs) && return cache1 
   ordfact, cdivfact = divfactorization(a, b)
   cache1[0] = cdivfact
  for k = 1:a.order-ordfact
    imin = max(0, k + ordfact - b.order)
    @inbounds cache1[k] = cache1[imin] * b[k+ordfact-imin]
    @inbounds for i = imin+1:k-1
      cache1[k] += cache1[i] * b[k+ordfact-i]
    end
    if k + ordfact ≤ b.order
      @inbounds cache1[k] = (a[k+ordfact] - cache1[k]) / b[ordfact]
    else
      @inbounds cache1[k] = -cache1[k] / b[ordfact]
    end
  end
  return cache1
end


function clearCache(cache::Vector{Taylor0}, ::Val{CS}, ::Val{2}) where {CS}
  for i = 1:CS
    cache[i][0] = 0.0
    cache[i][1] = 0.0
    cache[i][2] = 0.0
  end
end
function clearCache(cache::Vector{Taylor0}, ::Val{CS}, ::Val{1}) where {CS}
  for i = 1:CS
    cache[i][0] = 0.0
    cache[i][1] = 0.0
  end
end
