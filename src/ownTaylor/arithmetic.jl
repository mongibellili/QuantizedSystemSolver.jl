# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

# Arithmetic operations: +, -, *, /

## Equality ##
#= ==(a::Taylor0, b::Taylor0{S}) where {T<:Number,S<:Number} = ==(promote(a,b)...) =#

function ==(a::Taylor0, b::Taylor0) 
#=     if a.order != b.order
        a, b = fixorder(a, b)
    end =#
    return a.coeffs == b.coeffs
end
iszero(a::Taylor0) = iszero(a.coeffs)
## zero and one ##
 zero(a::Taylor0) = Taylor0(zero.(a.coeffs))
 function zero(a::Taylor0,cache::Taylor0)
    cache.coeffs .= 0.0
    return cache
 end
 function one(a::Taylor0,cache::Taylor0) 
    cache.coeffs .= 0.0
    cache[0]=1.0
    return cache
end
function one(a::Taylor0)
        b = zero(a)
        b[0] = one(b[0])
        return b

end
## Addition  fallback for case a+a+a+a+a+a+a+a+a+a+a+a+a+a##

#= (+)(a::Taylor0, b::Taylor0{S}) where {T<:Number,S<:Number} =
    (+)(promote(a,b)...) =#

function (+)(a::Taylor0, b::Taylor0) 
#=     if a.order != b.order
        a, b = fixorder(a, b)
    end =#
    v = similar(a.coeffs)
    @__dot__ v = (+)(a.coeffs, b.coeffs)
    return Taylor0(v, a.order)
end

function (+)(a::Taylor0)
    v = similar(a.coeffs)
    @__dot__ v = (+)(a.coeffs)
    return Taylor0(v, a.order)
end

#= (+)(a::Taylor0, b::S) where {T<:Number,S<:Number} =
(+)(promote(a,b)...) =#

function (+)(a::Taylor0, b::T) where {T<:Number}
    coeffs = copy(a.coeffs)
    @inbounds coeffs[1] = (+)(a[0], b)
    return Taylor0(coeffs, a.order)
end

#= (+)(b::S, a::Taylor0) where {T<:Number,S<:Number} =
    (+)(promote(b,a)...) =#

function (+)(b::T, a::Taylor0) where {T<:Number}
    coeffs = similar(a.coeffs)
    @__dot__ coeffs = (+)(a.coeffs)
    @inbounds coeffs[1] = (+)(b, a[0])
    return Taylor0(coeffs, a.order)
end

       
## substraction ##

#= (-)(a::Taylor0, b::Taylor0{S}) where {T<:Number,S<:Number} =
    (-)(promote(a,b)...) =#

function (-)(a::Taylor0, b::Taylor0) 
    if a.order != b.order
        a, b = fixorder(a, b)
    end
    v = similar(a.coeffs)
    @__dot__ v = (-)(a.coeffs, b.coeffs)
    return Taylor0(v, a.order)
end

function (-)(a::Taylor0)
    v = similar(a.coeffs)
    @__dot__ v = (-)(a.coeffs)
    return Taylor0(v, a.order)
end

#= (-)(a::Taylor0, b::S) where {T<:Number,S<:Number} =
(-)(promote(a,b)...) =#

function (-)(a::Taylor0, b::T) where {T<:Number}
    coeffs = copy(a.coeffs)
    @inbounds coeffs[1] = (-)(a[0], b)
    return Taylor0(coeffs, a.order)
end

#= (-)(b::S, a::Taylor0) where {T<:Number,S<:Number} =
    (-)(promote(b,a)...) =#

function (-)(b::T, a::Taylor0) where {T<:Number}
    coeffs = similar(a.coeffs)
    @__dot__ coeffs = (-)(a.coeffs)
    @inbounds coeffs[1] = (-)(b, a[0])
    return Taylor0(coeffs, a.order)
end       
    






## Multiplication fallback for case a*a*a*a*a*a*a*a*a*a*a*a*a##

function *(a::Taylor0, b::Taylor0) 
#=         if a.order != b.order
            a, b = fixorder(a, b)
        end =#
        c = Taylor0(zero(a[0]), a.order)
        for ord in eachindex(c)
            mul!(c, a, b, ord) # updates c[ord]
        end
        return c
end




#= function *(a::T, b::Taylor0{S}) where {T<:NumberNotSeries,S<:NumberNotSeries}
    @inbounds aux = a * b.coeffs[1]
    v = Array{typeof(aux)}(undef, length(b.coeffs))
    @__dot__ v = a * b.coeffs
    return Taylor0(v, b.order)
end

*(b::Taylor0{S}, a::T) where {T<:NumberNotSeries,S<:NumberNotSeries} = a * b =#

function (*)(a::T, b::Taylor0) where {T<:Number}
    v = Array{T}(undef, length(b.coeffs))
    @__dot__ v = a * b.coeffs
    return Taylor0(v, b.order)
end

*(b::Taylor0, a::T) where {T<:Number} = a * b
 

# Internal multiplication functions

     @inline function mul!(c::Taylor0, a::Taylor0, b::Taylor0, k::Int)        
            @inbounds c[k] = a[0] * b[k]       
        @inbounds for i = 1:k         
                c[k] += a[i] * b[k-i]         
        end
        return nothing
    end

#=     @inline function mul!(v::Taylor0, a::Taylor0, b::NumberNotSeries, k::Int)
        @inbounds v[k] = a[k] * b
        return nothing
    end
    @inline function mul!(v::Taylor0, a::NumberNotSeries, b::Taylor0, k::Int)
        @inbounds v[k] = a * b[k]
        return nothing
    end =#

#= 
@doc doc"""
    mul!(c, a, b, k::Int) --> nothing

Update the `k`-th expansion coefficient `c[k]` of `c = a * b`,
where all `c`, `a`, and `b` are either `Taylor0` or `TaylorN`.

The coefficients are given by

```math
c_k = \sum_{j=0}^k a_j b_{k-j}.
```

""" mul!


"""
    mul!(c, a, b) --> nothing

Return `c = a*b` with no allocation; all arguments are `HomogeneousPolynomial`.

"""
 =#


## Division ##
#= function /(a::Taylor0{Rational{T}}, b::S) where {T<:Integer,S<:NumberNotSeries}
    R = typeof( a[0] // b)
    v = Array{R}(undef, a.order+1)
    @__dot__ v = a.coeffs // b
    return Taylor0(v, a.order)
end =#

#= function /(a::Taylor0, b::S) where {T<:NumberNotSeries,S<:NumberNotSeries}
    @inbounds aux = a.coeffs[1] / b
    v = Array{typeof(aux)}(undef, length(a.coeffs))
    @__dot__ v = a.coeffs / b
    return Taylor0(v, a.order)
end =#

function /(a::Taylor0, b::T) where {T<:Number}
    @inbounds aux = a.coeffs[1] / b
    v = Array{typeof(aux)}(undef, length(a.coeffs))
    @__dot__ v = a.coeffs / b
    return Taylor0(v, a.order)
end

#= /(a::Taylor0, b::Taylor0{S}) where {T<:Number,S<:Number} = /(promote(a,b)...) =#

function /(a::Taylor0, b::Taylor0) 
    iszero(a) && !iszero(b) && return zero(a)
    if a.order != b.order
        a, b = fixorder(a, b)
    end
    # order and coefficient of first factorized term
    ordfact, cdivfact = divfactorization(a, b)
    c = Taylor0(cdivfact, a.order-ordfact)
    for ord in eachindex(c)
        div!(c, a, b, ord) # updates c[ord]
    end
    return c
end




@inline function divfactorization(a1::Taylor0, b1::Taylor0)
    # order of first factorized term; a1 and b1 assumed to be of the same order
    a1nz = findfirst(a1)
    b1nz = findfirst(b1)
    a1nz = a1nz ≥ 0 ? a1nz : a1.order
    b1nz = b1nz ≥ 0 ? b1nz : a1.order
    ordfact = min(a1nz, b1nz)
    cdivfact = a1[ordfact] / b1[ordfact]

    # Is the polynomial factorizable?
    iszero(b1[ordfact]) && throw( ArgumentError(
        """Division does not define a Taylor0 polynomial;
        order k=$(ordfact) => coeff[$(ordfact)]=$(cdivfact).""") )

    return ordfact, cdivfact
end


## TODO: Implement factorization (divfactorization) for TaylorN polynomials

#= 
# Homogeneous coefficient for the division
@doc doc"""
    div!(c, a, b, k::Int)

Compute the `k-th` expansion coefficient `c[k]` of `c = a / b`,
where all `c`, `a` and `b` are either `Taylor0` or `TaylorN`.

The coefficients are given by

```math
c_k =  \frac{1}{b_0} \big(a_k - \sum_{j=0}^{k-1} c_j b_{k-j}\big).
```

For `Taylor0` polynomials, a similar formula is implemented which
exploits `k_0`, the order of the first non-zero coefficient of `a`.
""" div! =#

@inline function div!(c::Taylor0, a::Taylor0, b::Taylor0, k::Int)

    # order and coefficient of first factorized term
    ordfact, cdivfact = divfactorization(a, b)
    if k == 0
        @inbounds c[0] = cdivfact
        return nothing
    end

    imin = max(0, k+ordfact-b.order)
    @inbounds c[k] = c[imin] * b[k+ordfact-imin]
    @inbounds for i = imin+1:k-1
        c[k] += c[i] * b[k+ordfact-i]
    end
    if k+ordfact ≤ b.order
        @inbounds c[k] = (a[k+ordfact]-c[k]) / b[ordfact]
    else
        @inbounds c[k] = - c[k] / b[ordfact]
    end
    return nothing
end

#= @inline function div!(v::Taylor0, a::Taylor0, b::NumberNotSeries, k::Int)
    @inbounds v[k] = a[k] / b
    return nothing
end

div!(v::Taylor0, b::NumberNotSeries, a::Taylor0, k::Int) =
    div!(v::Taylor0, Taylor0(b, a.order), a, k) =#





#= """
    mul!(Y, A, B)

Multiply A*B and save the result in Y.
""" =#
#= function mul!(y::Vector{Taylor0},
        a::Union{Matrix{T},SparseMatrixCSC{T}},
        b::Vector{Taylor0}) where {T<:Number}

    n, k = size(a)
    @assert (length(y)== n && length(b)== k)

    # determine the maximal order of b
    # order = maximum([b1.order for b1 in b])
    order = maximum(get_order.(b))

    # Use matrices of coefficients (of proper size) and mul!
    # B = zeros(T, k, order+1)
    B = Array{T}(undef, k, order+1)
    B = zero.(B)
    for i = 1:k
        @inbounds ord = b[i].order
        @inbounds for j = 1:ord+1
            B[i,j] = b[i][j-1]
        end
    end
    Y = Array{T}(undef, n, order+1)
    mul!(Y, a, B)
    @inbounds for i = 1:n
        # y[i] = Taylor0( collect(Y[i,:]), order)
        y[i] = Taylor0( Y[i,:], order)
    end

    return y
end =#


# Adapted from (Julia v1.2) stdlib/v1.2/LinearAlgebra/src/dense.jl#721-734,
# licensed under MIT "Expat".
# Specialize a method of `inv` for Matrix{Taylor0}. Simply, avoid pivoting,
# since the polynomial field is not an ordered one.
# function Base.inv(A::StridedMatrix{Taylor0}) where T
#     checksquare(A)
#     S = Taylor0{typeof((one(T)*zero(T) + one(T)*zero(T))/one(T))}
#     AA = convert(AbstractArray{S}, A)
#     if istriu(AA)
#         Ai = triu!(parent(inv(UpperTriangular(AA))))
#     elseif istril(AA)
#         Ai = tril!(parent(inv(LowerTriangular(AA))))
#     else
#         # Do not use pivoting !!
#         Ai = inv!(lu(AA, Val(false)))
#         Ai = convert(typeof(parent(Ai)), Ai)
#     end
#     return Ai
# end

#= # see https://github.com/JuliaLang/julia/pull/40623
const LU_RowMaximum = VERSION >= v"1.7.0-DEV.1188" ? RowMaximum() : Val(true)
const LU_NoPivot = VERSION >= v"1.7.0-DEV.1188" ? NoPivot() : Val(false)

# Adapted from (Julia v1.2) stdlib/v1.2/LinearAlgebra/src/lu.jl#240-253
# and (Julia v1.4.0-dev) stdlib/LinearAlgebra/v1.4/src/lu.jl#270-274,
# licensed under MIT "Expat".
# Specialize a method of `lu` for Matrix{Taylor0}, which avoids pivoting,
# since the polynomial field is not an ordered one.
# We can't assume an ordered field so we first try without pivoting
function lu(A::AbstractMatrix{Taylor0}; check::Bool = true) where {T<:Number}
    S = Taylor0{lutype(T)}
    F = lu!(copy_oftype(A, S), LU_NoPivot; check = false)
    if issuccess(F)
        return F
    else
        return lu!(copy_oftype(A, S), LU_RowMaximum; check = check)
    end
end =#
