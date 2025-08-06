
# ─────────────────────────────────────────────
# Linear Algebra Extensions for Vector{Taylor0}
# ─────────────────────────────────────────────

#constant_term(x::Taylor0) = x.coeffs[1]
constant_term(x::Real) = x


# Norm: vector of Taylor0 → take .coeffs[1] and compute norm
function norm(v::AbstractVector{Taylor0})
    return sqrt(sum(x.coeffs[1]^2 for x in v))
end


# Transpose: return matrix of first coeffs
function transpose(A::AbstractMatrix{Taylor0})
    return permutedims([x.coeffs[1] for x in A], (2, 1))
end

# Adjoint (Hermitian transpose): same as transpose for real coeffs
function adjoint(A::AbstractMatrix{Taylor0})
    return transpose(A)  # no conjugation needed for real numbers
end


function dot(u::AbstractVector{Taylor0}, v::AbstractVector{T}) where {T}
    sum(constant_term(x) * constant_term(y) for (x, y) in zip(u, v))
end

# Only missing case: vector of Real, vector of Taylor0
dot(u::AbstractVector{<:Real}, v::AbstractVector{Taylor0}) = dot(v, u)


function *(A::AbstractMatrix{Taylor0}, x::AbstractVector{T}) where {T}
    #allocates!!
    [sum(constant_term(A[i, j]) * constant_term(x[j]) for j in 1:size(A, 2)) for i in 1:size(A, 1)]
end

function *(A::AbstractMatrix{<:Real}, x::AbstractVector{Taylor0})
    #allocates
    [sum(constant_term(x[j]) * A[i, j] for j in 1:size(A, 2)) for i in 1:size(A, 1)]
end


function cross(u::AbstractVector{Taylor0}, v::AbstractVector{T}) where {T}
    @assert length(u) == 3 && length(v) == 3
    x1, y1, z1 = constant_term.(u)
    x2, y2, z2 = constant_term.(v)
    [
        y1*z2 - z1*y2,
        z1*x2 - x1*z2,
        x1*y2 - y1*x2,
    ]
end

cross(u::AbstractVector{<:Real}, v::AbstractVector{Taylor0}) = - cross(v, u)  # flipped
