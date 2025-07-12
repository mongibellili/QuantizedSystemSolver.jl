
"""
    CommonQSS_Data{Z}
 helper datastructures that hold settings, temporary saved data needed for simulation, and results. 
"""
 struct CommonQSS_Data{Z}
    quantum :: Vector{Float64} 
    x :: Vector{Taylor0}  #MVector cannot hold non-isbits
    q :: Vector{Taylor0}
    tx ::  Vector{Float64} 
    tq :: Vector{Float64} 
    d::Vector{Float64} 
    nextStateTime :: Vector{Float64}    
    nextInputTime :: Vector{Float64} # 
    nextEventTime :: MVector{Z,Float64}  
    t::Taylor0#  taylor var to be used with math functions to represent time
    taylorOpsCache::Vector{Taylor0}
    finalTime:: Float64    
    initialTime :: Float64    
    absQ ::Float64    
    relQ ::Float64  
    maxErr ::Float64  
    maxiters ::Int
    verbose::Bool
    savedTimes :: Vector{Vector{Float64}}
    savedVars:: Vector{Vector{Float64}}
end

"""
    AexprLiQSS_data{O,M} 
helper datastructures needed only for implicit case to store the Jacobian coefficients from an expression function.
The field variables are:
  - `cycleDetection::Val{M}` # cycle detection mechanism
  - `cacheA::MVector{1,Float64}` # the coef a is computed from a function and the result is saved in a chache.
  - `qaux::Vector{MVector{O,Float64}}` # to update u, q^- and dx^- are needed, so these 2 are to save the old values
  - `dxaux::Vector{MVector{O,Float64}}` # so these 2 are to save the old values
  - `olddx::Int64` # not needed for coefficient 'a' precomputed as expression, but kept for compatibility with parent type
"""
struct AexprLiQSS_data{O,M} <: LiQSS_Data{O,M} 
    cycleDetection::Val{M}               # cycle detection mechanism
    a::Int64                            # not needed, but kept for compatibility with parent type
    cacheA::MVector{1,Float64}      # the coef a is computed from a function and the result is saved in a chache.
    qaux::Vector{MVector{O,Float64}} # to update u, q^- and dx^- are needed
    dxaux::Vector{MVector{O,Float64}} # so these 2 are to save the old values
    olddx::Int64     # not needed for coefficient 'a' precomputed as expression, but kept for compatibility with parent type
end



"""
    AmanualLiQSS_data{O,M} 

A data structure for storing information specific to the manual implementation of the LiQSS (Linear Quantized State System) method.

# Type Parameters
- `O`: The order of the QSS method.
- `M`: The cycle detection mechanism number.

# Description
This struct extends `LiQSS_Data` and is intended for use in scenarios where a manual computation of the jacobian coefficient is recommended or required. It includes fields for storing the cycle detection mechanism, coefficient 'a', auxiliary vectors for quantized states and old values, and an additional vector for storing old dx values.

    The field variables are:
    - `cycleDetection::Val{M}`: Cycle detection mechanism.
    - `a::Vector{Vector{Float64}}`: Coefficient 'a' computed manually.
    - `cacheA::Int64`: Not needed when coefficient 'a' is computed manually, but kept for compatibility with parent type.
    - `qaux::Vector{MVector{O,Float64}}`: Auxiliary vectors to update `u`, `q^-`.
    - `dxaux::Vector{MVector{O,Float64}}`: Auxiliary vectors `dx^- to save the old values (and high order derivatives).
    - `olddx::Vector{MVector{1,Float64}}`: Needed to compute the coefficient 'a' manually, storing old dx values (only first derivative).   

# See Also
- [`LiQSS_Data`](@ref)
"""
struct AmanualLiQSS_data{O,M}<: LiQSS_Data{O,M}
    cycleDetection::Val{M}
    a::Vector{Vector{Float64}} # coefficient 'a' computed manually
    cacheA::Int64   # not needed when coefficient 'a' is computed manually, but kept for compatibility with parent type
    qaux::Vector{MVector{O,Float64}} # to update u, q^- and dx^- are needed
    dxaux::Vector{MVector{O,Float64}}# so these 2 are to save the old values
    olddx::Vector{MVector{1,Float64}} # needed to compute the coefficient 'a' manually
end

