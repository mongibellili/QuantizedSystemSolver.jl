
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
    LiQSS_Data{O,M}
helper datastructures needed only for implicit case
The field variables are:
  - cd::Val{M} #
  - cacheA::MVector{1,Float64}
  - qaux::Vector{MVector{O,Float64}}
  - dxaux::Vector{MVector{O,Float64}}
"""
struct LiQSS_Data{O,M}
    cd::Val{M}               # cycle detection mechanism
    cacheA::MVector{1,Float64}      # the coef a is computed from a function and the result is saved in a chache.
    qaux::Vector{MVector{O,Float64}} # to update u, q^- and dx^- are needed
    dxaux::Vector{MVector{O,Float64}} # so these 2 are to save the old values
end


