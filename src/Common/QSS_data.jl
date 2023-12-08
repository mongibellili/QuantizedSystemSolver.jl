#hold helper datastructures needed for simulation, can be seen as the model in the qss architecture (model-integrator-quantizer)
struct CommonQSS_data{Z}
    quantum :: Vector{Float64} 
    x :: Vector{Taylor0}  #MVector cannot hold non-isbits
    q :: Vector{Taylor0}
    tx ::  Vector{Float64} 
    tq :: Vector{Float64} 
    d::Vector{Float64} 
    nextStateTime :: Vector{Float64}    
    nextInputTime :: Vector{Float64} # 
    nextEventTime :: MVector{Z,Float64}  
    t::Taylor0# mute taylor var to be used with math functions to represent time
    integratorCache::Taylor0
    taylorOpsCache::Vector{Taylor0}
    finalTime:: Float64   
    savetimeincrement::Float64 
    initialTime :: Float64    
    dQmin ::Float64    
    dQrel ::Float64  
    maxErr ::Float64  
    savedTimes :: Vector{Vector{Float64}}
    savedVars:: Vector{Vector{Float64}}
    
end
struct LiQSS_data{O,Sparsity}
   vs::Val{Sparsity}
   a::Vector{Vector{Float64}}
   # u:: Vector{Vector{MVector{O,Float64}}}
    qaux::Vector{MVector{O,Float64}}
    olddx::Vector{MVector{O,Float64}}
    dxaux::Vector{MVector{O,Float64}}
    olddxSpec::Vector{MVector{O,Float64}}
end
struct LightSpecialQSS_data{O1,Lightness}<:SpecialQSS_data{O1,Lightness}
    ls::Val{Lightness}
    p::Val{O1}
    savedVars :: Vector{Vector{Float64}} #has to be vector (not SA) cuz to be resized in integrator
   
end

#= struct HeavySpecialQSS_data{T,O1,Lightness}<:SpecialQSS_data{T,O1,Lightness}
    ls::Val{Lightness}
    savedVars :: Vector{Array{Taylor0}} #has to be vector (not SA) cuz to be resized in integrator
    prevStepVal ::MVector{T,MVector{O1,Float64}} #use 
end =#

struct SpecialLiQSS_data<:SpecialLiqssQSS_data
    cacheA::MVector{1,Float64}
    direction::Vector{Float64}
    qminus::Vector{Float64}
    buddySimul::MVector{2,Int}
    prevStepVal ::Vector{Float64} 
end


#= function createSpecialQSS_data(savedVars :: Vector{Array{Taylor0}}, prevStepVal::MVector{T,MVector{O1,Float64}},cacheA::MVector{1,Int}) where {T,O1}
    hv=HeavySpecialQSS_data(savedVars,prevStepVal,cacheA)
 end
 =#
 