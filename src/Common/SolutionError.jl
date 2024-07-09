
"""getError(sol::Sol{T,O},index::Int,f::Function) where{T,O}

  This function calculates the relative error of the solution with respect to a reference function.
"""
function getError(sol::Sol{T,O},index::Int,f::Function) where{T,O}
  index>T &&  error("the system contains only $T variables!")
  numPoints=length(sol.savedTimes[index])
  sumTrueSqr=0.0
  sumDiffSqr=0.0
  relerror=0.0
  for i = 1:numPoints-1 #each point is a taylor
    ts=f(sol.savedTimes[index][i])
    sumDiffSqr+=(sol.savedVars[index][i]-ts)*(sol.savedVars[index][i]-ts)
    sumTrueSqr+=ts*ts
  end
  relerror=sqrt(sumDiffSqr/sumTrueSqr)
  return relerror
end

function getErrorByRefs(sol::Sol{T,O},index::Int,solRef::Vector{Any})where{T,O}
  #numVars=length(sol.savedVars)
  index>T &&  error("the system contains only $T variables!")
  numPoints=length(sol.savedTimes[index])
  sumTrueSqr=0.0
  sumDiffSqr=0.0
  relerror=0.0
  for i = 1:numPoints-1 #-1 because savedtimes holds init cond at begining
    ts=solRef[i][index]
    sumDiffSqr+=(sol.savedVars[index][i]-ts)*(sol.savedVars[index][i]-ts)
    sumTrueSqr+=ts*ts
  end
  
  if  abs(sumTrueSqr)>1e-12
    relerror=sqrt(sumDiffSqr/sumTrueSqr)
    else
      relerror=0.0
    end
  return relerror
end

#= @inline function getX_fromSavedVars(savedVars :: Vector{Array{Taylor0}},index::Int,i::Int)
  return savedVars[index][i].coeffs[1]
end =#
@inline function getX_fromSavedVars(savedVars :: Vector{Vector{Float64}},index::Int,i::Int)
  return savedVars[index][i]
end

"""getAverageErrorByRefs(sol::Sol{T,O},solRef::Vector{Any}) where{T,O}

  This function calculates the average relative error of the solution with respect to a reference solution.
  The relative error is calculated for each variable, and then it is averaged over all variables.
"""
function getAverageErrorByRefs(sol::Sol{T,O},solRef::Vector{Any}) where{T,O}
  numPoints=length(sol.savedTimes[1])
  allErrors=0.0
  for index=1:T
      sumTrueSqr=0.0
      sumDiffSqr=0.0
      relerror=0.0
      for i = 1:numPoints-1 #each point is a taylor
          ts=solRef[i][index]
          Ns=getX_fromSavedVars(sol.savedVars,index,i)
          sumDiffSqr+=(Ns-ts)*(Ns-ts)
          sumTrueSqr+=ts*ts
      end
      if  abs(sumTrueSqr)>1e-12
      relerror=sqrt(sumDiffSqr/sumTrueSqr)
      else
        relerror=0.0
      end
      
      allErrors+= relerror
  end
  return allErrors/T
end



