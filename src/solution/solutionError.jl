
"""
    getError(sol::Sol{T,O},index::Int,f::Function) where{T,O}

calculates the relative error of one variable of the solution with respect to a function of an analytic solution.
```math
\\begin{equation*}
err=\\sqrt{\\frac{\\sum(sol_{index}-f_{index})^2}{\\sum(f_{index})^2}}
\\end{equation*}
```

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
"""
    getAverageError(sol::Sol{T,O},f::Vector{Function}) where{T,O}

calculates the average relative error of all variables of the solution with respect to a vector of functions of analytic solutions.
```math
\\begin{align*}
& err_{index}=\\sqrt{\\frac{\\sum(sol_{index}-f_{index})^2}{\\sum(f_{index})^2}}\\
& avgError=\\frac{\\sum_{index=1}^{T}err_{index}}{T}
\\end{align*}
```
"""
function getAverageError(sol::Sol{T,O},f::Vector{Function}) where{T,O}
  numPoints=length(sol.savedTimes[1])
  allErrors=0.0
  for index=1:T
    sumTrueSqr=0.0
    sumDiffSqr=0.0
    relerror=0.0
    for i = 1:numPoints-1 #each point is a taylor
      ts=f[index](sol.savedTimes[index][i])
      sumDiffSqr+=(sol.savedVars[index][i]-ts)*(sol.savedVars[index][i]-ts)
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
"""
    getErrorByRefs(sol::Sol{T,O},index::Int,solRef::Vector{Any}) where{T,O}

calculates the relative error of one variable of the solution with respect to a vector of values received from a reference solution.
```math
\\begin{equation*}
err=\\sqrt{\\frac{\\sum(sol_{index}-solRef_{index})^2}{\\sum(solRef_{index})^2}}
\\end{equation*}
```

"""
function getErrorByRefs(sol::Sol{T,O},index::Int,solRef::Vector{Any}) where{T,O}
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


"""
    getAverageErrorByRefs(sol::Sol{T,O},solRef::Vector{Any}) where{T,O}

calculates the average relative error of the solution with respect to a reference solution.
The relative error is calculated for each variable, and then it is averaged over all variables.
```math
\\begin{align*}
& err_{index}=\\sqrt{\\frac{\\sum(sol_{index}-solRef_{index})^2}{\\sum(solRef_{index})^2}}\\
& avgError=\\frac{\\sum_{index=1}^{T}err_{index}}{T}
\\end{align*}
```
"""
function getAverageErrorByRefs(sol::Sol{T,O},solRef::Vector{Any}) where{T,O}
  # sol and solRef are interpolated using the same time points
  numPoints=length(sol.savedTimes[1])
  allErrors=0.0
  for index=1:T
      sumTrueSqr=0.0
      sumDiffSqr=0.0
      relerror=0.0
      for i = 1:numPoints-1 #each point is a taylor
          ts=solRef[i][index]
          Ns=sol.savedVars[index][i]
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



