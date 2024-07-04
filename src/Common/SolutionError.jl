
function getError(sol::Sol{T,O},index::Int,f::Function)where{T,O}
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

function getErrorByRefs(solRef::Vector{Any},solmliqss::Sol{T,O},index::Int)where{T,O}
  #numVars=length(solmliqss.savedVars)
  index>T &&  error("the system contains only $T variables!")
  numPoints=length(solmliqss.savedTimes[index])
  sumTrueSqr=0.0
  sumDiffSqr=0.0
  relerror=0.0
  for i = 1:numPoints-1 #-1 because savedtimes holds init cond at begining
    ts=solRef[i][index]
    sumDiffSqr+=(solmliqss.savedVars[index][i]-ts)*(solmliqss.savedVars[index][i]-ts)
    sumTrueSqr+=ts*ts
  end
  
  if  abs(sumTrueSqr)>1e-12
    relerror=sqrt(sumDiffSqr/sumTrueSqr)
    else
      relerror=0.0
    end
  return relerror
end


#= function getAllErrorsByRefs(solRef::Vector{Any},solmliqss::Sol{T,O})where{T,O}
  numPoints=length(solmliqss.savedTimes[1])
  allErrors=Array{Float64}(undef, T)
  for index=1:T
      sumTrueSqr=0.0
      sumDiffSqr=0.0
      relerror=0.0
      for i = 1:numPoints-1 #each point is a taylor
          ts=solRef[i][index]
          Ns=getX_fromSavedVars(solmliqss.savedVars,index,i)
          sumDiffSqr+=(Ns-ts)*(Ns-ts)
          sumTrueSqr+=ts*ts
      end
      relerror=sqrt(sumDiffSqr/sumTrueSqr)
      
      allErrors[index]= relerror
  end
  return allErrors
end =#
@inline function getX_fromSavedVars(savedVars :: Vector{Array{Taylor0}},index::Int,i::Int)
  return savedVars[index][i].coeffs[1]
end
@inline function getX_fromSavedVars(savedVars :: Vector{Vector{Float64}},index::Int,i::Int)
  return savedVars[index][i]
end

function getAverageErrorByRefs(solRef::Vector{Any},solmliqss::Sol{T,O})where{T,O}
  numPoints=length(solmliqss.savedTimes[1])
  allErrors=0.0
  for index=1:T
      sumTrueSqr=0.0
      sumDiffSqr=0.0
      relerror=0.0
      for i = 1:numPoints-1 #each point is a taylor
          ts=solRef[i][index]
          Ns=getX_fromSavedVars(solmliqss.savedVars,index,i)
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



#= 
function plotRelativeError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  if index<=numVars
    temp = []
    tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      if ft>1e-12 || ft < -1e-12
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)/ft))
        push!(tempt,sol.savedTimes[i])
      end
    end
    display(plot(tempt, temp,title="RelError:(S-T)/T for x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200),ylims=(0.0,0.0002) =#) )
    println("press enter to exit")
    readline() 
  else
    error("the system contains only $numVars variables!")
  end
end



function saveRelativeError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)

  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  if index<=numVars
    temp = []
    tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      if ft>1e-12 || ft < -1e-12
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)/ft))
        push!(tempt,sol.savedTimes[i])
      end
    end
  # display(plot!(tempt, temp,title="RelError:(S-T)/T for x$index",label="$(sol.algName)")) 
  p1=plot(tempt, temp,title="RelError:(S-T)/T for x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200),ylims=(0.0,0.0002) =#) 
  savefig(p1, "relError_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_x$(index)_$(timestamp).png")
  else
    error("the system contains only $numVars variables!")
  end
end

function plotAbsoluteError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  if index<=numVars
    temp = []
   # tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      #if ft>1e-2 || ft < -1e-2
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)))
       # push!(tempt,sol.savedTimes[i])
     # end
    end
   #display(plot!(sol.savedTimes, temp,title="AbsError:(S-T) for x$index",label="$(sol.algName)")) 
  display(plot(sol.savedTimes, temp,title="AbsError:(S-T) for x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200),ylims=(0.0,0.0002) =#))
  
    println("press enter to exit")
    readline() 
  else
    error("the system contains only $numVars variables!")
  end
end



function saveAbsoluteError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
 
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  if index<=numVars
    temp = []
   # tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      #if ft>1e-2 || ft < -1e-2
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)))
       # push!(tempt,sol.savedTimes[i])
     # end
    end
   #display(plot!(sol.savedTimes, temp,title="AbsError:(S-T) for x$index",label="$(sol.algName)")) 
   p1=plot(sol.savedTimes, temp,title="AbsError:(S-T) for x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200),ylims=(0.0,0.0002) =#)
   savefig(p1, "absError_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_x$(index)_$(timestamp).png")
     
  else
    error("the system contains only $numVars variables!")
  end
end

function plotCumulativeSquaredRelativeError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  sumTrueSqr=0.0
  sumDiffSqr=0.0
  
  if index<=numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      sumDiffSqr+=(sol.savedVars[index][i].coeffs[1]-ft)*(sol.savedVars[index][i].coeffs[1]-ft)
      sumTrueSqr+=ft*ft
        push!(temp, sqrt(sumDiffSqr/sumTrueSqr))
    end
   display(plot!(sol.savedTimes, temp,title="Error:sqrt(∑(S-T)^2/∑T^2) for x$index",label="$(sol.algName)")) 
  else
    error("the system contains only $numVars variables!")
  end
  println("press enter to exit")
  readline()
end
function plotMSE(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
 # sumTrueSqr=0.0
  sumDiffSqr=0.0
  
  if index<=numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      sumDiffSqr+=(sol.savedVars[index][i].coeffs[1]-ft)*(sol.savedVars[index][i].coeffs[1]-ft)
     # sumTrueSqr+=ft*ft
        push!(temp, (sumDiffSqr/i))
    end
   display(plot!(sol.savedTimes, temp,title="Error:(∑(S-T)^2/i) for x$index",label="$(sol.algName)")) 
  else
    error("the system contains only $numVars variables!")
  end
  println("press enter to exit")
  readline()
end
 =#