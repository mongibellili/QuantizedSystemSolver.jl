function prepareAii(i::Int,j::Int,a::Int,exactA::FU, q::Vector{Taylor0}, d::Vector{Float64}, cacheA::MVector{1,Float64} , simt::Float64, clF::F) where {F,FU}
  cacheA[1] = 0.0
  exactA(q, d, cacheA, i, j, simt,clF)
  return cacheA[1] 
end

function prepareAii(i::Int,j::Int,a::Vector{Vector{Float64}},exactA::FU, q::Vector{Taylor0}, d::Vector{Float64}, cacheA::Int, simt::Float64, clF::F) where {F,FU}
  return a[i][j]
end

function integrateOlddx(::Val{O},i::Int,x::Vector{Taylor0},tx::Vector{Float64},simt::Float64,olddx::Int64) where{O}
        #jac entries as expression...no need for olddx to manually compute entry
end

function integrateOlddx(::Val{1},i::Int,x::Vector{Taylor0},tx::Vector{Float64},simt::Float64,olddx::Vector{MVector{1,Float64}})
  olddx[i][1]=x[i].coeffs[2]
end
function integrateOlddx(::Val{2},i::Int,x::Vector{Taylor0},tx::Vector{Float64},simt::Float64,olddx::Vector{MVector{1,Float64}})
  elapsed = simt - tx[i]
  olddx[i][1]=x[i].coeffs[2]+elapsed*x[i].coeffs[3]*2
end


function updateOtherApprox(k::Int,j::Int,x::Vector{Taylor0},q::Vector{Taylor0},a::Vector{Vector{Float64}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{1,Float64}},simt::Float64)where{O}
  diffQ=q[j][0]-qaux[j][1]
  if diffQ != 0.0
    a[k][j]=(x[k][1]-olddx[k][1])/diffQ
  else
    a[k][j]=0.0
  end
end

function updateOtherApprox(k::Int,j::Int,x::Vector{Taylor0},q::Vector{Taylor0},a::Int64,qaux::Vector{MVector{O,Float64}},olddx::Int64,simt::Float64)where{O}

end

function updateLinearApprox(i::Int,x::Vector{Taylor0},q::Vector{Taylor0},a::Vector{Vector{Float64}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{1,Float64}},simt::Float64)where{O}
  diffQ=q[i][0]-qaux[i][1]
   if diffQ != 0.0
      a[i][i]=(x[i][1]-olddx[i][1])/diffQ
  else
      a[i][i]=0.0
  end
  
  return nothing
end


function updateLinearApprox(i::Int,x::Vector{Taylor0},q::Vector{Taylor0},a::Int64,qaux::Vector{MVector{O,Float64}},olddx::Int64,simt::Float64)where{O}

end
"""
    integrateState(::Val{0}, x::Taylor0,elapsed::Float64)

does nothing: created for elapse-updating q in order1 which does not happen. This is needed in order to have one integrator function for all orders.
"""
@inline function integrateState(::Val{0}, x::Taylor0,elapsed::Float64) 
  #nothing: created for elapse-updating q in order1 which does not happen
end

"""
    integrateState(::Val{1}, x::Taylor0,elapsed::Float64)

Integrates the state for a first-order quantized system.

# Arguments
- `::Val{1}`: A type parameter indicating the order of the quantized system.
- `x::Taylor0`: The current state variable x represented as a `Taylor0` object.
- `elapsed::Float64`: The elapsed time since the last integration step of this variable x.


"""
@inline function integrateState(::Val{1}, x::Taylor0,elapsed::Float64) 
  x.coeffs[1] = x(elapsed)
end

"""
    integrateState(::Val{2}, x::Taylor0,elapsed::Float64)

Integrates a state variable x and its first derivative using a second-order Taylor series approximation

# Arguments
- `::Val{2}`: A type parameter indicating the order of the Taylor series (second-order in this case).
- `x::Taylor0`: The current state variable x represented as a Taylor0 object.
- `elapsed::Float64`: The elapsed time over which to integrate the state x.

"""
@inline function integrateState(::Val{2}, x::Taylor0,elapsed::Float64) 
  x.coeffs[1] = x(elapsed)
  x.coeffs[2] = x.coeffs[2]+elapsed*x.coeffs[3]*2
end


"""
    computeDerivative(::Val{1}, x::Taylor0, f::Taylor0)

copies the derivative from `f` to the first derivative of x for the first order.

# Arguments
- `::Val{1}`: A type parameter indicating the order of the derivative.
- `x::Taylor0`: The state variable .
- `f::Taylor0`: The Taylor series function that corresponds to the derivative.

"""
function computeDerivative(::Val{1}, x::Taylor0, f::Taylor0)
    x.coeffs[2] =f.coeffs[1]
    return nothing
end

"""
    computeDerivative(::Val{2}, x::Taylor0, f::Taylor0)

copies the first and second derivatives from `f` to the derivatives of x for the second-order quantizer.

# Arguments
- `::Val{2}`: A type parameter indicating the order of the quantizer.
- `x::Taylor0`: The state variable .
- `f::Taylor0`: The Taylor series function that corresponds to the derivative.

"""
function computeDerivative(::Val{2}, x::Taylor0, f::Taylor0) 
    x.coeffs[2] =f.coeffs[1]
    x.coeffs[3] =f.coeffs[2]/2
    return nothing
end




"""
    computeNextTime(::Val{1}, i::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0}, quantum::Vector{Float64})

Compute the next time of change for a given state variable.

# Arguments
- `::Val{1}`: A type parameter indicating a first order quantization method.
- `i::Int`: The index of the current state variable.
- `simt::Float64`: The current simulation time.
- `nextTime::Vector{Float64}`: A vector containing the next time values for each state variable.
- `x::Vector{Taylor0}`: A vector of Taylor series coefficients representing the state variables and their derivatives.
- `quantum::Vector{Float64}`: A vector containing the quantum values for the state variables.

# Returns
- The function updates the `nextTime` vector with the computed next time values for the state variables.
"""
function computeNextTime(::Val{1}, i::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0}, quantum::Vector{Float64})
  absDeltaT=1e-9 # minimum deltaT to protect against der=Inf coming from sqrt(0) for example...similar to min ΔQ
    if (x[i].coeffs[2]) != 0
        tempTime=max(abs(quantum[i] /(x[i].coeffs[2])),absDeltaT)# i can avoid the use of max
        if tempTime!=absDeltaT #normal
            nextTime[i] = simt + tempTime#sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
            #println("nextTime=",nextTime, "x: ", x)
        else#usual (quant/der) is very small
          x[i].coeffs[2]=sign(x[i].coeffs[2])*(abs(quantum[i])/absDeltaT)# adjust  derivative if it is too high
          nextTime[i] = simt + tempTime
          println("tempTime=",tempTime)
        end
    else
      nextTime[i] = Inf
    end
    return nothing
end

"""
    computeNextTime(::Val{2}, i::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0}, quantum::Vector{Float64})

Compute the next time for a given state variable `i` in a second-order quantized system.

# Arguments
- `::Val{2}`: A type parameter indicating the order of the quantized system (second-order).
- `i::Int`: The index of variable for which the next time is being computed.
- `simt::Float64`: The current simulation time.
- `nextTime::Vector{Float64}`: A vector containing the next times for all variables.
- `x::Vector{Taylor0}`: A vector of Taylor series coefficients representing the state variables and their derivatives.
- `quantum::Vector{Float64}`: A vector of quantum values for the state variables.

# Returns
- `Nothing`: This function updates the `nextTime` vector in place.

"""
function computeNextTime(::Val{2}, i::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0}, quantum::Vector{Float64})
    absDeltaT=1e-12 # minimum deltaT to protect against der=Inf coming from sqrt(0) for example...similar to min ΔQ
      if (x[i].coeffs[3]) != 0
          tempTime=max(sqrt(abs(quantum[i] / ((x[i].coeffs[3])))),absDeltaT)
          if tempTime!=absDeltaT #normal
              nextTime[i] = simt + tempTime#sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
          else#usual sqrt(quant/der) is very small
            x[i].coeffs[3]=sign(x[i].coeffs[3])*(abs(quantum[i])/(absDeltaT*absDeltaT))/2# adjust second derivative if it is too high
            nextTime[i] = simt + tempTime
          end
      else
        if (x[i].coeffs[2]) != 0
          tempTime=max(abs(quantum[i] /(x[i].coeffs[2])),absDeltaT)# i can avoid the use of max
          if tempTime!=absDeltaT #normal
              nextTime[i] = simt + tempTime#sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
          else#usual (quant/der) is very small
            x[i].coeffs[2]=sign(x[i].coeffs[2])*(abs(quantum[i])/absDeltaT)# adjust  derivative if it is too high
            nextTime[i] = simt + tempTime
          end
        else
          nextTime[i] = Inf
        end
      end
      return nothing 
end

"""
    reComputeNextTime(::Val{1}, index::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64})

Recomputes the next time for a given state variable `i` in a first-order quantized system after the derivative has changed.
similar to [`computeNextTime`](@ref) but it also account for the first derivative change.

"""
function reComputeNextTime(::Val{1}, index::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64})
  absDeltaT=1e-12
  if abs(q[index].coeffs[1] - (x[index].coeffs[1])) >= quantum[index] # this happened when var i and j s turns are now...var i depends on j, j is asked here for next time
    nextTime[index] = simt+1e-16
  else
    time1 =  minPosRoot(q[index].coeffs[1] - (x[index].coeffs[1]) - quantum[index],-x[index].coeffs[2], Val(1))
    time2 =  minPosRoot(q[index].coeffs[1] - (x[index].coeffs[1]) + quantum[index],-x[index].coeffs[2], Val(1))
    timeTemp = time1 < time2 ? time1 : time2
    tempTime=max(timeTemp,absDeltaT) #guard against very small Δt 
    nextTime[index] = simt +tempTime
  end
end

"""
    reComputeNextTime(::Val{2}, index::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0}, q::Vector{Taylor0}, quantum::Vector{Float64})

Recomputes the next time for a given state variable `i` in a second-order quantized system after the derivatives have changed.
similar to [`computeNextTime`](@ref) but it also account for the first and second derivatives changes.

"""
function reComputeNextTime(::Val{2}, index::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0}, q::Vector{Taylor0}, quantum::Vector{Float64})
  absDeltaT=1e-15
  if abs(q[index].coeffs[1] - (x[index].coeffs[1])) >= quantum[index] # this happened when var i and j s turns are now...var i depends on j, j is asked here for next time
    nextTime[index] = simt+1e-12
  else
    time1 =  minPosRoot(q[index].coeffs[1] - (x[index].coeffs[1]) - quantum[index], q[index].coeffs[2]-x[index].coeffs[2],-x[index].coeffs[3], Val(2))
    time2 =  minPosRoot(q[index].coeffs[1] - (x[index].coeffs[1]) + quantum[index], q[index].coeffs[2]-x[index].coeffs[2],-x[index].coeffs[3], Val(2))
    timeTemp = time1 < time2 ? time1 : time2
    tempTime=max(timeTemp,absDeltaT)#guard against very small Δt 
    nextTime[index] = simt +tempTime
  end
  return nothing
end
  

"""
    computeNextInputTime(::Val{1}, i::Int, t::Taylor0,f::F,clF::FF,d::Vector{Float64}, taylorOpsCache::Vector{Taylor0} ,nextInputTime::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64}) where {F,FF}

Compute the next input time for a given state variable in a first order method. This is needed when the differential equation depends on time only (i.e. does not depend on other state variables). It uses a prediction of the derivatives.

# Arguments
- `::Val{1}`: A type parameter indicating the specific method to use.
- `i::Int`: The index of the current state variable.
- `simt::Float64`: The current simulation time.
- `elapsed::Float64`: The elapsed time since the last update.
- `tt::Taylor0`: The Taylor series expansion of the state variable in a small time advance.
- `nextInputTime::Vector{Float64}`: A vector to store the computed next input times.
- `x::Vector{Taylor0}`: A vector of Taylor series expansions of the state variables.
- `quantum::Vector{Float64}`: A vector of quantum values for the state variables.

# Returns
- Nothing, it updates the `nextInputTime` vector with the next input times for the state variables.
"""
function computeNextInputTime(::Val{1}, i::Int, t::Taylor0,f::F,clF::FF,d::Vector{Float64}, taylorOpsCache::Vector{Taylor0} ,nextInputTime::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64}) where {F,FF}
  df=0.0
  oldDerX=x[i].coeffs[2]
  simt=t[0]
    if oldDerX!=0 #  
      tempstep=abs(1*quantum[i] / oldDerX)
      t[0]=simt+tempstep
      f(i, q, d,t,taylorOpsCache,clF)
      predictedDer=taylorOpsCache[1].coeffs[1]  
      if abs(predictedDer-oldDerX)>abs(predictedDer+oldDerX)/2 #significant change-->bad prediction --> lower stepsize
        nextInputTime[i]=simt+tempstep/100  
      else
        nextInputTime[i]=simt+tempstep
      end
     # @show simt,t[0],predictedDer,oldDerX, nextInputTime[i]
    else # think over this case
      t[0]=simt+1e-3
      f(i, q, d,t,taylorOpsCache,clF)   
      predictedDer=taylorOpsCache[1].coeffs[1]  
      if predictedDer!=0.0
        nextInputTime[i]=simt+abs(1*quantum[i] / 10.0*predictedDer)
      else
        nextInputTime[i] = Inf
      end
      #@show "older==0 ",simt,t[0],predictedDer, nextInputTime[i]
    end
    t[0]=simt
  
  return nothing
end

function discrete_computeNextInputTime(::Val{1}, i::Int, t::Taylor0,f::F,clF::FF,d::Vector{Float64}, taylorOpsCache::Vector{Taylor0} ,nextInputTime::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64}) where {F,FF}
  df=0.0
  oldDerX=x[i].coeffs[2]
  simt=t[0]
    if oldDerX!=0 #
      tempstep=abs(1*quantum[i] / oldDerX)
      t[0]=simt+tempstep
      f(i,-1, -1, q, d,t,taylorOpsCache,clF)
      predictedDer=taylorOpsCache[1].coeffs[1]  
      if abs(predictedDer-oldDerX)>abs(predictedDer+oldDerX)/2 #significant change-->bad prediction --> lower stepsize
        nextInputTime[i]=simt+tempstep/100  
      else
        nextInputTime[i]=simt+tempstep
      end
    
    else # think over this case
      t[0]=simt+1e-3
      f(i,-1, -1, q, d,t,taylorOpsCache,clF)
      predictedDer=taylorOpsCache[1].coeffs[1]  
      if predictedDer!=0.0
        nextInputTime[i]=simt+abs(1*quantum[i] / 10.0*predictedDer)
      else
        nextInputTime[i] = Inf
      end
    end
    t[0]=simt
  
  return nothing
end

"""
    computeNextInputTime(::Val{2}, i::Int, t::Taylor0,f::F,clF::FF,d::Vector{Float64}, taylorOpsCache::Vector{Taylor0} ,nextInputTime::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64}) where {F,FF}

Compute the next input time for a given state variable in a second order method. This is needed when the differential equation depends on time only (i.e. does not depend on other state variables). It uses a prediction of the derivatives.

"""
function computeNextInputTime(::Val{2}, i::Int, t::Taylor0,f::F,clF::FF,d::Vector{Float64}, taylorOpsCache::Vector{Taylor0} ,nextInputTime::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64}) where {F,FF}
  df=0.0
  oldDerX=x[i].coeffs[3]*2.0

  simt=t[0]
    if oldDerX!=0 #
      
      tempstep=sqrt(abs(2*quantum[i] / oldDerX))
      
      t[0]=simt+tempstep
      f(i, q, d,t,taylorOpsCache,clF)
      
      predictedDer=taylorOpsCache[1].coeffs[2]*2 

      if abs(predictedDer-oldDerX)>abs(predictedDer+oldDerX)/2 #significant change-->bad prediction --> lower stepsize
        nextInputTime[i]=simt+tempstep/10
      else
        nextInputTime[i]=simt+tempstep
      end
    else # think over this case
      t[0]=simt+1e-3
      f(i, q,  d,t,taylorOpsCache,clF) 
      
      predictedDer=taylorOpsCache[1].coeffs[1]  
      if predictedDer!=0.0
        nextInputTime[i]=simt+sqrt(abs(2*quantum[i] / 1.0*predictedDer))
      else
        nextInputTime[i] = Inf
      end
    end
    t[0]=simt
  
  return nothing
end
function discrete_computeNextInputTime(::Val{2}, i::Int, t::Taylor0,f::F,clF::FF,d::Vector{Float64}, taylorOpsCache::Vector{Taylor0} ,nextInputTime::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64}) where {F,FF}
  df=0.0
  oldDerX=x[i].coeffs[3]*2.0
  simt=t[0]
    if oldDerX!=0 #
      tempstep=sqrt(abs(2*quantum[i] / oldDerX))
      t[0]=simt+tempstep
      f(i, -1, -1, q, d, t, taylorOpsCache,clF)
      predictedDer=taylorOpsCache[1].coeffs[2]*2 
     # @show predictedDer
      if abs(predictedDer-oldDerX)>abs(predictedDer+oldDerX)/2 #significant change-->bad prediction --> lower stepsize
        nextInputTime[i]=simt+tempstep/10
      else
        nextInputTime[i]=simt+tempstep
      end
    else # think over this case
      t[0]=simt+1e-3
      f(i,-1, -1, q, d,t,taylorOpsCache,clF) 
      predictedDer=taylorOpsCache[1].coeffs[1]  
      if predictedDer!=0.0
        nextInputTime[i]=simt+sqrt(abs(2*quantum[i] / 1.0*predictedDer))
      else
        nextInputTime[i] = Inf
      end
     # @show "older==0 ",simt,t[0],predictedDer, nextInputTime[i]
    end
    t[0]=simt
  
  return nothing
end


"""
    computeNextEventTime(::Val{O},j::Int,ZCFun::Taylor0,oldsignValue::MMatrix{Z,2} ,simt::Float64,  nextEventTime :: MVector{Z,Float64}, quantum::Vector{Float64},absQ::Float64) where {O, Z}

Compute the next event time for a given zero-crossing function.

# Arguments
- `::Val{O}`: A type parameter indicating the order of the quantizer.
- `j::Int`: The index of the zero-crossing function being processed.
- `ZCFun::Taylor0`: the value of the zero-crossing function of type `Taylor0`.
- `oldsignValue::MMatrix{Z,2}`: The previous sign and value of the zero-crossing function.
- `simt::Float64`: The current simulation time.
- `nextEventTime:: MVector{Z,Float64}`: Vector contains the next event time for all zero-crossing functions.
- `quantum::Vector{Float64}`: A vector of quantum values for the state variables.
- `absQ::Float64`: The absolute quantum value.

# Returns
- The computed next event time for the state variable `j`.
"""
function computeNextEventTime(::Val{O},j::Int,ZCFun::Taylor0,oldsignValue::MMatrix{Z,2} ,simt::Float64,  nextEventTime :: MVector{Z,Float64}, quantum::Vector{Float64},absQ::Float64) where {O, Z}
 
  if oldsignValue[j,2] ==0.0 && ZCFun[0] ==0.0  # initial value 0 --> do not do anything
    nextEventTime[j]=Inf
  elseif (oldsignValue[j,1] != sign(ZCFun[0])) && abs(oldsignValue[j,2]) >1e-9*absQ #prevent double tapping: when zcf is leaving zero it should be considered an event
    nextEventTime[j]=simt 
    #@show simt,oldsignValue[j,1], sign(ZCFun[0])
  else # old and new ZCF both pos or both neg
    mpr=minPosRoot(ZCFun, Val(O)) 
    if mpr<1e-13 # prevent very close events
      mpr=1e-10
     end
     #@show simt,mpr
    nextEventTime[j] =simt + mpr
    oldsignValue[j,1]=sign(ZCFun[0])#update the values
    oldsignValue[j,2]=ZCFun[0]
  end
  
end
