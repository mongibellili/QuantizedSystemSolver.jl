"""
    integrate(alg::QSSAlgorithm{:qss,O}, commonQssData::CommonQSS_Data{Z}, odep::ODEProblemData{JACMODE,T,D,Z,CS,F,JAC,CLS}, f::Function, jac::Function, SD::Function) where {JACMODE,T,D,Z,CS,F,JAC,CLS}

Integrates a nonlinear ordinary differential equation (ODE) problem with `events` using a Quantized State System (QSS) algorithm.

# Arguments
- `Al::QSSAlgorithm{:qss,O}`: The QSS algorithm to be used for integration.
- `commonQssData::CommonQSS_Data{Z}`: Common data structure for QSS algorithms.
- `odep::ODEProblemData{JACMODE,T,D,Z,CS,F,JAC,CLS}`: The nonlinear ODE problem to be solved.
- `f::Function`: The function defining the ODE system.
- `jac::Function`: The Jacobian dependency function of the ODE system.
- `SD::Function`: The state derivative dependency function.


# Type Parameters
- `JACMODE`: The type of the problem.
- `O`: The order of the QSS algorithm.
- `T`: The number of continuous variables.
- `Z`: The number of zero crossing functions.
- `D`: The number of discrete variables
- `CS`: The cache size.

# Returns
- A solution
"""
function integrate(alg::QSSAlgorithm{:qss,O}, commonQssData::CommonQSS_Data{Z,PT}, odep::ODEProblemData{JACMODE,T,D,Z,CS,F,JAC,CLS}, f::F, jac::JAC_CLS, SD::SD_CLS) where {O,JACMODE,T,D,Z,CS,F,JAC,CLS,PT,JAC_CLS,SD_CLS}

  VERBOSE,ft,initTime,relQ,absQ,relZ,absZ,maxErr,maxiters,quantum,nextStateTime,nextEventTime,nextInputTime,tx,tq,x,q,t,savedVars,savedTimes,
  taylorOpsCache,d,numStateSteps,numInputSteps,oldsignValue,clF,zc_SimpleJac,HZ,HD,SZ,evDep = initIntegrator(Val(O), Val(T), Val(CS), odep, commonQssData, f)

   f(1, -1, -1, q, d, t, taylorOpsCache,clF) # call once (for performance)
  
  for i = 1:T
    numStateSteps[i] = 0
    numInputSteps[i] = 0
    push!(savedVars[i], x[i][0])
    push!(savedTimes[i], initTime)
    quantum[i] = 2*relQ * abs(x[i].coeffs[1])
    quantum[i] = quantum[i] < absQ ? 2*absQ : quantum[i]
    quantum[i] = quantum[i] > maxErr ? maxErr : quantum[i]
    if isempty(jac(i))
        computeNextInputTime(Val(O), i, t, f,clF,d, taylorOpsCache, nextInputTime, x,q, quantum) 
    else
      computeNextTime(Val(O), i, initTime, nextStateTime, x, quantum)
    end
  end

  simt = initTime
  totalSteps = 0;inpuStepCount=0;simulStepCount = 0;evCount = 0;rejectedEvCount = 0
  
  while t[0] < ft && totalSteps < maxiters
    if totalSteps == maxiters - 1 @warn("The algorithm qss$O reached max iterations. The simulation will be stopped. Consider using a different algorithm!") end
    sch = updateScheduler(Val(T), nextStateTime, nextEventTime, nextInputTime)
    i = sch[1]; simt = sch[2]; stepType = sch[3]
    simt > ft && break  
    totalSteps += 1
    t[0] = simt
    if stepType == :ST_STATE
      #statestep += 1
      numStateSteps[i] += 1
      integrateState(Val(O), x[i], tx[i], simt) ;tx[i]=simt
      quantum[i] = 2*relQ * abs(x[i].coeffs[1]); quantum[i] = quantum[i] < absQ ? 2*absQ : quantum[i]; quantum[i] = quantum[i] > maxErr ? maxErr : quantum[i]
      for k = 1:O
        q[i].coeffs[k] = x[i].coeffs[k]
      end
      tq[i] = simt
      computeNextTime(Val(O), i, simt, nextStateTime, x, quantum)
      for j in (SD(i))
        elapsedx = simt - tx[j]
        if elapsedx > 0           x[j].coeffs[1] = x[j](elapsedx);          tx[j] = simt        end # 
        integrateState(Val(O - 1),q[j],tq[j],simt) ; tq[j]=simt 
        computeDerivatives(Val(O),Val(CS),f,j,q,tq,simt,d,t,taylorOpsCache,clF,x,jac)
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end
      for j in (SZ[i])
        updateZCF(Val(O),Val(CS),f,j,q,tq,simt,d,t,taylorOpsCache,clF,zc_SimpleJac)
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absZ, relZ)
      end
    elseif stepType == :ST_INPUT
      inpuStepCount+=1
      numInputSteps[i] += 1
      integrateState(Val(O), x[i], tx[i], simt) ;tx[i] = simt
      quantum[i] = 2*relQ * abs(x[i].coeffs[1]); quantum[i] = quantum[i] < absQ ? 2*absQ : quantum[i]; quantum[i] = quantum[i] > maxErr ? maxErr : quantum[i]
      for k = 1:O
        q[i].coeffs[k] = x[i].coeffs[k]
      end
      tq[i] = simt
      clearCache(taylorOpsCache, Val(CS), Val(O)); f(i, -1, -1, q, d, t, taylorOpsCache,clF)
      computeDerivative(Val(O), x[i], taylorOpsCache[1])
      computeNextInputTime(Val(O), i, t, f,clF,d, taylorOpsCache, nextInputTime, x,q, quantum)
 
  
      for j in (SD(i))
        elapsedx = simt - tx[j]
        if elapsedx > 0 
          x[j].coeffs[1] = x[j](elapsedx) # der later
          tx[j] = simt 
        end
        integrateState(Val(O - 1),q[j],tq[j],simt)  ;tq[j] = simt #q needs to be updated here for recomputeNext                 
        computeDerivatives(Val(O),Val(CS),f,j,q,tq,simt,d,t,taylorOpsCache,clF,x,jac)
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end
      for j in (SZ[i])
        updateZCF(Val(O),Val(CS),f,j,q,tq,simt,d,t,taylorOpsCache,clF,zc_SimpleJac)
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absZ,relZ)
      end
    else
        updateZCF(Val(O),Val(CS),f,i,q,tq,simt,d,t,taylorOpsCache,clF,zc_SimpleJac)


      if oldsignValue[i, 2] * taylorOpsCache[1][0] >= 0 # if computeNextEvent errored 
        #tol_zcf= absZ + relZ*abs(taylorOpsCache[1][0] )
        if abs(taylorOpsCache[1][0]) > absZ # if error is negligeable then ok consider as event, else reject....if both have same sign and zcf is not very small: zc==1e-9*absQ is allowed as an event
          
          oldsignValue[i, 2] = taylorOpsCache[1][0]; oldsignValue[i, 1] = sign(taylorOpsCache[1][0])
          computeNextEventTime(Val(O), i, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absZ,relZ)
          rejectedEvCount+=1
          continue #event rejected
        end
        if abs(oldsignValue[i, 2]) < absZ#*1e-6 # if oldsignValue very small, reject event , just already handled          
          rejectedEvCount+=1
          oldsignValue[i, 2] = taylorOpsCache[1][0]; oldsignValue[i, 1] = sign(taylorOpsCache[1][0])
          computeNextEventTime(Val(O), i, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absZ,relZ)
          continue #event rejected
        end
      end   
  
      if taylorOpsCache[1][0] > oldsignValue[i, 2] #scheduled rise
        modifiedIndex = 2 * i - 1
      elseif taylorOpsCache[1][0] < oldsignValue[i, 2] #scheduled drop
        modifiedIndex = 2 * i
      else # == ( zcf==oldZCF)
        nextEventTime[i]=Inf
        rejectedEvCount += 1
        continue
      end
      evCount += 1
      oldsignValue[i, 2] = taylorOpsCache[1][0]; oldsignValue[i, 1] = sign(taylorOpsCache[1][0])
      for b in evDep[modifiedIndex].evContRHS
          integrateState(Val(O - 1),q[b],tq[b],simt) ;tq[b] = simt 
      end
     #=  for k in zc_SimpleJac[i]  
        push!(savedVars[k], (q[k][0])) # save the variables that caused this event; (zcf depends on these)
        push!(savedTimes[k], simt)
      end =#
      for k in evDep[modifiedIndex].evCont
        push!(savedVars[k], q[k][0]) # save the variables that are affected by this event before the event; (var=...)
        push!(savedTimes[k], simt)
      end 
      f(-1, -1, modifiedIndex, q, d, t, taylorOpsCache,clF)# execute event----------------no need to clear cache; events touch vectors directly
      for k in evDep[modifiedIndex].evCont
         push!(savedVars[k], q[k][0]) # save the variables that are affected by this event before the event; (var=...)
        push!(savedTimes[k], simt)
        x[k][0]  =q[k][0] ;  tx[k] = simt# update x with q
        quantum[k] = 2*relQ * abs(x[k].coeffs[1]); quantum[k] = quantum[k] < absQ ? 2*absQ : quantum[k]; quantum[k] = quantum[k] > maxErr ? maxErr : quantum[k]
      end
      computeNextEventTime(Val(O), i, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absZ, relZ)
      for j in (HD[modifiedIndex]) # care about dependency to this event only
        elapsedx = simt - tx[j] ; if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx); tx[j] = simt end 
        integrateState(Val(O - 1),q[j],tq[j],simt)  ;tq[j] = simt #q needs to be updated here for recomputeNext                 
        computeDerivatives(Val(O),Val(CS),f,j,q,tq,simt,d,t,taylorOpsCache,clF,x,jac) 
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum) 
      end
      for j in (HZ[modifiedIndex])
         updateZCF(Val(O),Val(CS),f,j,q,tq,simt,d,t,taylorOpsCache,clF,zc_SimpleJac)
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absZ,relZ)
      end
    end
    if stepType != :ST_EVENT
      push!(savedVars[i], x[i][0])
      push!(savedTimes[i], simt)
    else
    #=   for j in (HD[modifiedIndex])
        push!(savedVars[j], x[j][0])
        push!(savedTimes[j], simt)
      end =#
    end
  end
  for i = 1:T
      elapsedx = ft- tx[i]  ;  
      computeDerivatives(Val(O),Val(CS),f,i,q,tq,simt,d,t,taylorOpsCache,clF,x,jac) 
      x[i][0] = x[i](elapsedx);
      push!(savedVars[i], x[i][0]); push!(savedTimes[i], ft) # final point
    end
  stats=Stats(totalSteps,simulStepCount,inpuStepCount,evCount,rejectedEvCount,numStateSteps,numInputSteps)
  createSol(Val(T), Val(O), savedTimes, savedVars, toString(alg), string(odep.prname), absQ, stats, ft)
end