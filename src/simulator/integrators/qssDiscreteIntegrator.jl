"""
    integrate(alg::QSSAlgorithm{:qss,O}, commonQssData::CommonQSS_Data{Z}, odep::ODEProblemData{F,PRTYPE,T,D,Z,CS}, f::Function, jac::Function, SD::Function) where {F,PRTYPE,O,T,D,Z,CS}

Integrates a nonlinear ordinary differential equation (ODE) problem with `events` using a Quantized State System (QSS) algorithm.

# Arguments
- `Al::QSSAlgorithm{:qss,O}`: The QSS algorithm to be used for integration.
- `commonQssData::CommonQSS_Data{Z}`: Common data structure for QSS algorithms.
- `odep::ODEProblemData{F,PRTYPE,T,D,Z,CS}`: The nonlinear ODE problem to be solved.
- `f::Function`: The function defining the ODE system.
- `jac::Function`: The Jacobian dependency function of the ODE system.
- `SD::Function`: The state derivative dependency function.


# Type Parameters
- `PRTYPE`: The type of the problem.
- `O`: The order of the QSS algorithm.
- `T`: The number of continuous variables.
- `Z`: The number of zero crossing functions.
- `D`: The number of discrete variables
- `CS`: The cache size.

# Returns
- A solution
"""
function integrate(alg::QSSAlgorithm{:qss,O}, commonQssData::CommonQSS_Data{Z}, odep::ODEProblemData{F,PRTYPE,T,D,Z,CS}, f::Function, jac::Function, SD::Function) where {F,PRTYPE,O,T,D,Z,CS}
  VERBOSE=commonQssData.verbose
  if VERBOSE println("integration...") end
  ft = commonQssData.finalTime
  initTime = commonQssData.initialTime
  relQ = commonQssData.relQ
  absQ = commonQssData.absQ
  maxErr = commonQssData.maxErr
  maxiters = commonQssData.maxiters
  quantum = commonQssData.quantum
  nextStateTime = commonQssData.nextStateTime
  nextEventTime = commonQssData.nextEventTime
  nextInputTime = commonQssData.nextInputTime
  tx = commonQssData.tx
  tq = commonQssData.tq
  x = commonQssData.x
  q = commonQssData.q
  t = commonQssData.t
  savedVars = commonQssData.savedVars
  savedTimes = commonQssData.savedTimes
  taylorOpsCache = commonQssData.taylorOpsCache
  d = commonQssData.d
  clF=odep.closureFuncs[1]
  zc_SimpleJac = odep.ZCjac
  HZ = odep.HZ
  HD = odep.HD
  SZ = odep.SZ
  evDep = odep.eventDependencies
  oldsignValue = MMatrix{Z,2}(zeros(Z * 2))
  numStateSteps = Vector{Int}(undef, T)
  numInputSteps = Vector{Int}(undef, T)
  n = 1
  for k = 1:O # compute initial derivatives for x and q (similar to a recursive way )
   n = n * k
   for i = 1:T
     q[i].coeffs[k] = x[i].coeffs[k]
   end # q computed from x and it is going to be used in the next x
   for i = 1:T
     clearCache(taylorOpsCache, Val(CS), Val(O))
      f(i, -1, -1, q, d, t, taylorOpsCache,clF)
       if k==1
         x[i].coeffs[k+1] = (taylorOpsCache[1].coeffs[1])  # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
       elseif k==2
         x[i].coeffs[k+1] = (taylorOpsCache[1].coeffs[2])/2
       elseif k==3
         x[i].coeffs[k+1] = (taylorOpsCache[1].coeffs[3])/6
       else
         println("order4 not implemented")
       end
     end
 end
  smallAdvance = ft / 100
  t[0] = initTime
  for i = 1:T
    numStateSteps[i] = 0
    numInputSteps[i] = 0
    push!(savedVars[i], x[i][0])
    push!(savedTimes[i], initTime)
    quantum[i] = relQ * abs(x[i].coeffs[1])
    quantum[i] = quantum[i] < absQ ? absQ : quantum[i]
    quantum[i] = quantum[i] > maxErr ? maxErr : quantum[i]
    if isempty(jac(i))
      #= computeNextTime(Val(O), i, initTime, nextInputTime, x, quantum)
      if nextInputTime[i] == Inf
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(i, -1, -1, q, d, t + smallAdvance, taylorOpsCache,clF) =#
       # discrete_computeNextInputTime(Val(O), i, initTime, smallAdvance, taylorOpsCache[1], nextInputTime, x, quantum)
        discrete_computeNextInputTime(Val(O), i, t, f,clF,d, taylorOpsCache, nextInputTime, x,q, quantum) 
        #= if nextInputTime[i] > initTime + 2 * smallAdvance
          nextInputTime[i] = initTime + 2 * smallAdvance
        end =#
      #end
    else
      computeNextTime(Val(O), i, initTime, nextStateTime, x, quantum)
    end
  end
  for i = 1:Z
    clearCache(taylorOpsCache, Val(CS), Val(O))
    f(-1, i, -1, x, d, t, taylorOpsCache,clF)
    oldsignValue[i, 2] = taylorOpsCache[1][0]
    oldsignValue[i, 1] = sign(taylorOpsCache[1][0])
    computeNextEventTime(Val(O), i, taylorOpsCache[1], oldsignValue, initTime, nextEventTime, quantum, absQ)
  end
  simt = initTime
  totalSteps = 0
  modifiedIndex = 0
  statestep = 0
  evCount = 0
   inpuStepCount=0
  ###################################################################################################################################################################
  ####################################################################################################################################################################
  #---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
  ###################################################################################################################################################################
  #################################################################################################################################################################### 
  
  while simt < ft && totalSteps < maxiters
    if totalSteps == maxiters - 1 @warn("The algorithm qss$O reached max iterations. The simulation will be stopped. Consider using a different algorithm!") end
    sch = updateScheduler(Val(T), nextStateTime, nextEventTime, nextInputTime)
    simt = sch[2]
    index = sch[1]
    stepType = sch[3]
    if simt > ft
      break
    end
    totalSteps += 1
    t[0] = simt
    if stepType == :ST_STATE
      statestep += 1
      numStateSteps[index] += 1
      elapsed = simt - tx[index]
      integrateState(Val(O), x[index], elapsed)
      tx[index] = simt
      quantum[index] = relQ * abs(x[index].coeffs[1])
      quantum[index] = quantum[index] < absQ ? absQ : quantum[index]
      quantum[index] = quantum[index] > maxErr ? maxErr : quantum[index]
      for k = 1:O
        q[index].coeffs[k] = x[index].coeffs[k]
      end
      tq[index] = simt
      computeNextTime(Val(O), index, simt, nextStateTime, x, quantum)
      for j in (SD(index))
        elapsedx = simt - tx[j]
        if elapsedx > 0
          x[j].coeffs[1] = x[j](elapsedx)
          tx[j] = simt
        end
        elapsedq = simt - tq[j]
        if elapsedq > 0
          integrateState(Val(O - 1), q[j], elapsedq)
          tq[j] = simt
        end
        for b in (jac(j))
          elapsedq = simt - tq[b]
          if elapsedq > 0
            integrateState(Val(O - 1), q[b], elapsedq)
            tq[b] = simt
          end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(j, -1, -1, q, d, t, taylorOpsCache,clF)
        computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end
      for j in (SZ[index])
        for b in zc_SimpleJac[j]
          elapsedq = simt - tq[b]
          if elapsedq > 0 integrateState(Val(O - 1), q[b], elapsedq) ;tq[b] = simt end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(-1, j, -1, q, d, t, taylorOpsCache,clF)
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
      end
    elseif stepType == :ST_INPUT
       inpuStepCount+=1
      numInputSteps[index] += 1
      elapsed = simt - tx[index]
      integrateState(Val(O), x[index], elapsed)
      tx[index] = simt
      quantum[index] = relQ * abs(x[index].coeffs[1])
      quantum[index] = quantum[index] < absQ ? absQ : quantum[index]
      quantum[index] = quantum[index] > maxErr ? maxErr : quantum[index]
      for k = 1:O
        q[index].coeffs[k] = x[index].coeffs[k]
      end
      tq[index] = simt
      clearCache(taylorOpsCache, Val(CS), Val(O))
      f(index, -1, -1, q, d, t, taylorOpsCache,clF)
      computeDerivative(Val(O), x[index], taylorOpsCache[1])
      discrete_computeNextInputTime(Val(O), index, t, f,clF,d, taylorOpsCache, nextInputTime, x,q, quantum)
     #=  computeNextTime(Val(O), index, simt, nextInputTime, x, quantum)
      if nextInputTime[index] > simt + 2 * elapsed
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(index, -1, -1, q, d, t + smallAdvance, taylorOpsCache,clF)
        discrete_computeNextInputTime(Val(O), index, simt, smallAdvance, taylorOpsCache[1], nextInputTime, x, quantum)
        if nextInputTime[index] > simt + 2 * elapsed
          nextInputTime[index] = simt + 2 * elapsed
        end
      end =#
  
      for j in (SD(index))
        elapsedx = simt - tx[j]
        if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx) ;tx[j] = simt  end
        elapsedq = simt - tq[j]
        if elapsedq > 0 integrateState(Val(O - 1), q[j], elapsedq) ;tq[j] = simt end
        for b in jac(j)
          elapsedq = simt - tq[b]
          if elapsedq > 0 integrateState(Val(O - 1), q[b], elapsedq) ;tq[b] = simt end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(j, -1, -1, q, d, t, taylorOpsCache,clF)
        computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end
      for j in (SZ[index])
        for b in zc_SimpleJac[j]
          elapsedq = simt - tq[b]
          if elapsedq > 0
            integrateState(Val(O - 1), q[b], elapsedq)
            tq[b] = simt
          end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(-1, j, -1, q, d, t, taylorOpsCache,clF)
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
      end
    else
      for b in zc_SimpleJac[index]
        elapsedq = simt - tq[b]
        if elapsedq > 0
          integrateState(Val(O - 1), q[b], elapsedq)
          tq[b] = simt
        end
      end
      clearCache(taylorOpsCache, Val(CS), Val(O))
      f(-1, index, -1, q, d, t, taylorOpsCache,clF)
      if oldsignValue[index, 2] * taylorOpsCache[1][0] >= 0
        if abs(taylorOpsCache[1][0]) > 1e-9 * absQ
          computeNextEventTime(Val(O), index, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
          continue
        end
      end
      if abs(oldsignValue[index, 2]) <= 1e-9 * absQ
        nextEventTime[index] = Inf
        continue
      end
      if taylorOpsCache[1][0] > oldsignValue[index, 2]
        modifiedIndex = 2 * index - 1
      elseif taylorOpsCache[1][0] < oldsignValue[index, 2]
        modifiedIndex = 2 * index
      else
        computeNextEventTime(Val(O), index, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
        continue
      end
      evCount += 1
      oldsignValue[index, 2] = taylorOpsCache[1][0]
      oldsignValue[index, 1] = sign(taylorOpsCache[1][0])
      for b in evDep[modifiedIndex].evContRHS
        elapsedq = simt - tq[b]
        if elapsedq > 0
          integrateState(Val(O - 1), q[b], elapsedq)
          tq[b] = simt
        end
      end
      f(-1, -1, modifiedIndex, x, d, t, taylorOpsCache,clF)
      for i in evDep[modifiedIndex].evCont
        quantum[i] = relQ * abs(x[i].coeffs[1])
        quantum[i] = quantum[i] < absQ ? absQ : quantum[i]
        quantum[i] = quantum[i] > maxErr ? maxErr : quantum[i]
        for k = 0:O-1 q[i][k] = x[i][k]  end
        tx[i] = simt
        tq[i] = simt
        computeNextTime(Val(O), i, simt, nextStateTime, x, quantum)
      end
      computeNextEventTime(Val(O), index, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
      for j in (HD[modifiedIndex])
        elapsedx = simt - tx[j]
        if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx) ;tx[j] = simt  end
        elapsedq = simt - tq[j]
        if elapsedq > 0
          integrateState(Val(O - 1), q[j], elapsedq)
          tq[j] = simt
        end
        for b = 1:T
          if b in jac(j)
            elapsedq = simt - tq[b]
            if elapsedq > 0
              integrateState(Val(O - 1), q[b], elapsedq)
              tq[b] = simt
            end
          end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(j, -1, -1, q, d, t, taylorOpsCache,clF)
        computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end
      for j in (HZ[modifiedIndex])
        for b in zc_SimpleJac[j]
          elapsedq = simt - tq[b]
          if elapsedq > 0
            integrateState(Val(O - 1), q[b], elapsedq)
            tq[b] = simt
          end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(-1, j, -1, q, d, t, taylorOpsCache,clF)
 
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
      end
    end
    if stepType != :ST_EVENT
      push!(savedVars[index], x[index][0])
      push!(savedTimes[index], simt)
    else
      for j in (HD[modifiedIndex])
        push!(savedVars[j], x[j][0])
        push!(savedTimes[j], simt)
      end
    end
  end
  stats=Stats(totalSteps,0,inpuStepCount,evCount,numStateSteps,numInputSteps)
  createSol(Val(T), Val(O), savedTimes, savedVars, toString(alg), string(odep.prname), absQ, stats, ft)
end