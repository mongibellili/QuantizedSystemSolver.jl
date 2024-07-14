function integrate(Al::QSSAlgorithm{:qss,O}, CommonqssData::CommonQSS_data{Z}, odep::NLODEProblem{PRTYPE,T,Z,Y,CS}, f::Function, jac::Function, SD::Function) where {PRTYPE,O,T,Z,Y,CS}
  ft = CommonqssData.finalTime
  initTime = CommonqssData.initialTime
  relQ = CommonqssData.dQrel
  absQ = CommonqssData.dQmin
  maxErr = CommonqssData.maxErr
  maxiters = CommonqssData.maxiters
  quantum = CommonqssData.quantum
  nextStateTime = CommonqssData.nextStateTime
  nextEventTime = CommonqssData.nextEventTime
  nextInputTime = CommonqssData.nextInputTime
  tx = CommonqssData.tx
  tq = CommonqssData.tq
  x = CommonqssData.x
  q = CommonqssData.q
  t = CommonqssData.t
  savedVars = CommonqssData.savedVars
  savedTimes = CommonqssData.savedTimes
  integratorCache = CommonqssData.integratorCache
  taylorOpsCache = CommonqssData.taylorOpsCache
  d = CommonqssData.d
  zc_SimpleJac = odep.ZCjac
  HZ = odep.HZ
  HD = odep.HD
  SZ = odep.SZ
  evDep = odep.eventDependencies
  oldsignValue = MMatrix{Z,2}(zeros(Z * 2))
  numSteps = Vector{Int}(undef, T)
  n = 1
  for k = 1:O
    n = n * k
    for i = 1:T
      q[i].coeffs[k] = x[i].coeffs[k]
    end
    for i = 1:T
      clearCache(taylorOpsCache, Val(CS), Val(O))
      f(i, -1, -1, q, d, t, taylorOpsCache)
      ndifferentiate!(integratorCache, taylorOpsCache[1], k - 1)
      x[i].coeffs[k+1] = (integratorCache.coeffs[1]) / n
    end
  end
  smallAdvance = ft / 100
  t[0] = initTime
  for i = 1:T
    numSteps[i] = 0
    push!(savedVars[i], x[i][0])
    push!(savedTimes[i], 0.0)
    quantum[i] = relQ * abs(x[i].coeffs[1])
    quantum[i] = quantum[i] < absQ ? absQ : quantum[i]
    quantum[i] = quantum[i] > maxErr ? maxErr : quantum[i]
    if isempty(jac(i))
      computeNextTime(Val(O), i, initTime, nextInputTime, x, quantum)
      if nextInputTime[i] == Inf
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(i, -1, -1, q, d, t + smallAdvance, taylorOpsCache)
        computeNextInputTime(Val(O), i, initTime, smallAdvance, taylorOpsCache[1], nextInputTime, x, quantum)
        if nextInputTime[i] > initTime + 2 * smallAdvance
          nextInputTime[i] = initTime + 2 * smallAdvance
        end
      end
    else
      computeNextTime(Val(O), i, initTime, nextStateTime, x, quantum)
    end
  end
  for i = 1:Z
    clearCache(taylorOpsCache, Val(CS), Val(O))
    f(-1, i, -1, x, d, t, taylorOpsCache)
    oldsignValue[i, 2] = taylorOpsCache[1][0]
    oldsignValue[i, 1] = sign(taylorOpsCache[1][0])
    computeNextEventTime(Val(O), i, taylorOpsCache[1], oldsignValue, initTime, nextEventTime, quantum, absQ)
  end
  simt = initTime
  totalSteps = 0
  prevStepTime = initTime
  modifiedIndex = 0
  statestep = 0
  evCount = 0
  while simt < ft && totalSteps < maxiters
    if totalSteps == maxiters - 1 @warn("The algorithm qss$O is taking too long to converge. The simulation will be stopped. Consider using a different algorithm!") end
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
      numSteps[index] += 1
      elapsed = simt - tx[index]
      integrateState(Val(O), x[index], elapsed)
      tx[index] = simt
      quantum[index] = relQ * abs(x[index].coeffs[1])
      quantum[index] = quantum[index] < absQ ? absQ : quantum[index]
      quantum[index] = quantum[index] > maxErr ? maxErr : quantum[index]
     #=  if abs(x[index].coeffs[2]) > 1e7
        quantum[index] = 10 * quantum[index]
      end =#
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
        f(j, -1, -1, q, d, t, taylorOpsCache)
        computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end
      for j in (SZ[index])
        for b in zc_SimpleJac[j]
          elapsedq = simt - tq[b]
          if elapsedq > 0 integrateState(Val(O - 1), q[b], elapsedq) ;tq[b] = simt end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(-1, j, -1, q, d, t, taylorOpsCache)
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
      end
    elseif stepType == :ST_INPUT
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
      #= for b in jac(index)
        elapsedq = simt - tq[b]
        if elapsedq > 0 ;integrateState(Val(O - 1), q[b], elapsedq) ;tq[b] = simt end
      end =#
      clearCache(taylorOpsCache, Val(CS), Val(O))
      f(index, -1, -1, q, d, t, taylorOpsCache)
      computeDerivative(Val(O), x[index], taylorOpsCache[1])
      computeNextTime(Val(O), index, simt, nextInputTime, x, quantum)
      if nextInputTime[index] > simt + 2 * elapsed
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(index, -1, -1, q, d, t + smallAdvance, taylorOpsCache)
        computeNextInputTime(Val(O), index, simt, smallAdvance, taylorOpsCache[1], nextInputTime, x, quantum)
        if nextInputTime[index] > simt + 2 * elapsed
          nextInputTime[index] = simt + 2 * elapsed
        end
      end
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
        for b in jac(j)
          elapsedq = simt - tq[b]
          if elapsedq > 0 integrateState(Val(O - 1), q[b], elapsedq) ;tq[b] = simt end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(j, -1, -1, q, d, t, taylorOpsCache)
        computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end
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
        f(j, -1, -1, q, d, t, taylorOpsCache)
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
        f(-1, j, -1, q, d, t, taylorOpsCache)
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
      f(-1, index, -1, q, d, t, taylorOpsCache)
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
      f(-1, -1, modifiedIndex, x, d, t, taylorOpsCache)
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
        f(j, -1, -1, q, d, t, taylorOpsCache)
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
        f(-1, j, -1, q, d, t, taylorOpsCache)
 
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
    #prevStepTime = simt
  end
  stats=Stats(totalSteps,0,evCount,numSteps)
  createSol(Val(T), Val(O), savedTimes, savedVars, "qss$O", string(odep.prname), absQ, stats, ft)
end