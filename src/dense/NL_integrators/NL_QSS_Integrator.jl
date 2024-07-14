function integrate(Al::QSSAlgorithm{:qss,O}, CommonqssData::CommonQSS_data{0}, odep::NLODEProblem{PRTYPE,T,0,0,CS}, f::Function, jac::Function, SD::Function) where {PRTYPE,O,T,CS}
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
  numSteps = Vector{Int}(undef, T)
  n = 1
  for k = 1:O
    n = n * k
    for i = 1:T
      q[i].coeffs[k] = x[i].coeffs[k]
    end
    for i = 1:T
      clearCache(taylorOpsCache, Val(CS), Val(O))
      f(i, q, t, taylorOpsCache)
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
        f(i, q, t + smallAdvance, taylorOpsCache)
        computeNextInputTime(Val(O), i, initTime, smallAdvance, taylorOpsCache[1], nextInputTime, x, quantum)
        if nextInputTime[i] > initTime + 2 * smallAdvance
          nextInputTime[i] = initTime + 2 * smallAdvance
        end
      end
    else
      computeNextTime(Val(O), i, initTime, nextStateTime, x, quantum)
    end
  end
  simt = initTime
  totalSteps = 0
  prevStepTime = initTime
  while simt < ft && totalSteps < maxiters
    if totalSteps == maxiters - 1 @warn("The algorithm qss$O is taking too long to converge. The simulation will be stopped. Consider using a different algorithm!") end
    sch = updateScheduler(Val(T), nextStateTime, nextEventTime, nextInputTime)
    simt = sch[2]
    index = sch[1]
    totalSteps += 1
    t[0] = simt
    if sch[3] == :ST_STATE
      numSteps[index] += 1
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
        f(j, q, t, taylorOpsCache)
        computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end
    elseif sch[3] == :ST_INPUT
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
     #=  for b in jac(index)
        elapsedq = simt - tq[b]
        if elapsedq > 0 ;integrateState(Val(O - 1), q[b], elapsedq) ;tq[b] = simt end
      end =#
      clearCache(taylorOpsCache, Val(CS), Val(O))
      f(index, q, t, taylorOpsCache)
      computeDerivative(Val(O), x[index], taylorOpsCache[1])
      computeNextTime(Val(O), index, simt, nextInputTime, x, quantum)
      if nextInputTime[index] > simt + 2 * elapsed
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(index, q, t + smallAdvance, taylorOpsCache)
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
        f(j, q, t, taylorOpsCache)
        computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end
    end
    push!(savedVars[index], x[index][0])
    push!(savedTimes[index], simt)
  end
  stats=Stats(totalSteps,0,0,numSteps)
  createSol(Val(T), Val(O), savedTimes, savedVars, "qss$O", string(odep.prname), absQ, stats, ft)
end
