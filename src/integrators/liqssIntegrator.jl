"""
    integrate(Al::QSSAlgorithm{:liqss,O}, CommonqssData::CommonQSS_Data{0}, liqssdata::LiQSS_Data{O,false}, odep::NLODEProblem{PRTYPE,T,D,0,CS}, f::Function, jac::Function, SD::Function, exactA::Function) where {PRTYPE,CS,O,T,D}

Integrates a nonlinear ordinary differential equation (ODE) problem (without events) using the LiQSS (Linearized Quantized State System) algorithm.

# Arguments
- `Al::QSSAlgorithm{:liqss,O}`: The QSS algorithm type with the `liqss` method.
- `CommonqssData::CommonQSS_Data{0}`: Common QSS data structure.
- `liqssdata::LiQSS_Data{O,false}`: LiQSS-specific data structure.
- `odep::NLODEProblem{PRTYPE,T,0,0,CS}`: Nonlinear ODE problem to be solved.
- `f::Function`: The function defining the ODE system.
- `jac::Function`: The Jacobian dependency function of the ODE system.
- `SD::Function`: The state derivative dependency function.
- `exactA::Function`: The exact jacobian expression function for the ODE system.

# Returns
- A solution.

"""
function integrate(Al::QSSAlgorithm{:liqss,O}, CommonqssData::CommonQSS_Data{0}, liqssdata::LiQSS_Data{O,false}, odep::NLODEProblem{PRTYPE,T,D,0,CS}, f::Function, jac::Function, SD::Function, exactA::Function) where {PRTYPE,CS,O,T,D}
  if VERBOSE println("integration...") end
  cacheA = liqssdata.cacheA
  ft = CommonqssData.finalTime
  initTime = CommonqssData.initialTime
  relQ = CommonqssData.relQ
  absQ = CommonqssData.absQ
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
  taylorOpsCache = CommonqssData.taylorOpsCache#Val(CS)=odep.Val(CS)
  d = CommonqssData.d 
  #***************************************************************  
  qaux = liqssdata.qaux
  dxaux = liqssdata.dxaux#olddxSpec=liqssdata.olddxSpec;olddx=liqssdata.olddx
  numSteps = Vector{Int}(undef, T)
  d = [0.0]# this is a dummy var used in updateQ and simulUpdate because in the discrete world exactA needs d
  exactA(q, d, cacheA, 1, 1, initTime + 1e-9)
  #######################################compute initial values##################################################
  n = 1
  for k = 1:O # compute initial derivatives for x and q (similar to a recursive way )
   n = n * k
   for i = 1:T
     q[i].coeffs[k] = x[i].coeffs[k]
   end # q computed from x and it is going to be used in the next x
   for i = 1:T
     clearCache(taylorOpsCache, Val(CS), Val(O))
     f(i, q, t, d,taylorOpsCache)
       if k==1
         x[i].coeffs[k+1] = (taylorOpsCache[1].coeffs[1]) 
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
        f(i, q, t + smallAdvance, d,taylorOpsCache)
        computeNextInputTime(Val(O), i, initTime, smallAdvance, taylorOpsCache[1], nextInputTime, x, quantum)
        if nextInputTime[i] == Inf
          nextInputTime[i] = initTime + 10 * smallAdvance
        end
      end
    else
      updateQ(Val(O), i, x, q, quantum, exactA, d, cacheA, dxaux, qaux, tx, tq, initTime + 1e-12, ft, nextStateTime) #1e-9 exactAfunc contains 1/t???
    end
  end
  ###################################################################################################################################################################
  ####################################################################################################################################################################
  #---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
  ###################################################################################################################################################################
  #################################################################################################################################################################### 
  simt = initTime
  totalSteps = 0
  while simt < ft && totalSteps < maxiters
    if totalSteps == maxiters - 1
      @warn("The algorithm liqss$O is taking too long to converge. The simulation will be stopped. Consider using a different algorithm!")
    end
    sch = updateScheduler(Val(T), nextStateTime, nextEventTime, nextInputTime)
    simt = sch[2]
    index = sch[1]
    if simt > ft
      break
    end
    numSteps[index] += 1
    totalSteps += 1
    t[0] = simt
    ##########################################state########################################
    if sch[3] == :ST_STATE
      elapsed = simt - tx[index]
      integrateState(Val(O), x[index], elapsed)
      tx[index] = simt
      quantum[index] = relQ * abs(x[index].coeffs[1])
      quantum[index] = quantum[index] < absQ ? absQ : quantum[index]
      quantum[index] = quantum[index] > maxErr ? maxErr : quantum[index]
      updateQ(Val(O), index, x, q, quantum, exactA, d, cacheA, dxaux, qaux, tx, tq, simt, ft, nextStateTime)
      tq[index] = simt
      for j in SD(index)
        elapsedx = simt - tx[j]
        if elapsedx > 0
          x[j].coeffs[1] = x[j](elapsedx)
          tx[j] = simt
          elapsedq = simt - tq[j]
          if elapsedq > 0
            integrateState(Val(O - 1), q[j], elapsedq)
            tq[j] = simt
          end
        end
        for b in jac(j)
          elapsedq = simt - tq[b]
          if elapsedq > 0
            integrateState(Val(O - 1), q[b], elapsedq)
            tq[b] = simt
          end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(j, q, t, d,taylorOpsCache)
        computeDerivative(Val(O), x[j], taylorOpsCache[1])
        Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end#end for SD
    ##################################input########################################
    elseif sch[3] == :ST_INPUT  # time of change has come to a state var that does not depend on anything...   
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
      f(index, q, t, d,taylorOpsCache)
      computeDerivative(Val(O), x[index], taylorOpsCache[1])
      computeNextTime(Val(O), index, simt, nextInputTime, x, quantum) #
      if nextInputTime[index] > simt + 2 * elapsed
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(index, q, t + smallAdvance, d,taylorOpsCache)
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
        end#q needs to be updated here for recomputeNext                 
        for b in jac(j)
          elapsedq = simt - tq[b]
          if elapsedq > 0
            integrateState(Val(O - 1), q[b], elapsedq)
            tq[b] = simt
          end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(j, q, t, d,taylorOpsCache)
        computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end#end for
    end#end state/input/event
    push!(savedVars[index], x[index][0])
    push!(savedTimes[index], simt)
  end#end while
  stats=Stats(totalSteps,0,0,numSteps)
  createSol(Val(T), Val(O), savedTimes, savedVars, "liqss$O", string(odep.prname), absQ, stats, ft)
end#end integrate




