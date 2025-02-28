"""
    integrate(Al::QSSAlgorithm{:liqss,O}, CommonqssData::CommonQSS_Data{Z}, liqssdata::LiQSS_Data{O,3}, odep::NLODEProblem{F,PRTYPE,T,D,Z,CS}, f::Function, jac::Function, SD::Function, exactA::Function) where {F,PRTYPE,O,T,D,Z,CS}

Integrates a nonlinear ordinary differential equation (ODE) problem with `events` using the LiQSS (Linearized Quantized State System) algorithm.

# Arguments
- `Al::QSSAlgorithm{:liqss,O}`: The QSS algorithm type for LiQSS.
- `CommonqssData::CommonQSS_Data{Z}`: Common data structure for QSS algorithms.
- `liqssdata::LiQSS_Data{O,1}`: Data structure specific to the LiQSS algorithm.
- `odep::NLODEProblem{F,PRTYPE,T,D,Z,CS}`: The nonlinear ODE problem to be solved.
- `f::Function`: The function defining the ODE system.
- `jac::Function`: The Jacobian dependency function of the ODE system.
- `SD::Function`: The state derivative dependency function.
- `exactA::Function`: The exact jacobian expression function for the ODE system.

# Returns
- A solution after the integration process.

"""
function integrate(Al::QSSAlgorithm{:liqss,O}, CommonqssData::CommonQSS_Data{Z}, liqssdata::LiQSS_Data{O,3}, odep::NLODEProblem{F,PRTYPE,T,D,Z,CS}, f::Function, jac::Function, SD::Function, exactA::Function) where {F,PRTYPE,O,T,D,Z,CS}
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
  qaux = liqssdata.qaux
  dxaux = liqssdata.dxaux
  savedVars = CommonqssData.savedVars
  savedTimes = CommonqssData.savedTimes
  taylorOpsCache = CommonqssData.taylorOpsCache
  #*********************************problem info*****************************************
  d = CommonqssData.d
  clF=odep.closureFuncs[1]
  zc_SimpleJac = odep.ZCjac
  HZ = odep.HZ
  HD = odep.HD
  SZ = odep.SZ
  evDep = odep.eventDependencies
 
  exactA(q, d, cacheA, 1, 1, initTime + 1e-9,clF)
  f(1, -1, -1, q, d, t, taylorOpsCache,clF)
  trackSimul = Vector{Int}(undef, 1)
  #********************************helper values*******************************  
  oldsignValue = MMatrix{Z,2}(zeros(Z * 2))  #usedto track if zc changed sign; each zc has a value and a sign 
  numSteps = Vector{Int}(undef, T)
  #######################################compute initial values##################################################
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


  smallAdvance = ft / 1000
  t[0] = initTime
  for i = 1:T
    numSteps[i] = 0
    push!(savedVars[i], x[i][0])
    push!(savedTimes[i], initTime)
    quantum[i] = relQ * abs(x[i].coeffs[1])
    quantum[i] = quantum[i] < absQ ? absQ : quantum[i]
    quantum[i] = quantum[i] > maxErr ? maxErr : quantum[i]
    if isempty(jac(i))
      computeNextTime(Val(O), i, initTime, nextInputTime, x, quantum)
      if nextInputTime[i] == Inf
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(i, -1, -1, q, d, t + smallAdvance, taylorOpsCache,clF)
        computeNextInputTime(Val(O), i, initTime, smallAdvance, taylorOpsCache[1], nextInputTime, x, quantum)
        if nextInputTime[i] > initTime + 2 * smallAdvance
          nextInputTime[i] = initTime + 2 * smallAdvance
        end
      end
    else
      updateQInit(Val(O), i, x, q, quantum, exactA, d, cacheA, dxaux, qaux, tx, tq, initTime + 1e-12, ft, nextStateTime,clF) 
    end
  end
  for i = 1:Z
    clearCache(taylorOpsCache, Val(CS), Val(O))
    f(-1, i, -1, x, d, t, taylorOpsCache,clF)
    oldsignValue[i, 2] = taylorOpsCache[1][0] #value
    oldsignValue[i, 1] = sign(taylorOpsCache[1][0]) #sign modify 
    computeNextEventTime(Val(O), i, taylorOpsCache[1], oldsignValue, initTime, nextEventTime, quantum, absQ)
  end
  ###################################################################################################################################################################
  ####################################################################################################################################################################
  #---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
  ###################################################################################################################################################################
  ####################################################################################################################################################################
  simt = initTime
  totalSteps = 0
  modifiedIndex = 0
  evCount = 0
  if VERBOSE println("start integration") end
  while simt < ft && totalSteps < maxiters
    if totalSteps == maxiters - 1 @warn("The algorithm liqss$O is taking too long to converge. The simulation will be stopped. Consider using a different algorithm!") end
    sch = updateScheduler(Val(T), nextStateTime, nextEventTime, nextInputTime)
    simt = sch[2]
    index = sch[1]
    stepType = sch[3]
    if simt > ft
      break   ###################################################break##########################################
    end
    totalSteps += 1
    t[0] = simt
    ##########################################state######################################## 
    if stepType == :ST_STATE
      xitemp = x[index][0]
      numSteps[index] += 1
      elapsed = simt - tx[index]
      integrateState(Val(O), x[index], elapsed)
      tx[index] = simt
      dirI = x[index][0] - xitemp
      quantum[index] = relQ * abs(x[index].coeffs[1])
      quantum[index] = quantum[index] < absQ ? absQ : quantum[index]
      quantum[index] = quantum[index] > maxErr ? maxErr : quantum[index]
      firstguess = updateQ(Val(O), index, x, q, quantum, exactA, d, cacheA, dxaux, qaux, tx, tq, simt, ft, nextStateTime,clF)
      tq[index] = simt
      for c in SD(index)   #index influences c  
        elapsedx = simt - tx[c]
        if elapsedx > 0
          x[c].coeffs[1] = x[c](elapsedx)
          tx[c] = simt
        end 
        elapsedq = simt - tq[c]
        if elapsedq > 0
          integrateState(Val(O - 1), q[c], elapsedq)
          tq[c] = simt
        end   # c never been visited 
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(c, -1, -1, q, d, t, taylorOpsCache,clF)
        computeDerivative(Val(O), x[c], taylorOpsCache[1])
        Liqss_reComputeNextTime(Val(O), c, simt, nextStateTime, x, q, quantum)
      end#end for SD
      for j in (SZ[index])
        for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
          elapsedq = simt - tq[b]
          if elapsedq > 0
            integrateState(Val(O - 1), q[b], elapsedq)
            tq[b] = simt
          end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(-1, j, -1, q, d, t, taylorOpsCache,clF)   # run ZCF--------      
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
      end#end for SZ
    ##################################input########################################
    elseif stepType == :ST_INPUT  # time of change has come to a state var that does not depend on anything...
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
      computeNextTime(Val(O), index, simt, nextInputTime, x, quantum) #
      if nextInputTime[index] > simt + 2 * elapsed
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(index, -1, -1, q, d, t + smallAdvance, taylorOpsCache,clF)
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
        # elapsed update all other vars that this derj depends upon.
        for b in jac(j)
          elapsedq = simt - tq[b]
          if elapsedq > 0 integrateState(Val(O - 1), q[b], elapsedq); tq[b] = simt end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(j, -1, -1, q, d, t, taylorOpsCache,clF)
        computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end#end for
      for j in (SD(index))
        elapsedx = simt - tx[j]
        if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx) ;tx[j] = simt end
        elapsedq = simt - tq[j]
        if elapsedq > 0 integrateState(Val(O - 1), q[j], elapsedq) ;tq[j] = simt end#q needs to be updated here for recomputeNext                 
        # elapsed update all other vars that this derj depends upon.
        for b in jac(j)
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
      end#end for
      for j in (SZ[index])
        for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
          elapsedq = simt - tq[b]
          if elapsedq > 0
            integrateState(Val(O - 1), q[b], elapsedq)
            tq[b] = simt
          end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(-1, j, -1, q, d, t, taylorOpsCache,clF)   # run ZCF--------      
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
      end#end for SZ
    #################################################################event########################################

    else
      for b in zc_SimpleJac[index] # elapsed update all other vars that this zc depends upon.
        elapsedq = simt - tq[b]
        if elapsedq > 0
          integrateState(Val(O - 1), q[b], elapsedq)
          tq[b] = simt
        end
      end
      #first we have a zc happened which corresponds to nexteventtime and index (one of zc) but we want also the sign to know ev+ or ev- 
      clearCache(taylorOpsCache, Val(CS), Val(O))
      f(-1, index, -1, q, d, t, taylorOpsCache,clF)    # run ZCF-------- 
      if oldsignValue[index, 2] * taylorOpsCache[1][0] >= 0
        if abs(taylorOpsCache[1][0]) > 1e-9 * absQ # if both have same sign and zcf is not very small: zc=1e-9*absQ is allowed as an event
          computeNextEventTime(Val(O), index, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
          continue
        end
      end
      if abs(oldsignValue[index, 2]) <= 1e-9 * absQ  #earlier zc=1e-9*absQ was considered event , so now it should be prevented from passing
        nextEventTime[index] = Inf # at this instant next zc will be triggered now, and this will lead to infinite events, so cannot computenextevent here
        continue
      end
      if taylorOpsCache[1][0] > oldsignValue[index, 2] #scheduled rise
        modifiedIndex = 2 * index - 1
      elseif taylorOpsCache[1][0] < oldsignValue[index, 2] #scheduled drop
        modifiedIndex = 2 * index
      else # == ( zcf==oldZCF)
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
      f(-1, -1, modifiedIndex, x, d, t, taylorOpsCache,clF)# execute event----------------no need to clear cache; events touch vectors directly
      for i in evDep[modifiedIndex].evCont
        quantum[i] = relQ * abs(x[i].coeffs[1])
        quantum[i] = quantum[i] < absQ ? absQ : quantum[i]
        quantum[i] = quantum[i] > maxErr ? maxErr : quantum[i]
        firstguess = updateQ(Val(O), i, x, q, quantum, exactA, d, cacheA, dxaux, qaux, tx, tq, simt, ft, nextStateTime,clF)
        tx[i] = simt;tq[i] = simt
        Liqss_reComputeNextTime(Val(O), i, simt, nextStateTime, x, q, quantum)
      end
      computeNextEventTime(Val(O), index, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ) #update zcf before thiscatch in qss quantizer to avoid infinite events
      for j in (HD[modifiedIndex]) # care about dependency to this event only
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
        for b = 1:T # elapsed update all other vars that this derj depends upon.
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
        Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end
      for j in (HZ[modifiedIndex])
        for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
          elapsedq = simt - tq[b]
          if elapsedq > 0
            integrateState(Val(O - 1), q[b], elapsedq)
            tq[b] = simt
          end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(-1, j, -1, q, d, t, taylorOpsCache,clF)  # run ZCF-------- 
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
      end
    end#end state/input/event
    if stepType != :ST_EVENT
      push!(savedVars[index], x[index][0])
      push!(savedTimes[index], simt)
    else
      for j in (HD[modifiedIndex])
        push!(savedVars[j], x[j][0])
        push!(savedTimes[j], simt)
      end
    end
  end#end while
  stats=Stats(totalSteps,0,evCount,numSteps)
  createSol(Val(T), Val(O), savedTimes, savedVars, "liqss$O", string(odep.prname), absQ, stats, ft) #= ,savedDers =#
end#end integrate