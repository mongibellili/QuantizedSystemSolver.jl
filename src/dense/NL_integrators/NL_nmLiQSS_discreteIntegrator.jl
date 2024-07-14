function integrate(Al::QSSAlgorithm{:nmliqss,O}, CommonqssData::CommonQSS_data{Z}, liqssdata::LiQSS_data{O,false}, odep::NLODEProblem{PRTYPE,T,Z,D,CS}, f::Function, jac::Function, SD::Function, exactA::Function) where {PRTYPE,O,T,Z,D,CS}
  if VERBOSE println("begining of intgrate function") end
  cacheA = liqssdata.cacheA
  ft = CommonqssData.finalTime
  initTime = CommonqssData.initialTime
  relQ = CommonqssData.dQrel
  absQ = CommonqssData.dQmin
  maxErr = CommonqssData.maxErr
  maxiters = CommonqssData.maxiters
  savetimeincrement = CommonqssData.savetimeincrement
  savetime = savetimeincrement
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
  taylorOpsCache = CommonqssData.taylorOpsCache#cacheSize=odep.cacheSize

  #*********************************problem info*****************************************
  d = CommonqssData.d
  zc_SimpleJac = odep.ZCjac
  HZ = odep.HZ
  HD = odep.HD
  SZ = odep.SZ
  evDep = odep.eventDependencies
  qaux = liqssdata.qaux
  dxaux = liqssdata.dxaux
  # to be deleted: used for root finding using Newton_Intrvals
  pp = pointer(Vector{NTuple{2,Float64}}(undef, 7))
  respp = pointer(Vector{Float64}(undef, 6))
  acceptedi = Vector{Vector{Float64}}(undef, 4 * O - 1)
  acceptedj = Vector{Vector{Float64}}(undef, 4 * O - 1)
  cacheRootsi = zeros(8 * O - 4) #vect of floats to hold roots for simul_analytic 
  cacheRootsj = zeros(8 * O - 4)
  for i = 1:4*O-1
    acceptedi[i] = [0.0, 0.0]#zeros(2)
    acceptedj[i] = [0.0, 0.0]#zeros(2)
  end
  exactA(q, d, cacheA, 1, 1, initTime + 1e-9)
  f(1, -1, -1, q, d, t, taylorOpsCache)
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
      f(i, -1, -1, q, d, t, taylorOpsCache)
      ndifferentiate!(integratorCache, taylorOpsCache[1], k - 1)
      x[i].coeffs[k+1] = (integratorCache.coeffs[1]) / n # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
    end
  end
  smallAdvance = ft / 1000
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
      updateQ(Val(O), i, x, q, quantum, exactA, d, cacheA, dxaux, qaux, tx, tq, initTime + 1e-12, ft, nextStateTime) #1e-9 exactAfunc contains 1/t???
    end
  end
  for i = 1:Z
    clearCache(taylorOpsCache, Val(CS), Val(O))
    f(-1, i, -1, x, d, t, taylorOpsCache)
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
  prevStepTime = initTime
  modifiedIndex = 0
  countEvents = 0
  inputstep = 0
  statestep = 0
  simulStepCount = 0
  ft < savetime && error("ft<savetime")
  if VERBOSE println("start integration") end
  while simt < ft && totalSteps < maxiters
    if totalSteps == maxiters - 1 @warn("The algorithm nmliqss$O is taking too long to converge. The simulation will be stopped. Consider using a different algorithm!") end
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
      statestep += 1
      xitemp = x[index][0]
      numSteps[index] += 1
      elapsed = simt - tx[index]
      integrateState(Val(O), x[index], elapsed)
      tx[index] = simt
      dirI = x[index][0] - xitemp
      quantum[index] = relQ * abs(x[index].coeffs[1])
      quantum[index] = quantum[index] < absQ ? absQ : quantum[index]
      quantum[index] = quantum[index] > maxErr ? maxErr : quantum[index]
      for b in (jac(index))    # update Qb : to be used to calculate exacte Aindexb
        elapsedq = simt - tq[b]
        if elapsedq > 0
          integrateState(Val(O - 1), q[b], elapsedq)
          tq[b] = simt
        end
      end
      firstguess = updateQ(Val(O), index, x, q, quantum, exactA, d, cacheA, dxaux, qaux, tx, tq, simt, ft, nextStateTime)
      tq[index] = simt
      #----------------------------------------------------check dependecy cycles---------------------------------------------  
      trackSimul[1] = 0
      for j in SD(index)
        for b in (jac(j))    # update Qb: to be used to calculate exacte Ajb
          elapsedq = simt - tq[b]
          if elapsedq > 0
            integrateState(Val(O - 1), q[b], elapsedq)
            tq[b] = simt
          end
        end
        cacheA[1] = 0.0
        exactA(q, d, cacheA, index, j, simt)
        aij = cacheA[1]# can be passed to simul so that i dont call exactfunc again
        cacheA[1] = 0.0
        exactA(q, d, cacheA, j, index, simt)
        aji = cacheA[1]
        if j != index && aij * aji != 0.0
          if nmisCycle_and_simulUpdate(cacheRootsi, cacheRootsj, acceptedi, acceptedj, aij, aji, respp, pp, trackSimul, Val(O), index, j, dirI, firstguess, x, q, quantum, exactA, d, cacheA, dxaux, qaux, tx, tq, simt, ft)
            simulStepCount += 1
            clearCache(taylorOpsCache, Val(CS), Val(O))
            f(index, -1, -1, q, d, t, taylorOpsCache)
            computeDerivative(Val(O), x[index], taylorOpsCache[1])
            for k in SD(j)  #j influences k
              if k != index && k != j
                elapsedx = simt - tx[k]
                x[k].coeffs[1] = x[k](elapsedx)
                tx[k] = simt
                elapsedq = simt - tq[k]
                if elapsedq > 0
                  integrateState(Val(O - 1), q[k], elapsedq)
                  tq[k] = simt
                end
                for b in (jac(k))
                  elapsedq = simt - tq[b]
                  if elapsedq > 0
                    integrateState(Val(O - 1), q[b], elapsedq)
                    tq[b] = simt
                  end
                end
                clearCache(taylorOpsCache, Val(CS), Val(O))
                f(k, -1, -1, q, d, t, taylorOpsCache)
                computeDerivative(Val(O), x[k], taylorOpsCache[1])
                Liqss_reComputeNextTime(Val(O), k, simt, nextStateTime, x, q, quantum)
              end#end if k!=0
            end#end for k depend on j     
            for k in (SZ[j]) # qj changed, so zcf should be checked
              for b in zc_SimpleJac[k] # elapsed update all other vars that this derj depends upon.
                elapsedq = simt - tq[b]
                if elapsedq > 0 integrateState(Val(O - 1), q[b], elapsedq) ;tq[b] = simt end
              end
              clearCache(taylorOpsCache, Val(CS), Val(O))
              f(-1, k, -1, q, d, t, taylorOpsCache)   # run ZCF--------      
              computeNextEventTime(Val(O), k, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
            end#end for SZ
          end#end ifcycle check
        end#end if j!=0
      end#end FOR_cycle check
      if trackSimul[1] != 0  #qi changed after throw
        Liqss_reComputeNextTime(Val(O), index, simt, nextStateTime, x, q, quantum)
      end

      #-------------------------------------------------------------------------------------
      #---------------------------------normal liqss: proceed--------------------------------
      #-------------------------------------------------------------------------------------
      for c in SD(index)   #index influences c  
        elapsedx = simt - tx[c]
        if elapsedx > 0
          x[c].coeffs[1] = x[c](elapsedx)
          tx[c] = simt
        end # 
        elapsedq = simt - tq[c]
        if elapsedq > 0 integrateState(Val(O - 1), q[c], elapsedq) ;tq[c] = simt end   # c never been visited 
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(c, -1, -1, q, d, t, taylorOpsCache)
        computeDerivative(Val(O), x[c], taylorOpsCache[1])
        Liqss_reComputeNextTime(Val(O), c, simt, nextStateTime, x, q, quantum)
      end#end for SD
      for j in (SZ[index])
        for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
          elapsedq = simt - tq[b]
          if elapsedq > 0  integrateState(Val(O - 1), q[b], elapsedq) ;tq[b] = simt end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(-1, j, -1, q, d, t, taylorOpsCache)   # run ZCF--------      
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
      end#end for SZ
    ##################################input########################################
    elseif stepType == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself    
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
      f(index, -1, -1, q, d, t, taylorOpsCache)
      computeDerivative(Val(O), x[index], taylorOpsCache[1])
      computeNextTime(Val(O), index, simt, nextInputTime, x, quantum) #
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
        end#q needs to be updated here for recomputeNext                 
        # elapsed update all other vars that this derj depends upon.
        for b in jac(j)
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
      end#end for
      for j in (SD(index))
        elapsedx = simt - tx[j]
        if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx) ;tx[j] = simt end
        elapsedq = simt - tq[j]
        if elapsedq > 0 integrateState(Val(O - 1), q[j], elapsedq) ;tq[j] = simt end#q needs to be updated here for recomputeNext                 
        # elapsed update all other vars that this derj depends upon.
        for b in jac(j)
          elapsedq = simt - tq[b]
          if elapsedq > 0 integrateState(Val(O - 1), q[b], elapsedq) ;tq[b] = simt end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(j, -1, -1, q, d, t, taylorOpsCache)
        computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end#end for
      for j in (SZ[index])
        for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
          elapsedq = simt - tq[b]
          if elapsedq > 0 integrateState(Val(O - 1), q[b], elapsedq); tq[b] = simt end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(-1, j, -1, q, d, t, taylorOpsCache)   # run ZCF--------      
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
      clearCache(taylorOpsCache, Val(CS), Val(O))
      f(-1, index, -1, q, d, t, taylorOpsCache)    # run ZCF again to verify-------- 
      if oldsignValue[index, 2] * taylorOpsCache[1][0] >= 0 # if computeNextEvent errored 
        if abs(taylorOpsCache[1][0]) > 1e-9 * absQ # if error is negligeable then ok consider as event, else reject....if both have same sign and zcf is not very small: zc==1e-9*absQ is allowed as an event
          computeNextEventTime(Val(O), index, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
          continue #event rejected
        end
      end
      if abs(oldsignValue[index, 2]) <= 1e-9 * absQ  #earlier zc==1e-9*absQ was considered event , so now it should be prevented from passing
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
      countEvents += 1
      oldsignValue[index, 2] = taylorOpsCache[1][0]
      oldsignValue[index, 1] = sign(taylorOpsCache[1][0])
      for b in evDep[modifiedIndex].evContRHS
        elapsedq = simt - tq[b]
        if elapsedq > 0
          integrateState(Val(O - 1), q[b], elapsedq)
          tq[b] = simt
        end
      end
      f(-1, -1, modifiedIndex, x, d, t, taylorOpsCache)# execute event----------------no need to clear cache; events touch vectors directly
      for i in evDep[modifiedIndex].evCont
        quantum[i] = relQ * abs(x[i].coeffs[1])
        quantum[i] = quantum[i] < absQ ? absQ : quantum[i]
        quantum[i] = quantum[i] > maxErr ? maxErr : quantum[i]
        firstguess = updateQ(Val(O), i, x, q, quantum, exactA, d, cacheA, dxaux, qaux, tx, tq, simt, ft, nextStateTime)
        tx[i] = simt
        tq[i] = simt
        Liqss_reComputeNextTime(Val(O), i, simt, nextStateTime, x, q, quantum)
      end
      computeNextEventTime(Val(O), index, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ) #update zcf before thiscatch in qss quantizer to avoid infinite events
      for j in (HD[modifiedIndex]) # care about dependency to this event only
        elapsedx = simt - tx[j]
        if elapsedx > 0
          x[j].coeffs[1] = x[j](elapsedx)
          tx[j] = simt
        end #= @show j,x[j] =#
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
            end #= @show q[b] =#
          end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(j, -1, -1, q, d, t, taylorOpsCache)
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
        f(-1, j, -1, q, d, t, taylorOpsCache)  # run ZCF--------  
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
      end
    end#end state/input/event
    savetime = simt + savetimeincrement
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
  stats=Stats(totalSteps,simulStepCount,evCount,numSteps)
  createSol(Val(T), Val(O), savedTimes, savedVars, "nmliqss$O", string(odep.prname), absQ, stats, ft)
end#end integrate