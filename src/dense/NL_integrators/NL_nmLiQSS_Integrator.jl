function integrate(Al::QSSAlgorithm{:nmliqss,O}, CommonqssData::CommonQSS_data{0}, liqssdata::LiQSS_data{O,false}, odep::NLODEProblem{PRTYPE,T,0,0,CS}, f::Function, jac::Function, SD::Function, exactA::Function) where {PRTYPE,CS,O,T}
  cacheA = liqssdata.cacheA
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
  taylorOpsCache = CommonqssData.taylorOpsCache#cacheSize=odep.cacheSize
  qaux = liqssdata.qaux
  dxaux = liqssdata.dxaux #= olddx=liqssdata.olddx; ; olddxSpec=liqssdata.olddxSpec =#
  d = [0.0]# this is a dummy var used in updateQ and simulUpdate because in the discrete world exactA needs d, this is better than creating new updateQ and simulUpdate functions
  pp = pointer(Vector{NTuple{2,Float64}}(undef, 7))
  respp = pointer(Vector{Float64}(undef, 6))
  acceptedi = Vector{Vector{Float64}}(undef, 4 * O - 1)
  acceptedj = Vector{Vector{Float64}}(undef, 4 * O - 1)
  cacheRootsi = zeros(8 * O - 4) #vect of floats to hold roots for simul_analytic 
  cacheRootsj = zeros(8 * O - 4)
  for i = 1:4*O-1 #3
    acceptedi[i] = [0.0, 0.0]#zeros(2)
    acceptedj[i] = [0.0, 0.0]#zeros(2)
  end
  exactA(q, d, cacheA, 1, 1, initTime + 1e-9)
  trackSimul = Vector{Int}(undef, 1)
  numSteps = Vector{Int}(undef, T)
  #######################################compute initial values##################################################
  n = 1
  for k = 1:O # compute initial derivatives for x and q (similar to a recursive way )
    n = n * k
    for i = 1:T
      q[i].coeffs[k] = x[i].coeffs[k]  # q computed from x and it is going to be used in the next x
    end
    for i = 1:T
      clearCache(taylorOpsCache, Val(CS), Val(O))
      f(i, q, t, taylorOpsCache)
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
      # #@show i
      computeNextTime(Val(O), i, initTime, nextInputTime, x, quantum)
      if nextInputTime[i] == Inf
        # #@show nextInputTime
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(i, q, t + smallAdvance, taylorOpsCache)
        computeNextInputTime(Val(O), i, initTime, smallAdvance, taylorOpsCache[1], nextInputTime, x, quantum)
        if nextInputTime[i] > initTime + 2 * smallAdvance
          nextInputTime[i] = initTime + 2 * smallAdvance
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
  simulStepCount = 0
  totalSteps = 0
  inputSteps = 0
  printonce = 0
  while simt < ft && totalSteps < maxiters
    if totalSteps == maxiters - 1
      @warn("The algorithm nmliqss$O is taking too long to converge. The simulation will be stopped. Consider using a different algorithm!")
    end
    sch = updateScheduler(Val(T), nextStateTime, nextEventTime, nextInputTime)
    simt = sch[2]
    index = sch[1]
    if simt > ft
      break # 
    end
    numSteps[index] += 1
    totalSteps += 1
    t[0] = simt
    ##########################################state########################################
    if sch[3] == :ST_STATE
      xitemp = x[index][0]
      elapsed = simt - tx[index]
      integrateState(Val(O), x[index], elapsed)
      tx[index] = simt
      quantum[index] = relQ * abs(x[index].coeffs[1])
      quantum[index] = quantum[index] < absQ ? absQ : quantum[index]
      quantum[index] = quantum[index] > maxErr ? maxErr : quantum[index]
      dirI = x[index][0] - xitemp
      for b in (jac(index))    # update Qb : to be used to calculate exacte Aindexb...move below updateQ
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
        aij = cacheA[1]# 
        cacheA[1] = 0.0
        exactA(q, d, cacheA, j, index, simt)
        aji = cacheA[1]
        if j != index && aij * aji != 0.0
          if nmisCycle_and_simulUpdate(cacheRootsi, cacheRootsj, acceptedi, acceptedj, aij, aji, respp, pp, trackSimul, Val(O), index, j, dirI, firstguess, x, q, quantum, exactA, d, cacheA, dxaux, qaux, tx, tq, simt, ft)
            simulStepCount += 1
            clearCache(taylorOpsCache, Val(CS), Val(O))
            f(index, q, t, taylorOpsCache)
            computeDerivative(Val(O), x[index], taylorOpsCache[1])
            for k in SD(j)  #j influences k
              if k != index && k != j
                elapsedx = simt - tx[k]
                x[k].coeffs[1] = x[k](elapsedx)
                tx[k] = simt
                elapsedq = simt - tq[k]
                if elapsedq > 0 integrateState(Val(O - 1), q[k], elapsedq) ;tq[k] = simt end
                for b in (jac(k))
                  elapsedq = simt - tq[b]
                  if elapsedq > 0
                    integrateState(Val(O - 1), q[b], elapsedq)
                    tq[b] = simt
                  end
                end
                clearCache(taylorOpsCache, Val(CS), Val(O))
                f(k, q, t, taylorOpsCache)
                computeDerivative(Val(O), x[k], taylorOpsCache[1])
                Liqss_reComputeNextTime(Val(O), k, simt, nextStateTime, x, q, quantum)
              end#end if k!=0
            end#end for k depend on j          
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
        if elapsedx > 0 x[c].coeffs[1] = x[c](elapsedx) ;tx[c] = simt end # 
        elapsedq = simt - tq[c]
        if elapsedq > 0 integrateState(Val(O - 1), q[c], elapsedq) ;tq[c] = simt end   # c never been visited 
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(c, q, t, taylorOpsCache)
        computeDerivative(Val(O), x[c], taylorOpsCache[1])
        Liqss_reComputeNextTime(Val(O), c, simt, nextStateTime, x, q, quantum)
      end#end for SD
    ##################################input########################################
    elseif sch[3] == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself    
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
      f(index, q, t, taylorOpsCache)
      computeDerivative(Val(O), x[index], taylorOpsCache[1])
      computeNextTime(Val(O), index, simt, nextInputTime, x, quantum) #
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
        end#q needs to be updated here for recomputeNext                 
        for b in jac(j)
          elapsedq = simt - tq[b]
          if elapsedq > 0 integrateState(Val(O - 1), q[b], elapsedq) ;tq[b] = simt
          end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(j, q, t, taylorOpsCache)
        computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end#end for
      #################################################################event########################################
    end#end state/input/event
    push!(savedVars[index], (x[index][0] + q[index][0]) / 2)
    push!(savedTimes[index], simt)
  end#end while
  stats=Stats(totalSteps,simulStepCount,evCount,numSteps)
  createSol(Val(T), Val(O), savedTimes, savedVars, "nmliqss$O", string(odep.prname), absQ, stats, ft)
end
