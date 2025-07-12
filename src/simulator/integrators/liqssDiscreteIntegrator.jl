"""
    integrate(alg::QSSAlgorithm{:liqss,O}, commonQssData::CommonQSS_Data{Z}, liqssdata::LiQSS_Data{O,M}, odep::ODEProblemData{F,PRTYPE,T,D,Z,CS}, f::Function, jac::Function, SD::Function, exactA::Function) where {F,PRTYPE,O,T,D,Z,CS,M}

Integrates a nonlinear ordinary differential equation (ODE) problem with `events` using the LiQSS (Linearized Quantized State System) algorithm.

# Arguments
- `Al::QSSAlgorithm{:liqss,O}`: The QSS algorithm type for LiQSS.
- `commonQssData::CommonQSS_Data{Z}`: Common data structure for QSS algorithms.
- `liqssdata::LiQSS_Data{O,1}`: Data structure specific to the LiQSS algorithm.
- `odep::ODEProblemData{F,PRTYPE,T,D,Z,CS}`: The nonlinear ODE problem to be solved.
- `f::Function`: The function defining the ODE system.
- `jac::Function`: The Jacobian dependency function of the ODE system.
- `SD::Function`: The state derivative dependency function.
- `exactA::Function`: The exact jacobian expression function for the ODE system.

# Returns
- A solution after the integration process.

"""
function integrate(alg::QSSAlgorithm{:liqss,O}, commonQssData::CommonQSS_Data{Z}, liqssdata::LiQSS_Data{O,M}, odep::ODEProblemData{F,PRTYPE,T,D,Z,CS}, f::Function, jac::Function, SD::Function, exactA::Function) where {F,PRTYPE,O,T,D,Z,CS,M}
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
  taylorOpsCache = commonQssData.taylorOpsCache#
  d = commonQssData.d 
  clF=odep.closureFuncs[1]

  # implicit integrator data
  a=liqssdata.a
  cacheA = liqssdata.cacheA
  qaux = liqssdata.qaux
  dxaux = liqssdata.dxaux
  olddx = liqssdata.olddx
  # discrete prob data
  zc_SimpleJac = odep.ZCjac
  HZ = odep.HZ
  HD = odep.HD
  SZ = odep.SZ
  evDep = odep.eventDependencies
  #@show typeof(olddx)

  #********************************helper values*******************************  
  oldsignValue = MMatrix{Z,2}(zeros(Z * 2))  #usedto track if zc changed sign; each zc has a value and a sign 
  numStateSteps = Vector{Int}(undef, T)
  numInputSteps = Vector{Int}(undef, T)

  # call once (performance)
  exactA(q, d, cacheA, 1, 1, initTime + 1e-9,clF)
  f(1, -1, -1, q, d, t, taylorOpsCache,clF)
  # track  simul_var update ...can add: allow user to pick coupled vars
  trackSimul = Vector{Int}(undef, 1)
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
  #smallAdvance = ft / 1000
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
        discrete_computeNextInputTime(Val(O), i, t, f,clF,d, taylorOpsCache, nextInputTime, x,q, quantum) 
    else
      aii=prepareAii(i,i,a,exactA, q, d, cacheA, initTime+1e-12, clF) # compute Aii for the first time
      updateQInit(Val(O), i, x, q, quantum, aii, dxaux, qaux, tx, tq, initTime + 1e-12, ft, nextStateTime) 
    end
  end
  for i = 1:Z
    clearCache(taylorOpsCache, Val(CS), Val(O))
     f(-1, i, -1, x, d, t, taylorOpsCache,clF)
    oldsignValue[i, 2] = taylorOpsCache[1][0] #value
    oldsignValue[i, 1] = sign(taylorOpsCache[1][0]) #sign modify 
    computeNextEventTime(Val(O), i, taylorOpsCache[1], oldsignValue, initTime, nextEventTime, quantum, absQ)
  end
  ####################################################################################################################################################################
  #---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
  ###################################################################################################################################################################
  
  simt = initTime
  totalSteps = 0
  modifiedIndex = 0
  evCount = 0
  simulStepCount = 0
  inpuStepCount=0

  while t[0] < ft && totalSteps < maxiters
   
    if totalSteps == maxiters - 1 @warn("The algorithm nmliqss$O reached max iterations. The simulation will be stopped. Consider using a different algorithm or a different cycle detection mechanism!") end
    sch = updateScheduler(Val(T), nextStateTime, nextEventTime, nextInputTime)
    simt = sch[2]
    i = sch[1]
    stepType = sch[3]
  
    if simt > ft
    #=   if VERBOSE 
        println("simulation ended ") 
        println("scheduler status: ", nextStateTime[i], " ", nextEventTime[i], " ", nextInputTime[i])
      end =#

      break   
    end
#=     if totalSteps> 3975
      @show simt, t[0],  nextStateTime[2]
     end =#
    totalSteps += 1
    t[0] = simt
    #@show x[i],q[i],quantum[i],nextStateTime[i]
    ##########################################state######################################## 
    if stepType == :ST_STATE
      #= if totalSteps> 3975
        @show simt, i, x[i],q[i] ,t[0], nextStateTime[i]
      end =#
     #=  if 23.5244<simt< 23.5244681977
        @show totalSteps,simt,i, t[0],x[i],q[i],  nextStateTime
       end  =#
      xitemp = x[i][0]
      numStateSteps[i] += 1
      elapsed = simt - tx[i]
      integrateState(Val(O), x[i], elapsed)
      integrateOlddx(Val(O),i,x,tx,simt,olddx)
      tx[i] = simt
      dirI = x[i][0] - xitemp
      quantum[i] = relQ * abs(x[i].coeffs[1])
      quantum[i] = quantum[i] < absQ ? absQ : quantum[i]
      quantum[i] = quantum[i] > maxErr ? maxErr : quantum[i]
      for b in (jac(i))    # update Qb : to be used to calculate exact Aib
        elapsedq = simt - tq[b]
        if elapsedq > 0
          integrateState(Val(O - 1), q[b], elapsedq)
          tq[b] = simt
        end
      end
      aii=prepareAii(i,i,a,exactA, q, d, cacheA, initTime, clF) 
      firstguess = updateQ(Val(O), i, x, q, quantum,aii,  dxaux, qaux, tx, tq, simt, ft, nextStateTime)

      tq[i] = simt
      #---------------------------------normal liqss: proceed--------------------------------
      #-------------------------------------------------------------------------------------

      for c in SD(i)   #i influences c  
        elapsedx = simt - tx[c]
        if elapsedx > 0
          x[c].coeffs[1] = x[c](elapsedx)
          integrateOlddx(Val(O),c,x,tx,simt,olddx)
          tx[c] = simt
        end # 

        elapsedq = simt - tq[c]
        if elapsedq > 0 integrateState(Val(O - 1), q[c], elapsedq) ;tq[c] = simt end   # c never been visited 
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(c, -1, -1, q, d, t, taylorOpsCache,clF)
        computeDerivative(Val(O), x[c], taylorOpsCache[1])
        Liqss_reComputeNextTime(Val(O), c, simt, nextStateTime, x, q, quantum)
        updateOtherApprox(c,i,x,q,a,qaux,olddx,simt)
      end#end for SD
       for j in (SZ[i])
        for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
          elapsedq = simt - tq[b]
          if elapsedq > 0  integrateState(Val(O - 1), q[b], elapsedq) ;tq[b] = simt end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(-1, j, -1, q, d, t, taylorOpsCache,clF)   # run ZCF--------      
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
      end#end for SZ
      updateLinearApprox(i,x,q,a,qaux,olddx,simt)
    ##################################input########################################
    elseif stepType == :ST_INPUT  # time of change has come to a state var that does not depend on anything..
      inpuStepCount+=1
      numInputSteps[i] += 1
      elapsed = simt - tx[i]
      integrateState(Val(O), x[i], elapsed)
      integrateOlddx(Val(O),i,x,tx,simt,olddx)
      tx[i] = simt
      quantum[i] = relQ * abs(x[i].coeffs[1])
      quantum[i] = quantum[i] < absQ ? absQ : quantum[i]
      quantum[i] = quantum[i] > maxErr ? maxErr : quantum[i]
    #=   for k = 1:O
        q[i].coeffs[k] = x[i].coeffs[k]  
      end =#
      aii=prepareAii(i,i,a,exactA, q, d, cacheA, initTime, clF) 
      updateQ(Val(O), i, x, q, quantum,aii,  dxaux, qaux, tx, tq, simt, ft, nextStateTime)
    
      tq[i] = simt
      #computeNextTime(Val(O), i, simt, nextStateTime, x, quantum) #
      clearCache(taylorOpsCache, Val(CS), Val(O))
      f(i, -1, -1, q, d, t, taylorOpsCache,clF)
  
      computeDerivative(Val(O), x[i], taylorOpsCache[1])
  
      #if nextInputTime[i] > simt + 2 * elapsed
        #= clearCache(taylorOpsCache, Val(CS), Val(O))
        f(i, -1, -1, q, d, t + 10*smallAdvance, taylorOpsCache,clF)
        discrete_computeNextInputTime(Val(O), i, simt, 10*smallAdvance, taylorOpsCache[1], nextInputTime, x, quantum) =#
        discrete_computeNextInputTime(Val(O), i, t, f,clF,d, taylorOpsCache, nextInputTime, x,q, quantum)
        #= if nextInputTime[i] > simt + 2 * elapsed
          nextInputTime[i] = simt + 2 * elapsed
        end =#
      #end

      for j in (SD(i))
        elapsedx = simt - tx[j]
        if elapsedx > 0 
          x[j].coeffs[1] = x[j](elapsedx)
          integrateOlddx(Val(O),j,x,tx,simt,olddx)
          tx[j] = simt 
        end
        elapsedq = simt - tq[j]
        if elapsedq > 0 integrateState(Val(O - 1), q[j], elapsedq) ;tq[j] = simt end#q needs to be updated here for recomputeNext                 
        # elapsed update all other vars that this derj depends upon.
        for b in jac(j)
          elapsedq = simt - tq[b]
          if elapsedq > 0 integrateState(Val(O - 1), q[b], elapsedq) ;tq[b] = simt end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(j, -1, -1, q, d, t, taylorOpsCache,clF)
        computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
        updateOtherApprox(j,i,x,q,a,qaux,olddx,simt)
      end#end for
      for j in (SZ[i]) 
        for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
          elapsedq = simt - tq[b]
          if elapsedq > 0 integrateState(Val(O - 1), q[b], elapsedq); tq[b] = simt end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(-1, j, -1, q, d, t, taylorOpsCache,clF)   # run ZCF--------      
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
      end#end for SZ
      updateLinearApprox(i,x,q,a,qaux,olddx,simt)
    #################################################################event########################################
    else
      evCount += 1
   
       for b in zc_SimpleJac[i] # elapsed update all other vars that this zc depends upon.
        elapsedq = simt - tq[b]
        if elapsedq > 0
          integrateState(Val(O - 1), q[b], elapsedq)
          tq[b] = simt
        end
      end
      clearCache(taylorOpsCache, Val(CS), Val(O))
      #@show taylorOpsCache,x[4],q[4]
      f(-1, i, -1, q, d, t, taylorOpsCache,clF)    # run ZCF again to verify-------- 
   #  @show simt, i, oldsignValue[i, 2], taylorOpsCache[1][0]
      if oldsignValue[i, 2] * taylorOpsCache[1][0] >= 0 # if computeNextEvent errored 
        if abs(taylorOpsCache[1][0]) > 1e-9 * absQ # if error is negligeable then ok consider as event, else reject....if both have same sign and zcf is not very small: zc==1e-9*absQ is allowed as an event
          #println("event rejected1: ",simt,"__",oldsignValue[i, 2] ," __ ", taylorOpsCache[1][0])
          computeNextEventTime(Val(O), i, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
          
          continue #event rejected
        end
      end
     #=  if abs(oldsignValue[i, 2]) < 1e-9 * absQ  #earlier zc==1e-9*absQ was considered event , so now it should be prevented from passing
        nextEventTime[i] = Inf # at this instant next zc will be triggered now, and this will lead to infinite events, so cannot computenextevent here
        println("event rejected2: ",oldsignValue[i, 2] ," __ ", taylorOpsCache[1][0])
        continue
      end =#
      if taylorOpsCache[1][0] > oldsignValue[i, 2] #scheduled rise
        modifiedIndex = 2 * i - 1
      elseif taylorOpsCache[1][0] < oldsignValue[i, 2] #scheduled drop
        modifiedIndex = 2 * i
      else # == ( zcf==oldZCF)
        #println("event rejected3: ",oldsignValue[i, 2] ," __ ", taylorOpsCache[1][0])
        computeNextEventTime(Val(O), i, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ)
       
        continue
      end
      #evCount += 1
      oldsignValue[i, 2] = taylorOpsCache[1][0]
      oldsignValue[i, 1] = sign(taylorOpsCache[1][0])
      for b in evDep[modifiedIndex].evContRHS
        elapsedq = simt - tq[b]
        if elapsedq > 0
          integrateState(Val(O - 1), q[b], elapsedq)
          tq[b] = simt
        end
      end
      #@show simt,d
      for i in evDep[modifiedIndex].evCont
        push!(savedVars[i], q[i][0])
        push!(savedTimes[i], simt)
      end
       f(-1, -1, modifiedIndex, x, d, t, taylorOpsCache,clF)# execute event----------------no need to clear cache; events touch vectors directly
      # @show d
      for i in evDep[modifiedIndex].evCont
        push!(savedVars[i], x[i][0])
        push!(savedTimes[i], simt)
        quantum[i] = relQ * abs(x[i].coeffs[1])
        quantum[i] = quantum[i] < absQ ? absQ : quantum[i]
        quantum[i] = quantum[i] > maxErr ? maxErr : quantum[i]
        aii=prepareAii(i,i,a,exactA, q, d, cacheA, initTime, clF) 
        firstguess = updateQ(Val(O), i, x, q, quantum,aii,  dxaux, qaux, tx, tq, simt, ft, nextStateTime)
        tx[i] = simt
        tq[i] = simt
        Liqss_reComputeNextTime(Val(O), i, simt, nextStateTime, x, q, quantum)
      end
      computeNextEventTime(Val(O), i, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absQ) #update zcf before thiscatch in qss quantizer to avoid infinite events
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
      push!(savedVars[i], x[i][0])
      push!(savedTimes[i], simt)
    else
      for j in (HD[modifiedIndex])
        push!(savedVars[j], x[j][0])
        push!(savedTimes[j], simt)
      end
    end
  end#end while
   stats=Stats(totalSteps,simulStepCount,inpuStepCount,evCount,numStateSteps,numInputSteps)
   #@show stats
   createSol(Val(T), Val(O), savedTimes, savedVars, toString(alg), string(odep.prname), absQ, stats, ft)
end#end integrate