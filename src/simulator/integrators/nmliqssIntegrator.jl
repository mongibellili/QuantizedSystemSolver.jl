"""
    integrate(alg::QSSAlgorithm{:nmliqss,O}, commonQssData::CommonQSS_Data{0}, liqssdata::LiQSS_Data{O,M}, odep::ODEProblemData{F,PRTYPE,T,D,0,CS}, f::Function, jac::Function, SD::Function, exactA::Function) where {F,PRTYPE,CS,O,T,D,M}

Integrates a nonlinear ordinary differential equation (ODE) problem (without events) using the nmLiqss (modified Liqss that detect events) algorithm.

# Arguments
- `Al::QSSAlgorithm{:nmliqss,O}`: The QSS algorithm to be used for integration.
- `commonQssData::CommonQSS_Data{0}`: Common data structure for QSS algorithms.
- `liqssdata::LiQSS_Data{O,M}`: Data specific to the LiQSS algorithm.
- `odep::ODEProblemData{F,PRTYPE,T,0,0,CS}`: The nonlinear ODE problem to be solved.
- `f::Function`: The function defining the ODE system.
- `jac::Function`: The Jacobian dependency function of the ODE system.
- `SD::Function`: The state derivative dependency function.
- `exactA::Function`: The exact jacobian expression function for the ODE system.

# Returns
- A solution 

"""
function integrate(alg::QSSAlgorithm{:nmliqss,O}, commonQssData::CommonQSS_Data{0}, liqssdata::LiQSS_Data{O,M}, odep::ODEProblemData{F,PRTYPE,T,D,0,CS}, f::Function, jac::Function, SD::Function, exactA::Function) where {F,PRTYPE,CS,O,T,D,M}
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


  #********************************helper values*******************************  
  
  numStateSteps = Vector{Int}(undef, T)
  numInputSteps = Vector{Int}(undef, T)

  # call once (performance)
  exactA(q, d, cacheA, 1, 1, initTime + 1e-9,clF)
  f(1, q, d, t, taylorOpsCache,clF)
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
     f(i, q, d, t,taylorOpsCache,clF)
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
      #= computeNextTime(Val(O), i, initTime, nextInputTime, x, quantum)
      if nextInputTime[i] == Inf
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(i, q, t , d,taylorOpsCache,clF)
        #computeNextInputTime(Val(O), i, initTime, smallAdvance, taylorOpsCache[1], nextInputTime, x, quantum)
        clearCache(taylorOpsCache, Val(CS), Val(O)) =#
        computeNextInputTime(Val(O), i, t, f,clF,d, taylorOpsCache, nextInputTime, x,q, quantum) 
       #=  if nextInputTime[i] > initTime + 2 * smallAdvance
          nextInputTime[i] = initTime + 2 * smallAdvance
        end  =#
      #end
    else
      #= cacheA[1] = 0.0
      exactA(q, d, cacheA, i, i, initTime+1e-12,clF)
      a = cacheA[1] =#
      aii=prepareAii(i, i,a,exactA, q, d, cacheA, initTime+1e-12, clF) # compute Aii for the first time
      updateQInit(Val(O), i, x, q, quantum,aii, dxaux, qaux, tx, tq, initTime + 1e-12, ft, nextStateTime) #1e-9 exactAfunc contains 1/t???
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
  inpuStepCount=0

  while t[0] < ft && totalSteps < maxiters
 
    if totalSteps == maxiters - 1
      @warn("The algorithm nmliqss$O reached max iterations. The simulation will be stopped. Consider using a different algorithm or a different cycle detection mechanism!")
    end
    sch = updateScheduler(Val(T), nextStateTime, nextEventTime, nextInputTime)
    simt = sch[2]
    i = sch[1]
    if simt > ft
      break # 
    end
    
    totalSteps += 1
    t[0] = simt

   # @show x[i],q[i],quantum[i],nextStateTime[i]

    ##########################################state########################################
    if sch[3] == :ST_STATE
      numStateSteps[i] += 1
      xitemp = x[i][0]
      elapsed = simt - tx[i]
      integrateState(Val(O), x[i], elapsed)
      integrateOlddx(Val(O),i,x,tx,simt,olddx)
      tx[i] = simt
      quantum[i] = relQ * abs(x[i].coeffs[1])
      quantum[i] = quantum[i] < absQ ? absQ : quantum[i]
      quantum[i] = quantum[i] > maxErr ? maxErr : quantum[i]
      dirI = x[i][0] - xitemp
      for b in (jac(i))    # update Qb : to be used to calculate exacte Aindexb...move below updateQ
        elapsedq = simt - tq[b]
        if elapsedq > 0
          integrateState(Val(O - 1), q[b], elapsedq)
          tq[b] = simt
        end
      end
     #=  cacheA[1] = 0.0
      exactA(q, d, cacheA, i, i, simt,clF)
      a = cacheA[1] =#
      aii=prepareAii(i, i,a,exactA, q, d, cacheA, simt, clF) # compute Aii for the first time
      updateQ(Val(O), i, x, q, quantum, aii, dxaux, qaux, tx, tq, simt, ft, nextStateTime)
      tq[i] = simt
      #----------------------------------------------------check dependecy cycles---------------------------------------------   
      trackSimul[1] = 0
      for j in SD(i)
       #=  if trackSimul[1] != j
          @show j
        end =#
          for b in (jac(j))    # update Qb: to be used to calculate exacte Ajb
            elapsedq = simt - tq[b]
            if elapsedq > 0
              integrateState(Val(O - 1), q[b], elapsedq)
              tq[b] = simt
            end
          end
         #=  cacheA[1] = 0.0
          exactA(q, d, cacheA, i, j, simt,clF)
          aij = cacheA[1]#  =#
          aij=prepareAii(i, j,a,exactA, q, d, cacheA, simt, clF)
          #= cacheA[1] = 0.0
          exactA(q, d, cacheA, j, i, simt,clF)# can have clF also
          aji = cacheA[1] =#
          aji=prepareAii(j, i,a,exactA, q, d, cacheA, simt, clF)
          ajj=prepareAii(j, j,a,exactA, q, d, cacheA, simt, clF)

       #=    if 0.000534422899295716<=simt<=0.00053442289929572
            @show simt,i,j
            @show  aii  , aij 
            @show  aji  , ajj 
          end =#
         # @show j,aij,aji
          if j != i && aij * aji != 0.0
            integrateOlddx(Val(O),j,x,tx,simt,olddx)
            if isCycle_simulUpdate(aii,ajj,aij, aji, trackSimul, Val(O),Val(M), i, j, dirI, x, q, quantum, dxaux, qaux, tx, tq, simt, ft)
              
              simulStepCount += 1
              clearCache(taylorOpsCache, Val(CS), Val(O))
              f(i, q, d, t,taylorOpsCache,clF)
              computeDerivative(Val(O), x[i], taylorOpsCache[1])
              Liqss_reComputeNextTime(Val(O), i, simt, nextStateTime, x, q, quantum)
              for k in SD(j)  #j influences k
                if k != i && k != j
                  elapsedx = simt - tx[k]
                  x[k].coeffs[1] = x[k](elapsedx)
                  integrateOlddx(Val(O),k,x,tx,simt,olddx)
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
                  f(k, q, d, t,taylorOpsCache,clF)
                  computeDerivative(Val(O), x[k], taylorOpsCache[1])
                  Liqss_reComputeNextTime(Val(O), k, simt, nextStateTime, x, q, quantum)
                  updateOtherApprox(k,j,x,q,a,qaux,olddx,simt)
                 # @show simt,k,x[k],q[k],quantum[k],nextStateTime[k]
                end#end if k!=0
              end#end for k depend on j   
              updateLinearApprox(j,x,q,a,qaux,olddx,simt)   
             #=  @show i,x[i],q[i],quantum[i],nextStateTime[i]
              @show j,x[j],q[j],quantum[j],nextStateTime[j]   =#    
            end#end ifcycle check
          end#end if j!=i
        #end#end if trackSimul[1] != j
      end#end FOR_cycle check
     #=  if trackSimul[1] != 0  #qi changed after throw
        Liqss_reComputeNextTime(Val(O), i, simt, nextStateTime, x, q, quantum)
      end =#
     #=  if trackSimul[1] > 1
        simulStepCount += 1
      end =#
      #-------------------------------------------------------------------------------------
      #---------------------------------normal liqss: proceed--------------------------------
      #-------------------------------------------------------------------------------------
      for c in SD(i)   #i influences c       
        elapsedx = simt - tx[c]
        if elapsedx > 0 x[c].coeffs[1] = x[c](elapsedx) ;integrateOlddx(Val(O),c,x,tx,simt,olddx);tx[c] = simt end # 
        elapsedq = simt - tq[c]
        if elapsedq > 0 integrateState(Val(O - 1), q[c], elapsedq) ;tq[c] = simt end   # c never been visited 
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(c, q, d, t,taylorOpsCache,clF)
       
        computeDerivative(Val(O), x[c], taylorOpsCache[1])
        Liqss_reComputeNextTime(Val(O), c, simt, nextStateTime, x, q, quantum)
        updateOtherApprox(c,i,x,q,a,qaux,olddx,simt)
      end#end for SD
      updateLinearApprox(i,x,q,a,qaux,olddx,simt)
    ##################################input########################################
    elseif sch[3] == :ST_INPUT  # time of change has come to a state var that does not depend on anything...
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
      #= cacheA[1] = 0.0
      exactA(q, d, cacheA, i, i, simt,clF)
      a = cacheA[1] =#
      aii=prepareAii(i, i,a,exactA, q, d, cacheA, simt, clF) # compute Aii for the first time
      updateQ(Val(O), i, x, q, quantum, aii, dxaux, qaux, tx, tq, simt, ft, nextStateTime)
      tq[i] = simt
      #computeNextTime(Val(O), i, simt, nextStateTime, x, quantum) # 
      clearCache(taylorOpsCache, Val(CS), Val(O))
      f(i, q, d, t,taylorOpsCache,clF)
      computeDerivative(Val(O), x[i], taylorOpsCache[1])
      #computeNextTime(Val(O), i, simt, nextInputTime, x, quantum) #v14
      #if nextInputTime[i] > simt + 2 * elapsed
       # clearCache(taylorOpsCache, Val(CS), Val(O))
       # f(i, q, t + 0.003*smallAdvance, d,taylorOpsCache,clF)
        computeNextInputTime(Val(O), i, t, f,clF,d, taylorOpsCache, nextInputTime, x,q, quantum)#v13
       # reComputeNextTime(Val(O), i, simt, nextInputTime, x, q, quantum) #v12
      # updateQ(Val(O), i, x, q, quantum, exactA, d, cacheA, dxaux, qaux, tx, tq, simt, ft, nextInputTime,clF)# 
      #  Liqss_reComputeNextTime(Val(O), i, simt, nextInputTime, x, q, quantum)#v11 
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
        f(j, q, d, t,taylorOpsCache,clF)
        computeDerivative(Val(O), x[j], taylorOpsCache[1])
       
        Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
        updateOtherApprox(j,i,x,q,a,qaux,olddx,simt)
      end#end for j
      updateLinearApprox(i,x,q,a,qaux,olddx,simt)
    end#end state/input/event
   # push!(savedVars[i], (x[i][0] + q[i][0]) / 2)
    push!(savedVars[i], (x[i][0] ) )
    push!(savedTimes[i], simt)
  end#end while
  stats=Stats(totalSteps,simulStepCount,inpuStepCount,0,numStateSteps,numInputSteps)
  createSol(Val(T), Val(O), savedTimes, savedVars, toString(alg), string(odep.prname), absQ, stats, ft)
end
