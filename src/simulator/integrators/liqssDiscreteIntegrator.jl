"""
    integrate(alg::QSSAlgorithm{:liqss,O}, commonQssData::CommonQSS_Data{Z}, liqssdata::LiQSS_Data{O,CDM}, odep::ODEProblemData{JACMODE,T,D,Z,CS,F,JAC,CLS}, f::Function, jac::Function, SD::Function, exactA::Function) where {JACMODE,T,D,Z,CS,F,JAC,CLS,CDM}

Integrates a nonlinear ordinary differential equation (ODE) problem with `events` using the LiQSS (Linearized Quantized State System) algorithm.

# Arguments
- `Al::QSSAlgorithm{:liqss,O}`: The QSS algorithm type for LiQSS.
- `commonQssData::CommonQSS_Data{Z}`: Common data structure for QSS algorithms.
- `liqssdata::LiQSS_Data{O,1}`: Data structure specific to the LiQSS algorithm.
- `odep::ODEProblemData{JACMODE,T,D,Z,CS,F,JAC,CLS}`: The nonlinear ODE problem to be solved.
- `f::Function`: The function defining the ODE system.
- `jac::Function`: The Jacobian dependency function of the ODE system.
- `SD::Function`: The state derivative dependency function.
- `exactA::Function`: The exact jacobian expression function for the ODE system.

# Returns
- A solution after the integration process.

"""
function integrate(alg::QSSAlgorithm{:liqss,O}, commonQssData::CommonQSS_Data{Z,PT},  odep::ODEProblemData{JACMODE,T,D,Z,CS,F,JAC,CLS}, f::F, jac::JAC_CLS, SD::SD_CLS,liqssdata::LiQSS_Data{O,CDM}, exactA::JAC) where {O,JACMODE,T,D,Z,CS,F,JAC,CLS,PT,JAC_CLS,SD_CLS,CDM}

  VERBOSE,ft,initTime,relQ,absQ,relZ,absZ,maxErr,maxiters,quantum,nextStateTime,nextEventTime,nextInputTime,tx,tq,x,q,t,savedVars,savedTimes,
  taylorOpsCache,d,numStateSteps,numInputSteps,oldsignValue,clF,zc_SimpleJac,HZ,HD,SZ,evDep = initIntegrator(Val(O), Val(T), Val(CS), odep, commonQssData, f)
  f(1, -1, -1, q, d, t, taylorOpsCache,clF) # call once (for performance)
  # implicit integrator data
  a=liqssdata.a; cacheA = liqssdata.cacheA; qaux = liqssdata.qaux; dxaux = liqssdata.dxaux; olddx = liqssdata.olddx
  exactA(q, d, cacheA, 1, 1, initTime + 1e-8,clF) # call once (for performance)
  for i = 1:T
    numStateSteps[i] = 0; numInputSteps[i] = 0
    push!(savedVars[i], x[i][0]); push!(savedTimes[i], initTime)
    quantum[i] = 2*relQ * abs(x[i].coeffs[1]);quantum[i] = quantum[i] < absQ ? 2*absQ : quantum[i];#quantum[i] = quantum[i] > maxErr ? maxErr : quantum[i]
    if isempty(jac(i))
        computeNextInputTime(Val(O), i, t, f,clF,d, taylorOpsCache, nextInputTime, x,q, quantum) 
    else
      aii=prepareAii(i,i,a,exactA, q, d, cacheA, initTime+1e-12, clF) # compute Aii for the first time
      updateQInit(Val(O), i, x, q, quantum, aii, dxaux, qaux, tx, tq, initTime + 1e-12, ft, nextStateTime) 
    end
  end
  #---------------------------------------------------------------------------------while loop------------------------------------------------------------------------- 
  simt = initTime
  totalSteps = 0;inpuStepCount=0;simulStepCount = 0;evCount = 0;rejectedEvCount = 0
  modifiedIndex = 0
  while t[0] < ft && totalSteps < maxiters
    if totalSteps == maxiters - 1 @warn("The algorithm nmliqss$O reached max iterations. The simulation will be stopped. Consider using a different algorithm!") end
    sch = updateScheduler(Val(T), nextStateTime, nextEventTime, nextInputTime)
    i = sch[1]; simt = sch[2]; stepType = sch[3]
    simt > ft && break   
    totalSteps += 1
    t[0] = simt
    ##########################################state######################################## 
    if stepType == :ST_STATE
      numStateSteps[i] += 1
      integrateState(Val(O), x[i], tx[i], simt) ;tx[i]=simt
      integrateOlddx(Val(O),i,x,tx,simt,olddx)
      quantum[i] = 2*relQ * abs(x[i].coeffs[1]); quantum[i] = quantum[i] < absQ ? 2*absQ : quantum[i]; #quantum[i] = quantum[i] > maxErr ? maxErr : quantum[i]
      for b in (jac(i))    # update Qb : to be used to calculate exact aii in case symbolic jac is used.
          integrateState(Val(O - 1),q[b],tq[b],simt) ;tq[b]=simt
      end
      aii=prepareAii(i,i,a,exactA, q, d, cacheA, initTime, clF) 
      updateQ(Val(O), i, x, q, quantum,aii,  dxaux, qaux, tx, tq, simt, ft, nextStateTime); tq[i] = simt
      #---------------------------------normal liqss: proceed--------------------------------
      #-------------------------------------------------------------------------------------
      for j in SD(i)   #i influences j  
        elapsedx = simt - tx[j]
        if elapsedx > 0
          x[j].coeffs[1] = x[j](elapsedx)
          integrateOlddx(Val(O),j,x,tx,simt,olddx)
          tx[j] = simt
        end # 
        integrateState(Val(O - 1),q[j],tq[j],simt) ; tq[j]=simt  # if j never been visited 
          computeDerivatives(Val(O),Val(CS),f,j,q,tq,simt,d,t,taylorOpsCache,clF,x,jac) 
        Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
        updateOtherApprox(j,i,x,q,a,qaux,olddx,simt)   # in case jac approximation mode is used.
      end#end for SD
      for j in (SZ[i])
        updateZCF(Val(O),Val(CS),f,j,q,tq,simt,d,t,taylorOpsCache,clF,zc_SimpleJac)
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absZ,relZ)
      end#end for SZ
      updateLinearApprox(i,x,q,a,qaux,olddx,simt)  # in case jac approximation mode is used.
    
    ##################################input########################################
    elseif stepType == :ST_INPUT  # time of change has come to a state var that does not depend on any variable..
      inpuStepCount+=1
      numInputSteps[i] += 1
      integrateState(Val(O), x[i], tx[i], simt) ;tx[i] = simt
      integrateOlddx(Val(O),i,x,tx,simt,olddx)
      quantum[i] = 2*relQ * abs(x[i].coeffs[1]); quantum[i] = quantum[i] < absQ ? 2*absQ : quantum[i]; quantum[i] = quantum[i] > maxErr ? maxErr : quantum[i]
      aii=prepareAii(i,i,a,exactA, q, d, cacheA, initTime, clF) 
      updateQ(Val(O), i, x, q, quantum,aii,  dxaux, qaux, tx, tq, simt, ft, nextStateTime)  ;tq[i] = simt
      clearCache(taylorOpsCache, Val(CS), Val(O)); f(i, -1, -1, q, d, t, taylorOpsCache,clF)
      computeDerivative(Val(O), x[i], taylorOpsCache[1])
      computeNextInputTime(Val(O), i, t, f,clF,d, taylorOpsCache, nextInputTime, x,q, quantum)
      for j in (SD(i))
        elapsedx = simt - tx[j]
        if elapsedx > 0 
          x[j].coeffs[1] = x[j](elapsedx) # der later
          integrateOlddx(Val(O),j,x,tx,simt,olddx)
          tx[j] = simt 
        end
        integrateState(Val(O - 1),q[j],tq[j],simt)  ;tq[j] = simt #q needs to be updated here for recomputeNext                 
        computeDerivatives(Val(O),Val(CS),f,j,q,tq,simt,d,t,taylorOpsCache,clF,x,jac)
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
        updateOtherApprox(j,i,x,q,a,qaux,olddx,simt)
      end#end for
      for j in (SZ[i]) 
         updateZCF(Val(O),Val(CS),f,j,q,tq,simt,d,t,taylorOpsCache,clF,zc_SimpleJac)
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absZ,relZ)
      end#end for SZ
      updateLinearApprox(i,x,q,a,qaux,olddx,simt)
    #################################################################event########################################
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
        if abs(oldsignValue[i, 2]) <= absZ#*1e-6 # if oldsignValue very small, reject event , just already handled     
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
        println("event rejected because zcf did not change")
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
        push!(savedVars[k], q[k][0]) # save the variables that are affected by this event; (var=...)
        push!(savedTimes[k], simt)
      dir=sign(q[k][0]-x[k][0]) # direction of event

        x[k][0]= q[k][0] ; tx[k] = simt # update x[k] with q[k] after event
        quantum[k] = 2*relQ * abs(x[k][0]); quantum[k] = quantum[k] < absQ ? 2*absQ : quantum[k]; #quantum[k] = quantum[k] > maxErr ? maxErr : quantum[k]
        #aii=prepareAii(k,k,a,exactA, q, d, cacheA, initTime, clF) 
        #updateQ(Val(O), k, x, q, quantum,aii,  dxaux, qaux, tx, tq, simt, ft, nextStateTime);  tq[k] = simt
        q[k][0]=x[k][0]+dir*quantum[k]  # update q[k] after event # here du[k] is outdated so  no need to call updateQ method. in fact, it may leave q as old value if old du was 0
      end
     computeNextEventTime(Val(O), i, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absZ,relZ) # in case zcf(i) depends on time (or unforseen case), its nexteventtime needs to be updated (without this, there is a possibililty nexteventtime never get updated). Warning: do not move this at the bottom, as taylorOpsCache will hold something else.
     for j in (HD[modifiedIndex]) # care about dependency to this event only
        elapsedx = simt - tx[j] ; if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx); tx[j] = simt end 
        integrateState(Val(O - 1),q[j],tq[j],simt)  ;tq[j] = simt #qj needs to be updated here for recomputeNext (xj updated above)            
        computeDerivatives(Val(O),Val(CS),f,j,q,tq,simt,d,t,taylorOpsCache,clF,x,jac)
        Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)  
      
      end
      for j in (HZ[modifiedIndex])
         updateZCF(Val(O),Val(CS),f,j,q,tq,simt,d,t,taylorOpsCache,clF,zc_SimpleJac)
        computeNextEventTime(Val(O), j, taylorOpsCache[1], oldsignValue, simt, nextEventTime, quantum, absZ,relZ)
      end
     
    end#end state/input/event
    if stepType != :ST_EVENT
      push!(savedVars[i], (x[i][0]+q[i][0])/2)
      push!(savedTimes[i], simt)
    else
     # for j in (HD[modifiedIndex])  # save also the variables that are affected by this event but that are outside the event (ie not in evCont)
       # push!(savedVars[j], (x[j][0]+q[j][0])/2)
       # push!(savedTimes[j], simt)
      #end
    end
  end#end while
    for i = 1:T
      elapsedx = ft- tx[i]  ;  
      computeDerivatives(Val(O),Val(CS),f,i,q,tq,simt,d,t,taylorOpsCache,clF,x,jac) 
      x[i][0] = x[i](elapsedx);
      push!(savedVars[i], x[i][0]); push!(savedTimes[i], ft) # final point
    end
     
   stats=Stats(totalSteps,simulStepCount,inpuStepCount,evCount,rejectedEvCount,numStateSteps,numInputSteps)
   createSol(Val(T), Val(O), savedTimes, savedVars, toString(alg), string(odep.prname), absQ, stats, ft)
end#end integrate