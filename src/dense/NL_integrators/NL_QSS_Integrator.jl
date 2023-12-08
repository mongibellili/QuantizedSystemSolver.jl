#using TimerOutputs
function integrate(Al::QSSAlgorithm{:qss,O},CommonqssData::CommonQSS_data{0}, odep::NLODEProblem{PRTYPE,T,0,0,CS},f::Function,jac::Function,SD::Function) where {PRTYPE,O,T,CS}

  ft = CommonqssData.finalTime;initTime = CommonqssData.initialTime;relQ = CommonqssData.dQrel;absQ = CommonqssData.dQmin;maxErr=CommonqssData.maxErr;
  #= savetimeincrement=CommonqssData.savetimeincrement;savetime = savetimeincrement =#
  quantum = CommonqssData.quantum;nextStateTime = CommonqssData.nextStateTime;nextEventTime = CommonqssData.nextEventTime;nextInputTime = CommonqssData.nextInputTime
  tx = CommonqssData.tx;tq = CommonqssData.tq;x = CommonqssData.x;q = CommonqssData.q;t=CommonqssData.t
   savedVars=CommonqssData.savedVars;
  savedTimes=CommonqssData.savedTimes;integratorCache=CommonqssData.integratorCache;taylorOpsCache=CommonqssData.taylorOpsCache;
 
  #a=deepcopy(odep.initJac);
    #********************************helper values*******************************  
 # qaux=CommonqssData.qaux;olddx=CommonqssData.olddx;olddxSpec = zeros(MVector{T,MVector{O,Float64}}) # later can only care about 1st der
 numSteps = Vector{Int}(undef, T)
 savedVarsQ = Vector{Vector{Float64}}(undef, T)  
 for i=1:T
  savedVarsQ[i]=Vector{Float64}()     
end
  #######################################compute initial values##################################################
n=1
for k = 1:O # compute initial derivatives for x and q (similar to a recursive way )
  n=n*k
   for i = 1:T q[i].coeffs[k] = x[i].coeffs[k] end # q computed from x and it is going to be used in the next x
   for i = 1:T
      clearCache(taylorOpsCache,Val(CS),Val(O));f(i,q, t ,taylorOpsCache)
      ndifferentiate!(integratorCache,taylorOpsCache[1] , k - 1)
      x[i].coeffs[k+1] = (integratorCache.coeffs[1]) / n # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
    end
end

for i = 1:T
  numSteps[i]=0
  push!(savedVarsQ[i],q[i][0])
  push!(savedVars[i],x[i][0])
     push!(savedTimes[i],0.0)
  quantum[i] = relQ * abs(x[i].coeffs[1]) ;quantum[i]=quantum[i] < absQ ? absQ : quantum[i];quantum[i]=quantum[i] > maxErr ? maxErr : quantum[i] 
  computeNextTime(Val(O), i, initTime, nextStateTime, x, quantum)
  initSmallAdvance=0.1
  #t[0]=initTime#initSmallAdvance
 # clearCache(taylorOpsCache,Val(CS),Val(O));
  #@timeit "f" 
  f(i,q,t,taylorOpsCache);#@show taylorOpsCache
 computeNextInputTime(Val(O), i, initTime, initSmallAdvance,taylorOpsCache[1] , nextInputTime, x,  quantum)
  #= assignXPrevStepVals(Val(O),prevStepVal,x,i) =#
end

###################################################################################################################################################################
####################################################################################################################################################################
#---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
###################################################################################################################################################################
####################################################################################################################################################################
simt = initTime ;totalSteps=0;prevStepTime=initTime
 # breakloop= zeros(MVector{1,Float64})
 #@timeit "qssintgrateWhile"
  while simt < ft && totalSteps < 200000000   
   #=  if breakloop[1]>5.0
      break
    end =#
     sch = updateScheduler(Val(T),nextStateTime,nextEventTime, nextInputTime)
    simt = sch[2]
   # @timeit "saveLast" 
    #=  if  simt>ft  
      saveLast!(Val(T),Val(O),savedVars, savedTimes,saveVarsHelper,ft,prevStepTime,integratorCache, x)
      break   ###################################################break##########################################
    end =#
    index = sch[1]
    totalSteps+=1
    t[0]=simt
  ##########################################state######################################## 
  if sch[3] == :ST_STATE
    numSteps[index]+=1;
    elapsed = simt - tx[index];integrateState(Val(O),x[index],elapsed);tx[index] = simt 
    quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index]   
    if abs(x[index].coeffs[2])>1e7 quantum[index]=10*quantum[index] end
    for k = 1:O q[index].coeffs[k] = x[index].coeffs[k] end; tq[index] = simt    
    computeNextTime(Val(O), index, simt, nextStateTime, x, quantum) #
    for j in (SD(index))
        elapsedx = simt - tx[j];if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt end
        # quantum[j] = relQ * abs(x[j].coeffs[1]) ;quantum[j]=quantum[j] < absQ ? absQ : quantum[j];quantum[j]=quantum[j] > maxErr ? maxErr : quantum[j]         
        elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],elapsedq);tq[j] = simt  end#q needs to be updated here for recomputeNext        
        for b in (jac(j)  )    
          elapsedq = simt - tq[b]
          if elapsedq>0
            integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt
          end
         end
        clearCache(taylorOpsCache,Val(CS),Val(O));f(j,q,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1])  
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
    end#end for SD
    ##################################input########################################
  elseif sch[3] == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself    
   @show index
    elapsed = simt - tx[index];integrateState(Val(O),x[index],elapsed);tx[index] = simt 
    quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index]   
    for k = 1:O q[index].coeffs[k] = x[index].coeffs[k] end; tq[index] = simt 
      for b in jac(index) 
        elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
      end
    clearCache(taylorOpsCache,Val(CS),Val(O));f(index,q,t,taylorOpsCache)
    computeNextInputTime(Val(O), index, simt, elapsed,taylorOpsCache[1] , nextInputTime, x,  quantum)
    computeDerivative(Val(O), x[index], taylorOpsCache[1])
    for j in(SD(index))  
        elapsedx = simt - tx[j];
        if elapsedx > 0 
          x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt 
          quantum[j] = relQ * abs(x[j].coeffs[1]) ;quantum[j]=quantum[j] < absQ ? absQ : quantum[j];quantum[j]=quantum[j] > maxErr ? maxErr : quantum[j]   
        end
        elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],elapsedq);tq[j] = simt  end#q needs to be updated here for recomputeNext                 
        # elapsed update all other vars that this derj depends upon.
          for b in jac(j) 
            elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
          end
        
        clearCache(taylorOpsCache,Val(CS),Val(O));f(j,q,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
    end#end for
  end#end state/input/event
 #=  if simt > savetime #|| sch[3] ==:ST_EVENT
    save!(Val(O),savedVars , savedTimes , saveVarsHelper,prevStepTime ,simt,tx ,tq , integratorCache,x , q,prevStepVal)
    savetime += savetimeincrement #next savetime 
  else#end if save
    for k = 1:T  
      
        elapsed = simt - tx[k];integrateState(Val(O),x[k],elapsed);tx[k] = simt #in case this point did not get updated.  
        elapsedq = simt - tq[k];integrateState(Val(O-1),q[k],elapsedq);tq[k]=simt        
      
      assignXPrevStepVals(Val(O),prevStepVal,x,k)
    end
  end
  prevStepTime=simt =#
     push!(savedVars[index],x[index][0])
    push!(savedTimes[index],simt)
  #  push!(savedVarsQ[index],q[index][0])



   #=  for i=1:T
      push!(savedVars[i],x[i][0])
      push!(savedTimes[i],simt)
      push!(savedVarsQ[i],q[i][0])
    end =#




end#end while

#= for i=1:T# throw away empty points
  resize!(savedVars[i],saveVarsHelper[1])
end
resize!(savedTimes,saveVarsHelper[1]) =#


createSol(Val(T),Val(O),savedTimes,savedVars, "qss$O",string(odep.prname),absQ,totalSteps,0,0,numSteps,ft)

end#end integrate

