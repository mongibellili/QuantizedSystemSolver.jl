#using TimerOutputs
function integrate(Al::QSSAlgorithm{:qss,O},CommonqssData::CommonQSS_data{Z}, odep::NLODEProblem{PRTYPE,T,Z,Y,CS},f::Function,jac::Function,SD::Function) where {PRTYPE,O,T,Z,Y,CS}

  ft = CommonqssData.finalTime;initTime = CommonqssData.initialTime;relQ = CommonqssData.dQrel;absQ = CommonqssData.dQmin;maxErr=CommonqssData.maxErr;;maxStepsAllowed=CommonqssData.maxStepsAllowed

  #savetimeincrement=CommonqssData.savetimeincrement;savetime = savetimeincrement
  quantum = CommonqssData.quantum;nextStateTime = CommonqssData.nextStateTime;nextEventTime = CommonqssData.nextEventTime;nextInputTime = CommonqssData.nextInputTime
  tx = CommonqssData.tx;tq = CommonqssData.tq;x = CommonqssData.x;q = CommonqssData.q;t=CommonqssData.t
  savedVars=CommonqssData.savedVars;
  savedTimes=CommonqssData.savedTimes;integratorCache=CommonqssData.integratorCache;taylorOpsCache=CommonqssData.taylorOpsCache;#cacheSize=odep.cacheSize
  #prevStepVal = specialLiqssData.prevStepVal
  #*********************************problem info*****************************************
  d = odep.discreteVars
  zc_SimpleJac=odep.ZCjac
  HZ=odep.HZ
  HD=odep.HD
  SZ=odep.SZ
  evDep = odep.eventDependencies

  if DEBUG @show HD,HZ,SZ,zc_SimpleJac,d end
  if DEBUG @show evDep end
  if DEBUG @show f end
  
  #********************************helper values*******************************  

  oldsignValue = MMatrix{Z,2}(zeros(Z*2))  #usedto track if zc changed sign; each zc has a value and a sign 
  numSteps = Vector{Int}(undef, T)
#######################################compute initial values##################################################
n=1
for k = 1:O # compute initial derivatives for x and q (similar to a recursive way )
  n=n*k
   for i = 1:T q[i].coeffs[k] = x[i].coeffs[k] end # q computed from x and it is going to be used in the next x
   for i = 1:T
      clearCache(taylorOpsCache,Val(CS),Val(O));f(i,-1,-1,q,d, t ,taylorOpsCache)
      ndifferentiate!(integratorCache,taylorOpsCache[1] , k - 1)
      x[i].coeffs[k+1] = (integratorCache.coeffs[1]) / n # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
    end
end

for i = 1:T
  numSteps[i]=0
   push!(savedVars[i],x[i][0])
   push!(savedTimes[i],0.0)
  quantum[i] = relQ * abs(x[i].coeffs[1]) ;quantum[i]=quantum[i] < absQ ? absQ : quantum[i];quantum[i]=quantum[i] > maxErr ? maxErr : quantum[i] 
  computeNextTime(Val(O), i, initTime, nextStateTime, x, quantum)
  #updateQ(Val(O),i,x,q,quantum,exacteA,cacheA,dxaux,qaux,tx,tq,initTime,ft,nextStateTime) 
end

for i=1:Z
  clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,i,-1,x,d,t,taylorOpsCache)                     
  oldsignValue[i,2]=taylorOpsCache[1][0] #value
  oldsignValue[i,1]=sign(taylorOpsCache[1][0]) #sign modify 
  computeNextEventTime(Val(O),i,taylorOpsCache[1],oldsignValue,initTime,  nextEventTime, quantum,absQ)
end

###################################################################################################################################################################
####################################################################################################################################################################
#---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
###################################################################################################################################################################
####################################################################################################################################################################
simt = initTime ;totalSteps=0;prevStepTime=initTime;modifiedIndex=0;statestep=0; countEvents=0;#needSaveEvent=false
  
while simt < ft && totalSteps < maxStepsAllowed
  if totalSteps==maxStepsAllowed-1
    @warn("The algorithm qss$O is taking too long to converge. The simulation will be stopped. Consider using a different algorithm!")
  end
  sch = updateScheduler(Val(T),nextStateTime,nextEventTime, nextInputTime)
  simt = sch[2];index = sch[1];stepType=sch[3]
   if  simt>ft  
    #saveLast!(Val(T),Val(O),savedVars, savedTimes,saveVarsHelper,ft,prevStepTime, x)
    break   ###################################################break##########################################
  end
  totalSteps+=1
  t[0]=simt
  DEBUG_time=DEBUG  && 0.0<=simt<=ft
  ##########################################state######################################## 
  if stepType == :ST_STATE
    statestep+=1
    numSteps[index]+=1;
   
    elapsed = simt - tx[index];integrateState(Val(O),x[index],elapsed);tx[index] = simt 
    quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index]   
    if abs(x[index].coeffs[2])>1e7 quantum[index]=10*quantum[index] end
    for k = 1:O q[index].coeffs[k] = x[index].coeffs[k] end; tq[index] = simt    
    computeNextTime(Val(O), index, simt, nextStateTime, x, quantum) #
    for j in (SD(index))
        elapsedx = simt - tx[j];if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt end
        elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],elapsedq);tq[j] = simt  end#q needs to be updated here for recomputeNext        
        for b in (jac(j)  )    
          elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
         end
        clearCache(taylorOpsCache,Val(CS),Val(O));f(j,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1])  
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
    end#end for SD

    for j in (SZ[index])
      for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
          elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
      end            
      clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,j,-1,q,d,t,taylorOpsCache)   # run ZCF--------      
      computeNextEventTime(Val(O),j,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum,absQ)
  end#end for SZ

   

    ##################################input########################################
  elseif stepType == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself    
  #################################################################event########################################
  else
      if DEBUG println("at start of event simt=$simt index=$index") end
      for b in zc_SimpleJac[index] # elapsed update all other vars that this zc depends upon.
          elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
      end    
      clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,index,-1,q,d,t,taylorOpsCache)    # run ZCF-------- 
      if DEBUG  @show oldsignValue[index,2],taylorOpsCache[1][0]  end
      
   #=    if oldsignValue[index,2]*taylorOpsCache[1][0]>=0 && abs(taylorOpsCache[1][0])>1e-9*absQ # if both have same sign and zcf is not very small
        computeNextEventTime(Val(O),index,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum,absQ) 
        if DEBUG  println("wrong estimation of event at $simt") end
        continue
      end
      # needSaveEvent=true
      countEvents+=1
      if taylorOpsCache[1][0]>oldsignValue[index,2] # rise
        modifiedIndex=2*index-1 
      elseif taylorOpsCache[1][0]<oldsignValue[index,2] # drop
        modifiedIndex=2*index
      else # == ( zcf==oldZCF)
        if DEBUG  println("ZCF==oldZCF at $simt") end
        continue
      end
      oldsignValue[index,2]=taylorOpsCache[1][0]
      oldsignValue[index,1]=sign(taylorOpsCache[1][0]) =#



      if oldsignValue[index,2]*taylorOpsCache[1][0]>=0  
        if abs(taylorOpsCache[1][0])>1e-9*absQ # if both have same sign and zcf is not very small: zc=1e-9*absQ is allowed as an event
                computeNextEventTime(Val(O),index,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum,absQ) 
              #  @show index,countEvents
              #  @show oldsignValue[index,2],taylorOpsCache[1][0]
              if DEBUG_time  println("wrong estimation of event at $simt") end
                continue
        end
      end
      if abs(oldsignValue[index,2]) <=1e-9*absQ  #earlier zc=1e-9*absQ was considered event , so now it should be prevented from passing
        nextEventTime[index]=Inf # at this instant next zc will be triggered now, and this will lead to infinite events, so cannot computenextevent here
        continue
      end
      # needSaveEvent=true
      
      if taylorOpsCache[1][0]>oldsignValue[index,2] #scheduled rise
        modifiedIndex=2*index-1 
      elseif taylorOpsCache[1][0]<oldsignValue[index,2] #scheduled drop
        modifiedIndex=2*index
      else # == ( zcf==oldZCF)
        if DEBUG_time  println("this should never be reached ZCF==oldZCF at $simt cuz small old not allowed and large zc not allowed!!") end
        computeNextEventTime(Val(O),index,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum,absQ) 
        continue
      end
      countEvents+=1
      oldsignValue[index,2]=taylorOpsCache[1][0]
      oldsignValue[index,1]=sign(taylorOpsCache[1][0])
      
     #=  if DEBUG 
        @show modifiedIndex,x 
      @show q
      end =#
          
      for b in evDep[modifiedIndex].evContRHS
          elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
      end
      f(-1,-1,modifiedIndex,x,d,t,taylorOpsCache)# execute event--------
      for i in evDep[modifiedIndex].evCont
        #------------event influences a Continete var: x already updated in event...here update quantum and q and computenext
            quantum[i] = relQ * abs(x[i].coeffs[1]) ;quantum[i]=quantum[i] < absQ ? absQ : quantum[i];quantum[i]=quantum[i] > maxErr ? maxErr : quantum[i] 
            for k = 0:O-1  q[i][k]=x[i][k] end;tx[i] = simt;tq[i] = simt # for liqss updateQ?
          #   firstguess=updateQ(Val(O),i,x,q,quantum,exacteA,cacheA,dxaux,qaux,tx,tq,simt,ft,nextStateTime) ;tq[i] = simt   
            computeNextTime(Val(O), i, simt, nextStateTime, x, quantum) 
      end
     # nextEventTime[index]=Inf   #investigate more 
      computeNextEventTime(Val(O),index,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum,absQ) # it could happen zcf=0.0 then infinite event

      for j in (HD[modifiedIndex]) # care about dependency to this event only     
          elapsedx = simt - tx[j];if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt;#= @show j,x[j] =# end
          elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],elapsedq);tq[j] = simt;#= @show q[j] =#  end#q needs to be updated here for recomputeNext                 
          for b = 1:T # elapsed update all other vars that this derj depends upon.
            if b in jac(j)   
              elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt;#= @show q[b] =# end
            end
          end
          clearCache(taylorOpsCache,Val(CS),Val(O));f(j,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1])
          reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
          # @show simt,j,x
      end
      for j in (HZ[modifiedIndex])
                for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
                    elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
                end            
              clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,j,-1,q,d,t,taylorOpsCache)  # run ZCF-------- 
              if VERBOSE @show j,oldsignValue[j,2],taylorOpsCache[1][0] end     
              computeNextEventTime(Val(O),j,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum,absQ)  
    
             
     end
     if DEBUG
      println("x at end of event simt=$simt x=$x") 
      println("q at end of event simt=$simt q=$q")
      @show countEvents,totalSteps,statestep
      @show nextStateTime,quantum
     end

  end#end state/input/event
  if stepType != :ST_EVENT
      push!(savedVars[index],x[index][0])
      push!(savedTimes[index],simt)

     #=  for i =1:T 
        push!(savedVars[i],x[i][0])
        
        push!(savedTimes[i],simt)
      end =#
  else
    #if needSaveEvent
    for j in (HD[modifiedIndex])
      push!(savedVars[j],x[j][0])
      push!(savedTimes[j],simt)
    end
    # end
  end

  prevStepTime=simt

end#end while
 
 #@show countEvents,totalSteps


createSol(Val(T),Val(O),savedTimes,savedVars, "qss$O",string(odep.prname),absQ,totalSteps,0,countEvents,numSteps,ft)
end#end integrate