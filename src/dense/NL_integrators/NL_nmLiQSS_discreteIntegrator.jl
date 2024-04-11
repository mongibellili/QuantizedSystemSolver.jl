#using TimerOutputs
function integrate(Al::QSSAlgorithm{:nmliqss,O},CommonqssData::CommonQSS_data{Z},liqssdata::LiQSS_data{O,false},specialLiqssData::SpecialLiqssQSS_data, odep::NLODEProblem{PRTYPE,T,Z,D,CS},f::Function,jac::Function,SD::Function,exactA::Function) where {PRTYPE,O,T,Z,D,CS}
  cacheA=specialLiqssData.cacheA
  ft = CommonqssData.finalTime;initTime = CommonqssData.initialTime;relQ = CommonqssData.dQrel;absQ = CommonqssData.dQmin;maxErr=CommonqssData.maxErr;

  savetimeincrement=CommonqssData.savetimeincrement;savetime = savetimeincrement
  quantum = CommonqssData.quantum;nextStateTime = CommonqssData.nextStateTime;nextEventTime = CommonqssData.nextEventTime;nextInputTime = CommonqssData.nextInputTime
  tx = CommonqssData.tx;tq = CommonqssData.tq;x = CommonqssData.x;q = CommonqssData.q;t=CommonqssData.t
  savedVars=CommonqssData.savedVars;
  savedTimes=CommonqssData.savedTimes;integratorCache=CommonqssData.integratorCache;taylorOpsCache=CommonqssData.taylorOpsCache;#cacheSize=odep.cacheSize
  #prevStepVal = specialLiqssData.prevStepVal
  #*********************************problem info*****************************************
  d = CommonqssData.d
  

  zc_SimpleJac=odep.ZCjac

  HZ=odep.HZ
  HD=odep.HD
  SZ=odep.SZ
 
  evDep = odep.eventDependencies

  if DEBUG2 @show HD,HZ end
 

  qaux=liqssdata.qaux;dxaux=liqssdata.dxaux#= olddx=liqssdata.olddx; ; olddxSpec=liqssdata.olddxSpec =#

  savedDers = Vector{Vector{Float64}}(undef, T)
  # savedVarsQ = Vector{Vector{Float64}}(undef, T) 
 
#=  setprecision(BigFloat,80)
  pp=pointer(Vector{NTuple{2,BigFloat}}(undef, 7))
 respp = pointer(Vector{BigFloat}(undef, 6))
 acceptedi=Vector{Vector{BigFloat}}(undef,4*O-1);acceptedj=Vector{Vector{BigFloat}}(undef,4*O-1); #inner vector always of size 2....interval...low and high...later optimize maybe
  cacheRootsi=Vector{BigFloat}(undef,8*O-4);
  cacheRootsj=Vector{BigFloat}(undef,8*O-4); =#
                            

  
  pp=pointer(Vector{NTuple{2,Float64}}(undef, 7))
  respp = pointer(Vector{Float64}(undef, 6))
  acceptedi=Vector{Vector{Float64}}(undef,4*O-1);acceptedj=Vector{Vector{Float64}}(undef,4*O-1);
  cacheRootsi=zeros(8*O-4) #vect of floats to hold roots for simul_analytic 
  cacheRootsj=zeros(8*O-4)



for i =1:4*O-1 #3
  acceptedi[i]=[0.0,0.0]#zeros(2)
  acceptedj[i]=[0.0,0.0]#zeros(2)
end
  exactA(q,d,cacheA,1,1,initTime+1e-9)
  trackSimul = Vector{Int}(undef, 1)
 # cacheRatio=zeros(5);cacheQ=zeros(5)
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
  savedDers[i]=Vector{Float64}()
  push!(savedVars[i],x[i][0])
   push!(savedDers[i],x[i][1])
   push!(savedTimes[i],0.0)
  
  quantum[i] = relQ * abs(x[i].coeffs[1]) ;quantum[i]=quantum[i] < absQ ? absQ : quantum[i];quantum[i]=quantum[i] > maxErr ? maxErr : quantum[i] 
  
  updateQ(Val(O),i,x,q,quantum,exactA,d,cacheA,dxaux,qaux,tx,tq,initTime+1e-9,ft,nextStateTime) #1e-9 exactAfunc contains 1/t
end
#= for i = 1:T
   clearCache(taylorOpsCache,Val(CS),Val(O));f(i,-1,-1,q,d,t,taylorOpsCache);
   computeDerivative(Val(O), x[i], taylorOpsCache[1])#0.0 used to be elapsed...even down below not neeeded anymore
  Liqss_reComputeNextTime(Val(O), i, initTime, nextStateTime, x, q, quantum)
  computeNextInputTime(Val(O), i, initTime, 0.1,taylorOpsCache[1] , nextInputTime, x,  quantum)#not complete, currently elapsed=0.1 is temp until fixed
 #prevStepVal[i]=x[i][0]#assignXPrevStepVals(Val(O),prevStepVal,x,i)
end =#
  

  #assignXPrevStepVals(Val(O),prevStepVal,x,i)
  


#@show nextStateTime,nextInputTime
for i=1:Z
  clearCache(taylorOpsCache,Val(CS),Val(O));  f(-1,i,-1,x,d,t,taylorOpsCache)                
  oldsignValue[i,2]=taylorOpsCache[1][0] #value
  oldsignValue[i,1]=sign(taylorOpsCache[1][0]) #sign modify 
  computeNextEventTime(Val(O),i,taylorOpsCache[1],oldsignValue,initTime,  nextEventTime, quantum,absQ)
end

###################################################################################################################################################################
####################################################################################################################################################################
#---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
###################################################################################################################################################################
####################################################################################################################################################################
simt = initTime ;totalSteps=0;prevStepTime=initTime;modifiedIndex=0; countEvents=0;inputstep=0;statestep=0;simulStepCount=0
ft<savetime && error("ft<savetime")
while simt< ft && totalSteps < 50000000
  
  sch = updateScheduler(Val(T),nextStateTime,nextEventTime, nextInputTime)
  simt = sch[2];index = sch[1];stepType=sch[3]
 # @timeit "saveLast" 
   if  simt>ft  
    #saveLast!(Val(T),Val(O),savedVars, savedTimes,saveVarsHelper,ft,prevStepTime, x)
    break   ###################################################break##########################################
  end
  totalSteps+=1
 
  t[0]=simt

  DEBUG_time=DEBUG  && 0.0002<=simt<=0.00022
  ##########################################state######################################## 
  if stepType == :ST_STATE
    statestep+=1
  
 
    xitemp=x[index][0]
    numSteps[index]+=1;
    
    elapsed = simt - tx[index];integrateState(Val(O),x[index],elapsed);tx[index] = simt ; 
    dirI=x[index][0]-xitemp
    if abs(dirI)>3*quantum[index] x[index][0]= 2*quantum[index] *sign(dirI) end # this is a rare case where dxi gets changed a lot by an event
   
    quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index]   
   
    if abs(x[index].coeffs[2])>1e9 quantum[index]=10*quantum[index] end  # i added this for the case a function is climbing (up/down) fast     

    
   
    for b in (jac(index)  )    # update Qb : to be used to calculate exacte Aindexb
      elapsedq = simt - tq[b] ;
      if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
    end
    firstguess=updateQ(Val(O),index,x,q,quantum,exactA,d,cacheA,dxaux,qaux,tx,tq,simt,ft,nextStateTime) ;tq[index] = simt   
    #----------------------------------------------------check dependecy cycles---------------------------------------------  
   
    trackSimul[1]=0 
    #= for i =1:5
      cacheRatio[i]=0.0; cacheQ[i]=0.0; 
    end       =#       
    for j in SD(index)
      for b in (jac(j)  )    # update Qb: to be used to calculate exacte Ajb
        elapsedq = simt - tq[b] ;
        if elapsedq>0  integrateState(Val(O-1),q[b],elapsedq); tq[b]=simt  end
      end
      cacheA[1]=0.0; exactA(q,d,cacheA,index,j,simt);aij=cacheA[1]# can be passed to simul so that i dont call exactfunc again
      cacheA[1]=0.0;exactA(q,d,cacheA,j,index,simt);aji=cacheA[1]
     
      
    
      if j!=index && aij*aji!=0.0
          #prvStepValj= savedVars[j][end]#getPrevStepVal(prevStepVal,j) 
         #=  if simt==0.0005977635736228422
              @show j,aij,aji
          end =#
         #=  for i =1:3
            cacherealPosi[i][1]=0.0; cacherealPosi[i][2]=0.0
            cacherealPosj[i][1]=0.0; cacherealPosj[i][2]=0.0
          end  =#
         # @show aij,aji
          if nmisCycle_and_simulUpdate(cacheRootsi,cacheRootsj,acceptedi,acceptedj,aij,aji,respp,pp,trackSimul,Val(O),index,j,dirI,firstguess,x,q,quantum,exactA,d,cacheA,dxaux,qaux,tx,tq,simt,ft)
            simulStepCount+=1


         
           clearCache(taylorOpsCache,Val(CS),Val(O));f(index,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[index], taylorOpsCache[1])
          #  clearCache(taylorOpsCache,Val(CS),Val(O));f(j,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1])
          # Liqss_reComputeNextTime(Val(O), index, simt, nextStateTime, x, q, quantum)
          #  Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
         

            for k in SD(j)  #j influences k
                if k!=index && k!=j
                    elapsedx = simt - tx[k]; x[k].coeffs[1] = x[k](elapsedx);tx[k] = simt
                    elapsedq = simt - tq[k]; if elapsedq > 0  ;integrateState(Val(O-1),q[k],elapsedq); tq[k] = simt end
                    for b in (jac(k)  )   
                      elapsedq = simt - tq[b]
                      if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
                    end  
                    clearCache(taylorOpsCache,Val(CS),Val(O));f(k,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[k], taylorOpsCache[1])
                   
                    Liqss_reComputeNextTime(Val(O), k, simt, nextStateTime, x, q, quantum)
                end#end if k!=0
            end#end for k depend on j     
            
            for k in (SZ[j]) # qj changed, so zcf should be checked
              for b in zc_SimpleJac[k] # elapsed update all other vars that this derj depends upon.
                  elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
                  #elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
              end            
              clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,k,-1,q,d,t,taylorOpsCache)   # run ZCF--------      
              computeNextEventTime(Val(O),k,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum,absQ)
            end#end for SZ

          end#end ifcycle check


      # end #end if allow one simulupdate
      end#end if j!=0
    end#end FOR_cycle check
        

  if trackSimul[1]!=0  #qi changed after throw
   # clearCache(taylorOpsCache,Val(CS),Val(O));f(index,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[index], taylorOpsCache[1])
    Liqss_reComputeNextTime(Val(O), index, simt, nextStateTime, x, q, quantum)
  end

  #-------------------------------------------------------------------------------------
  #---------------------------------normal liqss: proceed--------------------------------
  #-------------------------------------------------------------------------------------
 #=  if simt==0.0007352244633292483
    println("before normal dependency")
    @show x,q
    @show tx,tq
  end =#
    for c in SD(index)   #index influences c  
         
      elapsedx = simt - tx[c] ;
      if elapsedx>0 
       
        x[c].coeffs[1] = x[c](elapsedx);
       # if 0.0003>simt>0.00029 @show index, c,simt, tx[c], x[c].coeffs[1] end
        tx[c] = simt
        
       end # 
      elapsedq = simt - tq[c];if elapsedq > 0 ;integrateState(Val(O-1),q[c],elapsedq);tq[c] = simt    end   # c never been visited 
      #= for b in (jac(c)  )    # update other influences
          elapsedq = simt - tq[b] ;if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt  end
      end  =#
      clearCache(taylorOpsCache,Val(CS),Val(O)); f(c,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[c], taylorOpsCache[1])
      
      Liqss_reComputeNextTime(Val(O), c, simt, nextStateTime, x, q, quantum)
    end#end for SD
   #=  if 0.0003>simt>0.00029
      println("---after normal SD---------------", q,x,totalSteps)
    end =#
    for j in (SZ[index])
      for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
          elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
          #elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
      end            
      clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,j,-1,q,d,t,taylorOpsCache)   # run ZCF--------      
      computeNextEventTime(Val(O),j,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum,absQ)
  end#end for SZ
 
 #=  if DEBUG_time
     println("at simt=$simt x end of state step  = $x") 
     println("-------------q end of state step  = $q") 
     
    end =#

    if DEBUG_time && (index==3 || index==4)
      println("========end state=======")
      @show index,simt
      @show x,q
      @show nextStateTime,quantum
       end

    ##################################input########################################
  elseif stepType == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself    
     inputstep+=1
    #=   if VERBOSE println("nmliqss discreteintgrator under input, index= $index, totalsteps= $totalSteps")  end
    elapsed = simt - tx[index];integrateState(Val(O),x[index],elapsed);tx[index] = simt 
    quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index]   
    for k = 1:O q[index].coeffs[k] = x[index].coeffs[k] end; tq[index] = simt 
      for b in jac(index) 
        elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
      end
    clearCache(taylorOpsCache,Val(CS),Val(O));f(index,-1,-1,q,d,t,taylorOpsCache)
    computeNextInputTime(Val(O), index, simt, elapsed,taylorOpsCache[1] , nextInputTime, x,  quantum)
    computeDerivative(Val(O), x[index], taylorOpsCache[1])

    for j in(SD(index))  
      elapsedx = simt - tx[j];
      if elapsedx > 0 
        x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt 
       # quantum[j] = relQ * abs(x[j].coeffs[1]) ;quantum[j]=quantum[j] < absQ ? absQ : quantum[j];quantum[j]=quantum[j] > maxErr ? maxErr : quantum[j]   
      end
      elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],elapsedq);tq[j] = simt  end#q needs to be updated here for recomputeNext                 
      # elapsed update all other vars that this derj depends upon.
        for b in jac(j) 
          elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
        end
      
        clearCache(taylorOpsCache,Val(CS),Val(O));f(j,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
    end#end for
    for j in (SZ[index])
      
        for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
             
            elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
           # elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
       
        end              
       #=  clearCache(taylorOpsCache,Val(CS),Val(O))#normally and later i should update x,q (integrate q=q+e derQ  for higher orders)
        computeNextEventTime(Val(O),j,zcf[j](x,d,t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum)#,maxIterer)  =#
        clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,j,-1,x,d,t,taylorOpsCache)        
         computeNextEventTime(Val(O),j,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum)
     
    end =#
  #################################################################event########################################
  
else
  
        
        if DEBUG_time
           println("x at start of event simt=$simt index=$index") 
          
         # println("-------------q start of event  = $q") 
         #@show index,d
        end

       

          for b in zc_SimpleJac[index] # elapsed update all other vars that this zc depends upon.
             elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
          end    
         # modifiedIndex=0#first we have a zc happened which corresponds to nexteventtime and index (one of zc) but we want also the sign in O to know ev+ or ev- 
        
          clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,index,-1,q,d,t,taylorOpsCache)    # run ZCF-------- 
          if DEBUG_time
             @show index,oldsignValue[index,2],taylorOpsCache[1][0] 
             end
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
          
          #= if DEBUG_time 
            @show modifiedIndex,x
            @show q
          end =#
              
          for b in evDep[modifiedIndex].evContRHS
              elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
          end
          f(-1,-1,modifiedIndex,x,d,t,taylorOpsCache)# execute event----------------no need to clear cache; events touch vectors directly
       
          if VERBOSE 
            @show d
          end

          for i in evDep[modifiedIndex].evCont
            #------------event influences a Continete var: ...here update quantum and q and computenext
          #=   if DEBUG_time
             
              println("q begin evcont simt=$simt q2=$(q[2])")
           
             @show nextStateTime[2],quantum[2]
             end =#
                quantum[i] = relQ * abs(x[i].coeffs[1]) ;quantum[i]=quantum[i] < absQ ? absQ : quantum[i];quantum[i]=quantum[i] > maxErr ? maxErr : quantum[i] 
               # q[i][0]=x[i][0]; # for liqss updateQ?

              #=  if abs(x[i].coeffs[2])>1e9
               # @show quantum[i]
                quantum[i]=10*quantum[i]
              end =#

               firstguess=updateQ(Val(O),i,x,q,quantum,exactA,d,cacheA,dxaux,qaux,tx,tq,simt,ft,nextStateTime)   
              #  computeNextTime(Val(O), i, simt, nextStateTime, x, quantum) 
              tx[i] = simt;tq[i] = simt
              #Liqss_reComputeNextTime(Val(O), i, simt, nextStateTime, x, q, quantum)
          #=     if DEBUG_time
                println("x end evcont simt=$simt x2=$(x[2])") 
                println("q end evcont simt=$simt q2=$(q[2])")
               @show countEvents,totalSteps,statestep
               @show nextStateTime[2],quantum[2]
               end =#
           
          end
         #=  if VERBOSE 
            @show x,q,nextStateTime
          end =#
         # nextEventTime[index]=Inf   #investigate more 
          computeNextEventTime(Val(O),index,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum,absQ) #update zcf before thiscatch in qss quantizer to avoid infinite events
         
          for j in (HD[modifiedIndex]) # care about dependency to this event only
                 
              elapsedx = simt - tx[j];if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt;#= @show j,x[j] =# end
              elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],elapsedq);tq[j] = simt;#= @show q[j] =#  end#q needs to be updated here for recomputeNext                 
              for b = 1:T # elapsed update all other vars that this derj depends upon.
                if b in jac(j)   
                  elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt;#= @show q[b] =# end
                end
              end
          

              clearCache(taylorOpsCache,Val(CS),Val(O));f(j,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1])
 
              #firstguess=updateQ(Val(O),j,x,q,quantum,exactA,d,cacheA,dxaux,qaux,tx,tq,simt,ft,nextStateTime)  
              Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
        
           
          end
          for j in (HZ[modifiedIndex])
                
                    for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
                          
                        #elapsedx = simt - tx[b];if elapsedx>0 integrateState(Val(O),x[b],elapsedx);tx[b]=simt end
                       elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
                    
                    end            
                  #= clearCache(taylorOpsCache,Val(CS),Val(O)) #normally and later i should update q (integrate q=q+e derQ  for higher orders)          
                  computeNextEventTime(Val(O),j,zcf[j](x,d,t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum)#,maxIterer) =#
                  
                  clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,j,-1,q,d,t,taylorOpsCache)  # run ZCF-------- 
                  if VERBOSE @show j,oldsignValue[j,2],taylorOpsCache[1][0] end     
                 computeNextEventTime(Val(O),j,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum,absQ)
                
              # if 0.022>simt > 0.018  println("$index $j at simt=$simt nexteventtime from HZ= ",nextEventTime)   end   
          end
       
        
          if DEBUG_time
             println("x at end of event simt=$simt x2=$(x)") 
             println("q at end of event simt=$simt q2=$(q)")
            @show countEvents,totalSteps,statestep
            @show nextStateTime,quantum
            end
          #=   if 9.236846540089048e-5<=simt<9.237519926276279e-5 
              println("-------------end of event------------")
              @show simt,nextStateTime[5]
              @show x[5],q[5],quantum[5]
            
              end =#
        
  end#end state/input/event
  #for i=1:T
 # if simt > savetime  #limit saveat
    savetime = simt+savetimeincrement
      if stepType != :ST_EVENT
      
          push!(savedVars[index],x[index][0])
          #push!(savedDers[index],x[index][0])
          push!(savedTimes[index],simt)
    
         #=  for i =1:T 
            push!(savedVars[i],x[i][0])
            push!(savedDers[i],x[i][1])
            push!(savedTimes[i],simt)
          end =#
      else
        #countEvents+=1
        # if 1e-3<simt<1e-2 || 1e-5<simt<1e-4
        for j in (HD[modifiedIndex])
          
          push!(savedVars[j],x[j][0])
          #push!(savedDers[j],x[j][1])
          push!(savedTimes[j],simt)
          
        end
      # end
      end
      
  #  end
    #push!(savedVarsQ[i],q[i][0])
    prevStepTime=simt
# end
end#end while
 
#@show countEvents,inputstep,statestep,simulStepCount
#@show savedVars
#createSol(Val(T),Val(O),savedTimes,savedVars, "qss$O",string(nameof(f)),absQ,totalSteps,0)#0 I track simulSteps 
#createSol(Val(T),Val(O),savedTimes,savedVars, "nmLiqss$O",string(odep.prname),absQ,totalSteps,simulStepCount,countEvents,numSteps,ft)
createSol(Val(T),Val(O),savedTimes,savedVars#= ,savedDers =#, "nmLiqss$O",string(odep.prname),absQ,totalSteps,simulStepCount,countEvents,numSteps,ft)
end#end integrate