
function initIntegrator(::Val{O},::Val{T},::Val{CS},odep, commonQssData::CommonQSS_Data{Z}, f::PRFUN) where {O,T,Z,CS,PRFUN}
  VERBOSE=commonQssData.verbose
  if VERBOSE println("integrating...") end


   


  ft = commonQssData.finalTime
  initTime = commonQssData.initialTime
  relQ = commonQssData.relQ
  absQ = commonQssData.absQ
  relZ = commonQssData.relZ
  absZ = commonQssData.absZ
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
 
 # discrete prob data
     clF=odep.closureFuncs[1]
  zc_SimpleJac = odep.ZCjac
  HZ = odep.HZ
  HD = odep.HD
  SZ = odep.SZ
  evDep = odep.eventDependencies

  #********************************helper values*******************************  
  oldsignValue = MMatrix{Z,2}(zeros(Z * 2))  #usedto track if zc changed sign; each zc has a value and a sign 
  numStateSteps = Vector{Int}(undef, T)
  numInputSteps = Vector{Int}(undef, T)

  # call once (for performance)
 
  f(1, -1, -1, q, d, t, taylorOpsCache,clF)

  #######################################compute initial values##################################################
  n = 1
   for k = 1:O # compute initial derivatives for x and q (similar to a recursive way )
    n = n * k
    for i = 1:T
      q[i].coeffs[k] = x[i].coeffs[k]
    end # q computed from x and it is going to be used in the next x
    for i = 1:T
      clearCache(taylorOpsCache, Val(CS), Val(O)); f(i, -1, -1, q, d, t, taylorOpsCache,clF)    
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
  t[0] = initTime



  for i = 1:Z
    clearCache(taylorOpsCache, Val(CS), Val(O))
    f(-1, i, -1, x, d, t, taylorOpsCache,clF)
    oldsignValue[i, 2] = taylorOpsCache[1][0] #value
    oldsignValue[i, 1] = sign(taylorOpsCache[1][0]) #sign modify 
    computeNextEventTime(Val(O), i, taylorOpsCache[1], oldsignValue, initTime, nextEventTime, quantum, absZ, relZ)
  end
  return VERBOSE,ft,initTime,relQ,absQ,relZ,absZ,maxErr,maxiters,quantum,nextStateTime,nextEventTime,nextInputTime,tx,tq,x,q,t,
        savedVars,savedTimes,taylorOpsCache,d,numStateSteps,numInputSteps,oldsignValue,clF,zc_SimpleJac,HZ,HD,SZ,evDep 
  
end#end init