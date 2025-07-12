"""
    integrate(alg::QSSAlgorithm{:qss,O}, commonQssData::CommonQSS_Data{0}, odep::ODEProblemData{F,PRTYPE,T,D,0,CS}, f::Function, jac::Function, SD::Function) where {F,PRTYPE,O,T,CS,D}

Integrates a nonlinear ordinary differential equation (ODE) problem using a Quantized State System (QSS) algorithm.

# Arguments
- `Al::QSSAlgorithm{:qss,O}`: The QSS algorithm to be used for integration.
- `commonQssData::CommonQSS_Data{0}`: Common data structure for QSS integration.
- `odep::ODEProblemData{F,PRTYPE,T,0,0,CS}`: The nonlinear ODE problem to be solved.
- `f::Function`: The function defining the ODE system.
- `jac::Function`: The Jacobian dependency function of the ODE system.
- `SD::Function`: The state derivative dependency function.


# Returns
- A solution


"""
function integrate(alg::QSSAlgorithm{:qss,O}, commonQssData::CommonQSS_Data{0}, odep::ODEProblemData{F,PRTYPE,T,D,0,CS}, f::Function, jac::Function, SD::Function) where {F,PRTYPE,O,T,CS,D}
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
  taylorOpsCache = commonQssData.taylorOpsCache
  d = commonQssData.d 
  clF=odep.closureFuncs[1]
  numStateSteps = Vector{Int}(undef, T)
  numInputSteps = Vector{Int}(undef, T)
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
  smallAdvance = ft / 1000
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
        f(i, -1, -1, q, d, t + smallAdvance, taylorOpsCache,clF) =#
       # computeNextInputTime(Val(O), i, initTime, smallAdvance, taylorOpsCache[1], nextInputTime, x, quantum)
        computeNextInputTime(Val(O), i, t, f,clF,d, taylorOpsCache, nextInputTime, x,q, quantum) 
        #= if nextInputTime[i] > initTime + 2 * smallAdvance
          nextInputTime[i] = initTime + 2 * smallAdvance
        end =#
      #end
    else
      computeNextTime(Val(O), i, initTime, nextStateTime, x, quantum)
    end
  end
  simt = initTime
  totalSteps = 0
   inpuStepCount=0
  ###################################################################################################################################################################
  ####################################################################################################################################################################
  #---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
  ###################################################################################################################################################################
  #################################################################################################################################################################### 
  
  while simt < ft && totalSteps < maxiters
    if totalSteps == maxiters - 1 @warn("The algorithm qss$O reached max iterations. The simulation will be stopped. Consider using a different algorithm!") end
    sch = updateScheduler(Val(T), nextStateTime, nextEventTime, nextInputTime)
    simt = sch[2]
    index = sch[1]
    totalSteps += 1
    t[0] = simt
    if simt > ft
      break
    end
    if sch[3] == :ST_STATE
      numStateSteps[index] += 1
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
      computeNextTime(Val(O), index, simt, nextStateTime, x, quantum)
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
        end
        for b in (jac(j))
          elapsedq = simt - tq[b]
          if elapsedq > 0
            integrateState(Val(O - 1), q[b], elapsedq)
            tq[b] = simt
          end
        end
        clearCache(taylorOpsCache, Val(CS), Val(O))
        f(j, q, d, t,taylorOpsCache,clF)
        computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end
    elseif sch[3] == :ST_INPUT
      numInputSteps[index] += 1
       inpuStepCount+=1
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
      #computeNextTime(Val(O), index, simt, nextStateTime, x, quantum) # 
      clearCache(taylorOpsCache, Val(CS), Val(O))
      f(index, q, d, t,taylorOpsCache,clF)
      computeDerivative(Val(O), x[index], taylorOpsCache[1])
      #reComputeNextTime(Val(O), index, simt, nextStateTime, x, q, quantum)
      #if nextInputTime[index] > simt + 2 * elapsed
        #= clearCache(taylorOpsCache, Val(CS), Val(O))
        f(index, q, t + 0.01*smallAdvance, d,taylorOpsCache,clF)
        computeNextInputTime(Val(O), index, simt, 0.01*smallAdvance, taylorOpsCache[1], nextInputTime, x, quantum) =#
        computeNextInputTime(Val(O), index, t, f,clF,d, taylorOpsCache, nextInputTime, x,q, quantum)
        #= if nextInputTime[index] > simt + 2 * elapsed
          nextInputTime[index] = simt + 2 * elapsed
        end =#
      #end
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
        f(j, q, d, t,taylorOpsCache,clF)
        computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end#end for
    end
    push!(savedVars[index], x[index][0])
    push!(savedTimes[index], simt)
  end
  stats=Stats(totalSteps,0,inpuStepCount,0,numStateSteps,numInputSteps)
  createSol(Val(T), Val(O), savedTimes, savedVars, toString(alg), string(odep.prname), absQ, stats, ft)
end
