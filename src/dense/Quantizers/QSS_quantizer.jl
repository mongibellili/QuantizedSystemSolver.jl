
######################################################################################################################################"
function computeNextTime(::Val{1}, i::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0}, quantum::Vector{Float64})#i can be absorbed
  absDeltaT=1e-12 # minimum deltaT to protect against der=Inf coming from sqrt(0) for example...similar to min ΔQ
    if (x[i].coeffs[2]) != 0
        tempTime=max(abs(quantum[i] /(x[i].coeffs[2])),absDeltaT)# i can avoid the use of max
        if tempTime!=absDeltaT #normal
            nextTime[i] = simt + tempTime#sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
        else#usual (quant/der) is very small
          x[i].coeffs[2]=sign(x[i].coeffs[2])*(abs(quantum[i])/absDeltaT)# adjust  derivative if it is too high
          nextTime[i] = simt + tempTime
          if DEBUG  println("qssQuantizer: smalldelta in compute next") end
        end
    else
      nextTime[i] = Inf
    end
    return nothing
end

function computeNextTime(::Val{2}, i::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0}, quantum::Vector{Float64})
    absDeltaT=1e-12 # minimum deltaT to protect against der=Inf coming from sqrt(0) for example...similar to min ΔQ
      if (x[i].coeffs[3]) != 0
          tempTime=max(sqrt(abs(quantum[i] / ((x[i].coeffs[3])))),absDeltaT)
          
          if tempTime!=absDeltaT #normal
              nextTime[i] = simt + tempTime#sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
          else#usual sqrt(quant/der) is very small
            x[i].coeffs[3]=sign(x[i].coeffs[3])*(abs(quantum[i])/(absDeltaT*absDeltaT))/2# adjust second derivative if it is too high
            nextTime[i] = simt + tempTime
            if DEBUG  println("qssQuantizer: smalldelta in compute next") end
          end
      else
        if (x[i].coeffs[2]) != 0
          #quantum[i]=2*quantum[i]
          tempTime=max(abs(quantum[i] /(x[i].coeffs[2])),absDeltaT)# i can avoid the use of max
          if tempTime!=absDeltaT #normal
              nextTime[i] = simt + tempTime#sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
          else#usual (quant/der) is very small
            x[i].coeffs[2]=sign(x[i].coeffs[2])*(abs(quantum[i])/absDeltaT)# adjust  derivative if it is too high
            nextTime[i] = simt + tempTime
            if DEBUG  println("qssQuantizer: smalldelta in compute next") end
          end
        else
          nextTime[i] = Inf
        end
      end
      return nothing
end
#= function computeNextTime(::Val{3}, i::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0}, quantum::Vector{Float64})
  absDeltaT=1e-12 # minimum deltaT to protect against der=Inf coming from sqrt(0) for example...similar to min ΔQ
    if (x[i].coeffs[4]) != 0
        tempTime=max(cbrt(abs(quantum[i] / ((x[i].coeffs[4])))),absDeltaT)   #6/6
        if tempTime!=absDeltaT #normal
            nextTime[i] = simt + tempTime#sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
        else#usual sqrt(quant/der) is very small
          x[i].coeffs[4]=sign(x[i].coeffs[4])*(abs(quantum[i])/(absDeltaT*absDeltaT*absDeltaT))/6# adjust third derivative if it is too high
          nextTime[i] = simt + tempTime
          if DEBUG  println("qssQuantizer: smalldelta in compute next") end
        end
    else
      nextTime[i] = Inf
    end
    return nothing
end =#
######################################################################################################################################"
function reComputeNextTime(::Val{1}, index::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64})
  absDeltaT=1e-12
  coef=@SVector [q[index].coeffs[1] - (x[index].coeffs[1]) - quantum[index], -x[index].coeffs[2]]
  time1 = simt + minPosRoot(coef, Val(1))
  coef=setindex(coef,q[index].coeffs[1] - (x[index].coeffs[1]) + quantum[index],1)
  time2 = simt + minPosRoot(coef, Val(1))
  timeTemp = time1 < time2 ? time1 : time2
  tempTime=max(timeTemp,absDeltaT) #guard against very small Δt 

  #if tempTime==absDeltaT #normal
  # if DEBUG  println("smalldelta in recompute next, simt= ",simt) end
 # end



  nextTime[index] = simt +tempTime
end

function reComputeNextTime(::Val{2}, index::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64})
  absDeltaT=1e-15

  if abs(q[index].coeffs[1] - (x[index].coeffs[1])) >= quantum[index] # this happened when var i and j s turns are now...var i depends on j, j is asked here for next time
    nextTime[index] = simt+1e-15
  else
    coef=@SVector [q[index].coeffs[1] - (x[index].coeffs[1]) - quantum[index], q[index].coeffs[2]-x[index].coeffs[2],-(x[index].coeffs[3])]#not *2 because i am solving c+bt+a/2*t^2
    time1 =  minPosRoot(coef, Val(2))
    coef=setindex(coef,q[index].coeffs[1] - (x[index].coeffs[1]) + quantum[index],1)
    time2 =  minPosRoot(coef, Val(2))
    timeTemp = time1 < time2 ? time1 : time2
    tempTime=max(timeTemp,absDeltaT)#guard against very small Δt 

  # if tempTime==absDeltaT #normal
    # x[index].coeffs[3]=sign(x[index].coeffs[3])*(abs(quantum[index])/(absDeltaT*absDeltaT))/2
    # if DEBUG  println("smalldelta in recompute next, simt= ",simt) end
  # end



    nextTime[index] = simt +tempTime

  end
 
  return nothing
end
  
#= function reComputeNextTime(::Val{3}, index::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64})
  #coef=@SVector [q[index].coeffs[1] - (x[index].coeffs[1]) - quantum[index], q[index].coeffs[2]-x[index].coeffs[2],(q[index].coeffs[3])-(x[index].coeffs[3]),-(x[index].coeffs[4])]
  #time1 = simt + minPosRoot(coef, Val(3))
  #pp=pointer(Vector{NTuple{2,Float64}}(undef, 5))
  a=-x[index][3];b=q[index][2]-x[index][2];c=q[index][1]-x[index][1];d=q[index][0]-x[index][0]-quantum[index]
  #GC.enable(false)
 #time1 = simt+smallestpositiverootintervalnewtonregulafalsi((a, b, c, d))
  time1 = simt+cubic5(a, b, c, d)
 # GC.enable(true)
  #coef=setindex(coef,q[index].coeffs[1] - (x[index].coeffs[1]) + quantum[index],1)
  d=q[index][0]-x[index][0]+quantum[index]
 # GC.enable(false)
  #time2 = simt+smallestpositiverootintervalnewtonregulafalsi((a, b, c, d))
  time2 = simt+cubic5(a, b, c, d)
 # GC.enable(true)
  #time2 = simt + minPosRoot(coef, Val(3))
  timeTemp = time1 < time2 ? time1 : time2
  absDeltaT=1e-12
  tempTime=max(timeTemp,absDeltaT)#guard against very small Δt 

  nextTime[index] = simt +tempTime

  return nothing
end =#
######################################################################################################################################"
#= function computeNextInputTime(::Val{1}, i::Int, simt::Float64,elapsed::Float64, tt::Taylor0 ,nextInputTime::Vector{Float64}, x::Vector{Taylor0}, quantum::Vector{Float64})
    df=0.0
    oldDerX=x[i].coeffs[2]
    newDerX=tt.coeffs[1] 
    if elapsed > 0
      df=(newDerX-oldDerX)/elapsed #df mimics second der
    else
      df= quantum[i]*1e6
    end     
    if df!=0.0
       nextInputTime[i]=simt+sqrt(abs(2*quantum[i] / df))
    else
      nextInputTime[i] = Inf
    end
    return nothing
end

  

function computeNextInputTime(::Val{2}, i::Int, simt::Float64,elapsed::Float64, tt::Taylor0 ,nextInputTime::Vector{Float64}, x::Vector{Taylor0}, quantum::Vector{Float64})
  ddf=0.0;df=0.0
  oldDerDerX=((x[i].coeffs[3])*2.0)
 # @show x
  newDerDerX=(tt.coeffs[2])# 1st der of tt cuz tt itself is derx=f
  #@show tt
  if elapsed > 0.0
      ddf=(newDerDerX-oldDerDerX)/(elapsed)
    #  println("df=new-old= ",df)
  else
      ddf= quantum[i]*1e6#*1e12
      if DEBUG  println("QSS quantizer elapsed=0!") end
  end       
  if ddf!=0.0
      nextInputTime[i]=simt+cbrt((abs(quantum[i]/df)))    #ddf mimics 3rd der 
     # println("usedcbrt")
  else #df=0->newddx=oldddx ->
    if DEBUG  println("QSS quantizer ddf=0 ") end
    oldDerX=((x[i].coeffs[2]))
    # @show x
     newDerX=(tt.coeffs[1])# 1st der of tt cuz tt itself is derx=f
    # @show tt
     if elapsed > 0.0
         df=(newDerX-oldDerX)/(elapsed) #df mimic second derivative
       #  println("df=new-old= ",df)
     else
         df= quantum[i]*1e6#*1e12
         if DEBUG  println("QSS quantizer elapsed=0!") end
     end       
     if df!=0.0
        
         nextInputTime[i]=simt+sqrt((abs(quantum[i]/df))) 
        
     else
      if x[i][1]!=0 #predicted second derivative is 0 & should not be used to determine nexttime. use 1st der
        #nextInputTime[i]=simt+(abs(2*quantum[i] / x[i][1]))  #*2 is not analytic: is just there to increase stepsize
        nextInputTime[i]=simt+sqrt(abs(1*quantum[i] / x[i][1]))  #I used the same formulae even with 1st der so that it is fair to other vars
     else
         nextInputTime[i] = Inf
     end
     end
  end
 # @show nextInputTime
    return nothing
end
function computeNextInputTime(::Val{3}, i::Int, simt::Float64,elapsed::Float64, tt::Taylor0 ,nextInputTime::Vector{Float64}, x::Vector{Taylor0}, quantum::Vector{Float64})
  df=0.0
  oldDerDerX=((x[i].coeffs[3])*2.0)
  #@show x
  newDerDerX=(tt.coeffs[2])# 1st der of tt cuz tt itself is derx=f
  #@show newDerDerX
  if elapsed > 0.0
      df=(newDerDerX-oldDerDerX)/(elapsed)
     # println("df=new-old= ",df)
  else
      df= quantum[i]*1e6#*1e12
  end       
  if df!=0.0
      nextInputTime[i]=simt+cbrt((abs(quantum[i]/df)))    #df mimics 3rd der 
     # println("usedcbrt")
  else
      if newDerDerX==0 && x[i][1]!=0 #predicted second derivative is 0 & should be not be used to determine nexttime. use 1st der
         #nextInputTime[i]=simt+(abs(2*quantum[i] / x[i][1]))  #*2 is not analytic: is just there to increase stepsize
         nextInputTime[i]=simt+sqrt(abs(1*quantum[i] / x[i][1]))  #I used the same formulae even with 1st der so that it is fair to other vars
      else
          nextInputTime[i] = Inf
      end
  end
    return nothing
end
 =#




function computeNextEventTime(::Val{O},j::Int,ZCFun::Taylor0,oldsignValue,simt,  nextEventTime, quantum::Vector{Float64},absQ::Float64) where {O}
 #=  if ZCFun[0]==0.0 && oldsignValue[j,1] !=0.0#abs(oldsignValue[j,1])>1e-15
    if DEBUG2 println("qss quantizer:zcf$j ZCF=0.0 (rare) immediate event at simt= $simt oldzcf value= $(oldsignValue[j,2])  ") end
    nextEventTime[j]=simt 
  else =#
  if (oldsignValue[j,1] != sign(ZCFun[0])) && abs(oldsignValue[j,2]) >1e-9*absQ #prevent double tapping: when zcf is leaving zero it should be considered an event
    if DEBUG  println("qss quantizer:zcf$j immediate event at simt= $simt oldzcf value= $(oldsignValue[j,2])  newZCF value= $(ZCFun[0])");  end
    nextEventTime[j]=simt 
  elseif oldsignValue[j,2] ==0.0 && ZCFun[0] ==0.0  # initial value 0 --> do not do anything
    nextEventTime[j]=Inf
  else # old and new ZCF both pos or both neg
   # coef=@SVector [ZCFun[0],ZCFun[1],ZCFun[2]]

    mpr=minPosRoot(ZCFun, Val(O)) 
    if mpr<1e-14 # prevent very close events
      mpr=1e-12
      #sl=ZCFun[0]+mpr*ZCFun[1]+mpr*mpr*ZCFun[2]/2
    end
    nextEventTime[j] =simt + mpr
      if DEBUG2
         println("qss quantizer:zcf$j at simt= $simt ****scheduled**** event at $(simt + mpr) oldzcf value= $(oldsignValue[j,2])  newZCF value= $(ZCFun[0])") 
        end
    oldsignValue[j,1]=sign(ZCFun[0])#update the values
    oldsignValue[j,2]=ZCFun[0]
  end
end