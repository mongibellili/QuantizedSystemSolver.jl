"""
    updateScheduler(::Val{T}, nextStateTime::Vector{Float64}, nextEventTime::MVector{Z,Float64}, nextInputTime::Vector{Float64}) where {T, Z}

Updates the scheduler by finding the minimum times among state, event, and input transitions.

# Arguments
- `::Val{T}`: A type parameter indicating the number of state transitions.
- `nextStateTime::Vector{Float64}`: A vector of times for the next state transitions.
- `nextEventTime::MVector{Z,Float64}`: A mutable vector of times for the next event transitions.
- `nextInputTime::Vector{Float64}`: A vector of times for the next input transitions.

# Returns
- A tuple containing the minimum (state index,state time), (event index, event tim, or (input index, input time).
"""
function updateScheduler(::Val{T},nextStateTime::Vector{Float64},nextEventTime :: MVector{Z,Float64},nextInputTime :: Vector{Float64})where{T,Z}   #later MVect for nextInput
    minStateTime=Inf
    minState_index=0  # what if all nextstateTime= Inf ...especially at begining????? min_index stays 0!!!
    minEventTime=Inf
    minEvent_index=0
    minInputTime=Inf
    minInput_index=0

    returnedVar=() #  used to print something if something is bad
    for i=1:T
        if nextStateTime[i]<minStateTime
            minStateTime=nextStateTime[i]
            minState_index=i
        end
        if nextInputTime[i] < minInputTime
            minInputTime=nextInputTime[i]
            minInput_index=i
        end
    end
    for i=1:Z
        if nextEventTime[i] < minEventTime
            minEventTime=nextEventTime[i]
            minEvent_index=i
        end
    end
  
    if minEventTime<minStateTime
       if minInputTime<minEventTime
         returnedVar=(minInput_index,minInputTime,:ST_INPUT)
       else
          returnedVar=(minEvent_index,minEventTime,:ST_EVENT)     
       end
    else      
        if minInputTime<minStateTime
            returnedVar=(minInput_index,minInputTime,:ST_INPUT)
        else
             returnedVar=(minState_index,minStateTime,:ST_STATE)
        end
    end
    if returnedVar[1]==0
        returnedVar=(1,Inf,:ST_STATE)
        println("null step made state step")
    end
    return returnedVar 
end