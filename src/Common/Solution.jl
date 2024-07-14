struct Stats
  totalSteps::Int
  simulStepCount::Int
  evCount::Int
  numSteps ::Vector{Int}
end

"""LightSol{T,O}
A struct that holds the solution of a system of ODEs. It has the following fields:\n
    - size: The number of continuous variables T\n
    - order: The order of the algorithm O\n
    - savedTimes: A vector of vectors of Float64 that holds the times at which the continuous variables were saved\n
    - savedVars: A vector of vectors of Float64 that holds the values of the continuous variables at the times they were saved\n
    - algName: The name of the algorithm used to solve the system\n
    - sysName: The name of the system\n
    - absQ: The absolute tolerance used in the simulation\n
    - totalSteps: The total number of steps taken by the algorithm\n
    - simulStepCount: The number of simultaneous updates during the simulation\n
    - evCount: The number of events that occurred during the simulation\n
    - numSteps: A vector of Int that holds the number of steps taken by the algorithm for each continuous variable\n
    - ft: The final time of the simulation
"""
struct LightSol{T,O}<:Sol{T,O}
  size::Val{T}
  order::Val{O}
  savedTimes::Vector{Vector{Float64}}
  savedVars::Vector{Vector{Float64}}
  #savedVarsQ::Vector{Vector{Float64}}
  algName::String
  sysName::String
  absQ::Float64
  stats::Stats
  ft::Float64
end

@inline function createSol(::Val{T},::Val{O}, savedTimes:: Vector{Vector{Float64}},savedVars :: Vector{Vector{Float64}},solver::String,nameof_F::String,absQ::Float64,stats::Stats,ft::Float64)where {T,O}
 # println("light")
  sol=LightSol(Val(T),Val(O),savedTimes, savedVars,solver,nameof_F,absQ,stats,ft#= ,simulStepsVals,simulStepsDers,simulStepsTimes =#)
end
function getindex(s::Sol, i::Int64)
  if i==1
     return s.savedTimes
  elseif i==2
     return s.savedVars
  else
     error("sol has 2 attributes: time and states")
  end
end

@inline function evaluateSol(sol::Sol{T,O},index::Int,t::Float64)where {T,O}
  (t>sol.ft) && error("given point is outside the solution range! Verify where you want to evaluate the solution")

  x=sol[2][index][end] 
  #integratorCache=Taylor0(zeros(O+1),O)
  for i=2:length(sol[1][index])#savedTimes after the init time...init time is at index i=1
      if sol[1][index][i]>t # i-1 is closest lower point
        f1=sol[2][index][i-1];f2=sol[2][index][i];t1=sol[1][index][i-1] ;t2=sol[1][index][i]# find x=f(t)=at+b...linear interpolation
        a=(f2-f1)/(t2-t1)
        b=(f1*t2-f2*t1)/(t2-t1)
        x=a*t+b
       # println("1st case")
        return x#taylor evaluation after small elapsed with the point before (i-1)
      elseif sol[1][index][i]==t # i-1 is closest lower point
        x=sol[2][index][i]
      #  println("2nd case")
        return x
      end
  end
 # println("3rd case")
  return x #if var never changed then return init cond or if t>lastSavedTime for this var then return last value
end
function solInterpolated(sol::Sol{T,O},step::Float64)where {T,O}
  interpTimes=Float64[]
  allInterpTimes=Vector{Vector{Float64}}(undef, T)
  t=0.0  #later can change to init_time which could be diff than zero
  push!(interpTimes,t)
  while t+step<sol.ft
    t=t+step
    push!(interpTimes,t) 
  end
  push!(interpTimes,sol.ft)
  numInterpPoints=length(interpTimes)
  interpValues=nothing
  if sol isa LightSol
    interpValues=Vector{Vector{Float64}}(undef, T)
  end
  for index=1:T
          interpValues[index]=[]
          push!(interpValues[index],sol[2][index][1]) #1st element is the init cond (true value)
        for i=2:numInterpPoints-1
          push!(interpValues[index],evaluateSol(sol,index,interpTimes[i]))
        end
          push!(interpValues[index],sol[2][index][end]) #last pt @ft
          allInterpTimes[index]=interpTimes
  end
  createSol(Val(T),Val(O),allInterpTimes,interpValues,sol.algName,sol.sysName,sol.absQ,sol.stats,sol.ft#= ,sol.simulStepsVals,sol.simulStepsDers,sol.simulStepsVals =#)
end
(sol::Sol)(index::Int,t::Float64) = evaluateSol(sol,index,t)
(sol::Sol)(t::Float64;idxs=1::Int) = evaluateSol(sol,idxs,t)

function show(io::IO, a::Stats)
  println("The total simulation steps: $(a.totalSteps)")
  println("The simultaneous  steps: $(a.simulStepCount)")
  println("The number of events: $(a.evCount)")
  println("The simulation steps per variable: $(a.numSteps)")
  
end