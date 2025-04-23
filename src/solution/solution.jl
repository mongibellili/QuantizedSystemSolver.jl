"""
    Stats

A struct that holds the statistics of a simulation. It has the following fields:
  - `totalSteps::Int`: The total number of simulation steps.
  - `simulStepCount::Int`: The number of simultaneous steps.
  - `evCount::Int`: The number of events.
  - `numStateSteps::Vector{Int}`: A vector holding the number of steps for each state update.
  - `numInputSteps::Vector{Int}`: A vector holding the number of steps for each input update.
"""
struct Stats
  totalSteps::Int
  simulStepCount::Int
  evCount::Int
  numStateSteps ::Vector{Int}
  numInputSteps ::Vector{Int}
end

"""
    LightSol{T,O}
A struct that holds the solution of a system of ODEs. It has the following fields:
  - `size:` The number of continuous variables T
  - `order:` The order of the algorithm O
  - `savedTimes:` A vector of vectors of Float64 that holds the times at which the continuous variables were saved
  - `savedVars:` A vector of vectors of Float64 that holds the values of the continuous variables at the times they were saved
  - `algName:` The name of the algorithm used to solve the system
  - `sysName:` The name of the system
  - `absQ:` The absolute tolerance used in the simulation
  - `stats:` A Stats struct that holds the statistics of the simulation
  - `ft:` The final time of the simulation
"""
struct LightSol{T,O}<:Sol{T,O}
  size::Val{T}
  order::Val{O}
  savedTimes::Vector{Vector{Float64}}
  savedVars::Vector{Vector{Float64}}
  algName::String
  sysName::String
  absQ::Float64
  stats::Stats
  ft::Float64
end
"""
    createSol(::Val{T}, ::Val{O}, savedTimes::Vector{Vector{Float64}}, savedVars::Vector{Vector{Float64}}, solver::String, nameof_F::String, absQ::Float64, stats::Stats, ft::Float64) where {T,O}

Creates a LightSol struct with the given parameters.

# Arguments
- `::Val{T}`: The number of continuous variables.
- `::Val{O}`: The order of the algorithm.
- `savedTimes::Vector{Vector{Float64}}`: A vector of vectors of times at which the continuous variables were saved.
- `savedVars::Vector{Vector{Float64}}`: A vector of vectors of values of the continuous variables at the times they were saved.
- `solver::String`: The name of the algorithm used to solve the system.
- `nameof_F::String`: The name of the system.
- `absQ::Float64`: The absolute tolerance used in the simulation.
- `stats::Stats`: A Stats struct that holds the statistics of the simulation.
- `ft::Float64`: The final time of the simulation.

# Returns
- A LightSol struct.
"""
@inline function createSol(::Val{T},::Val{O}, savedTimes:: Vector{Vector{Float64}},savedVars :: Vector{Vector{Float64}},solver::String,nameof_F::String,absQ::Float64,stats::Stats,ft::Float64)where {T,O}
  sol=LightSol(Val(T),Val(O),savedTimes, savedVars,solver,nameof_F,absQ,stats,ft)
end
function getindex(s::Sol, i::Int64)#helper for calling savedTimes & savedVars
  if i==1
     return s.savedTimes
  else
     return s.savedVars
  end
end
"""
    evaluateSol(sol::Sol{T,O}, index::Int, t::Float64) where {T,O}

Evaluates the solution at a given time `t` for a specified variable index.

# Arguments
- `sol::Sol{T,O}`: The solution struct.
- `index::Int`: The index of the variable to evaluate. If `index` is 0, evaluates all variables.
- `t::Float64`: The time at which to evaluate the solution.

# Returns
- The value of the specified variable at time `t`, or a vector of values if `index` is 0.

# Throws
- An error if the given time `t` is outside the solution range.
"""
@inline function evaluateSol(sol::Sol{T,O},index::Int,t::Float64)where {T,O}
  (t>sol.ft) && error("given point is outside the solution range! Verify where you want to evaluate the solution")
  if index==0
    eval_All=Vector{Float64}(undef, T)
    for index=1:T
      eval_All[index] =evaluateSol(sol,index,t)
    end
    eval_All
  else
  x=sol[2][index][end] 
  for i=2:length(sol[1][index])#savedTimes after the init time...init time is at index i=1
      if sol[1][index][i]>t # i-1 is closest lower point
        f1=sol[2][index][i-1];f2=sol[2][index][i];t1=sol[1][index][i-1] ;t2=sol[1][index][i]# find x=f(t)=at+b...linear interpolation
        a=(f2-f1)/(t2-t1)
        b=(f1*t2-f2*t1)/(t2-t1)
        x=a*t+b
        return x#taylor evaluation after small elapsed with the point before (i-1)
      elseif sol[1][index][i]==t # i-1 is closest lower point
        x=sol[2][index][i]
        return x
      end
  end
  return x #if var never changed then return init cond, or if t>lastSavedTime for this var then return last value
end
end
"""
    solInterpolated(sol::Sol{T,O},step::Float64) where {T,O}

Construct a new solution by interpolating the current solution at each step for all variables.

# Arguments
- `sol::Sol{T,O}`: The solution struct.
- `step::Float64`: the step size at which to generate the new solution.

# Returns
- A new solution that contains information at each step size..

"""
function solInterpolated(sol::Sol{T,O},step::Float64) where {T,O}
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
  interpValues=Vector{Vector{Float64}}(undef, T)
  for index=1:T
      interpValues[index]=Float64[]
      push!(interpValues[index],sol[2][index][1]) #1st element is the init cond (true value)
    for i=2:numInterpPoints-1
      push!(interpValues[index],evaluateSol(sol,index,interpTimes[i]))
    end
      push!(interpValues[index],sol[2][index][end]) #last pt @ft
      allInterpTimes[index]=interpTimes
  end
  createSol(Val(T),Val(O),allInterpTimes,interpValues,sol.algName,sol.sysName,sol.absQ,sol.stats,sol.ft#= ,sol.simulStepsVals,sol.simulStepsDers,sol.simulStepsVals =#)
end
"""
    solInterpolated(sol::Sol{T,O},index::Int,step::Float64) where {T,O}

Constructs a new solution by interpolating the current solution at each step for one variable.

# Arguments
- `sol::Sol{T,O}`: The solution struct.
- `index::Int`: the index of the variable to interpolate.
- `step::Float64`: the step size at which to generate the new solution.

# Returns
- A new solution that contains the information of one variable at each step size.

"""
function solInterpolated(sol::Sol{T,O},index::Int,step::Float64)where {T,O}
  interpTimes=Float64[]
  allInterpTimes=Vector{Vector{Float64}}(undef, 1)
  t=0.0  #later can change to init_time which could be diff than zero
  push!(interpTimes,t)
  while t+step<sol.ft
    t=t+step
    push!(interpTimes,t) 
  end
  push!(interpTimes,sol.ft)
  numInterpPoints=length(interpTimes)
  #interpValues=nothing
 # if sol isa LightSol
    interpValues=Vector{Vector{Float64}}(undef, 1)
  #end
 # for index=1:T
          interpValues[1]=[]
          push!(interpValues[1],sol[2][index][1]) #1st element is the init cond (true value)
        for i=2:numInterpPoints-1
          push!(interpValues[1],evaluateSol(sol,index,interpTimes[i]))
        end
          push!(interpValues[1],sol[2][index][end]) #last pt @ft
          allInterpTimes[1]=interpTimes
  #end
  createSol(Val(T),Val(O),allInterpTimes,interpValues,sol.algName,sol.sysName,sol.absQ,sol.stats,sol.ft#= ,sol.simulStepsVals,sol.simulStepsDers,sol.simulStepsVals =#)
end


(sol::Sol)(t::Float64,index::Int) = evaluateSol(sol,index,t)

(sol::Sol)(t::Float64;idxs=0::Int) = evaluateSol(sol,idxs,t)


"""
    show(io::IO, a::Stats)

Displays the statistics of the simulation.

# Arguments
- `io::IO`: The IO stream to write to.
- `a::Stats`: The Stats struct containing the simulation statistics.

# Prints
- The total simulation steps.
- The simultaneous steps.
- The number of events.
- The number of state steps.
- The number of input steps.
"""
function show(io::IO, a::Stats)
  println(io,"")
  println(io,"The total simulation steps: $(a.totalSteps)")
  println(io,"The simultaneous  steps: $(a.simulStepCount)")
  println(io,"The number of events: $(a.evCount)") 
  println(io,"The number of state steps: $(a.numStateSteps)")
  println(io,"The number of input steps: $(a.numInputSteps)")

end

