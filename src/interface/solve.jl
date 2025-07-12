"""
    solve(prob::ODEProblemData{F,PRTYPE,T,D,Z,CS},al::QSSAlgorithm{SolverType, OrderType};detection::Val{M}=Val(2),saveat::Float64=Inf,abstol::Float64=1e-3,reltol::Float64=1e-3,maxErr::Float64=1.0,maxiters::Int=Int(1e7),verbose=false::Bool) where{F,PRTYPE,T,D,Z,CS,SolverType,OrderType,M}

dispatches on a specific integrator based on the algorithm provided and sends a nonlinear ODE problem to the integrator. With the exception of the argument prob, all other arguments are optional and have default values:

# Arguments
- `prob::ODEProblemData{F,PRTYPE,T,D,Z,CS}`: The nonlinear ODE problem to solve.
- `al::QSSAlgorithm{SolverType, OrderType}`: The QSS algorithm to use for solving the problem.
- `detection::Val{M}`: A type parameter indicating which detection mechanism to use.
- `saveat::Float64`: The time interval at which to save the solution (default: `Inf`).
- `abstol::Float64`: The absolute tolerance for the solver (default: `1e-4`).
- `reltol::Float64`: The relative tolerance for the solver (default: `1e-3`).
- `maxErr::Float64`: The maximum allowable error (default: `Inf`).
- `maxiters::Int`: The maximum number of iterations (default: `Int(1e7)`).

After the simulation, the solution is returned as a Solution object.
"""
function solve(prob::ODEProblemData{F,PRTYPE,T,D,Z,CS},al::QSSAlgorithm{SolverType, OrderType};detection::Val{M}=Val(2),saveat::Float64=Inf,abstol::Float64=1e-3,reltol::Float64=1e-3,maxErr::Float64=1.0,maxiters::Int=Int(1e7),verbose=false::Bool) where{F,PRTYPE,T,D,Z,CS,SolverType,OrderType,M}    
   tspan=prob.tspan  # tspan separated before custom_solve to allow user to enter tspan with solve (not with prob)
   custom_Solve(prob,al,Val(M),tspan[2],saveat,tspan[1],abstol,reltol,maxErr,maxiters,verbose)
end

"""
    solve(prob::ODEProblemData{F,PRTYPE,T,D,Z,CS},al::QSSAlgorithm{SolverType, OrderType},tspan::Tuple{Float64, Float64};detection::Val{M}=Val(2),saveat::Float64=Inf,abstol::Float64=1e-3,reltol::Float64=1e-3,maxErr::Float64=1.0,maxiters::Int=Int(1e7),verbose=false::Bool) where{F,PRTYPE,T,D,Z,CS,SolverType,OrderType,M}     

same as the previous solve method, but with a specified time span. This method is useful when the user wants to specify the time span separately from the problem definition.
"""
function solve(prob::ODEProblemData{F,PRTYPE,T,D,Z,CS},al::QSSAlgorithm{SolverType, OrderType},tspan::Tuple{Float64, Float64};detection::Val{M}=Val(2),saveat::Float64=Inf,abstol::Float64=1e-3,reltol::Float64=1e-3,maxErr::Float64=1.0,maxiters::Int=Int(1e7),verbose=false::Bool) where{F,PRTYPE,T,D,Z,CS,SolverType,OrderType,M}    
   custom_Solve(prob,al,Val(M),tspan[2],saveat,tspan[1],abstol,reltol,maxErr,maxiters,verbose)
end

 # user does not provide solver. default mliqss2 is used. tspan provided in solve
function solve(prob::ODEProblemData{F,PRTYPE,T,D,Z,CS},tspan::Tuple{Float64, Float64};detection::Val{M}=Val(2),saveat::Float64=Inf,abstol::Float64=1e-3,reltol::Float64=1e-3,maxErr::Float64=1.0,maxiters::Int=Int(1e7),verbose::Bool=false) where {F,PRTYPE,T,D,Z,CS,M}    
   solve(prob,QSSAlgorithm(Val(:nmliqss),Val(2)),tspan;detection=detection,saveat=saveat,abstol=abstol,reltol=reltol,maxErr=maxErr,maxiters=maxiters,verbose=verbose)  
end
# user does not provide solver. default mliqss2 is used. tspan provided in prob
function solve(prob::ODEProblemData{F,PRTYPE,T,D,Z,CS};detection::Val{M}=Val(2),saveat::Float64=Inf,abstol::Float64=1e-3,reltol::Float64=1e-3,maxErr::Float64=1.0,maxiters::Int=Int(1e7)) where{F,PRTYPE,T,D,Z,CS,M}    
   solve(prob,QSSAlgorithm(Val(:nmliqss),Val(2));detection=detection,saveat=saveat,abstol=abstol,reltol=reltol,maxErr=maxErr,maxiters=maxiters) 
end


"""
    custom_Solve(prob::ODEProblemData{F,PRTYPE,T,D,Z,CS},al::QSSAlgorithm{Solver, Order},::Val{M},finalTime::Float64,saveat::Float64,initialTime::Float64,abstol::Float64,reltol::Float64,maxErr::Float64,maxiters::Int,verbose::Bool) where{F,PRTYPE,T,D,Z,CS,Solver,Order,M}

default solve method: calls the integrator to solve the nonlinear ODE problem.

# Arguments
- `prob::ODEProblemData{F,PRTYPE,T,D,Z,CS}`: The nonlinear ODE problem to solve.
- `al::QSSAlgorithm{Solver, Order}`: The QSS algorithm to use for solving the problem.
- `::Val{M}`: A type parameter indicating which detection mechanism to use.
- `finalTime::Float64`: The final time for the simulation.
- `saveat::Float64`: The time interval at which to save the solution.
- `initialTime::Float64`: The initial time for the simulation.
- `abstol::Float64`: The absolute tolerance for the solver.
- `reltol::Float64`: The relative tolerance for the solver.
- `maxErr::Float64`: The maximum allowable error.
- `maxiters::Int`: The maximum number of iterations.

"""
function custom_Solve(prob::ODEProblemData{F,PRTYPE,T,D,Z,CS},al::QSSAlgorithm{Solver, Order},::Val{M},finalTime::Float64,saveat::Float64,initialTime::Float64,abstol::Float64,reltol::Float64,maxErr::Float64,maxiters::Int,verbose::Bool) where{F,PRTYPE,T,D,Z,CS,Solver,Order,M}
    if Order==1 && abstol<=1e-4 && reltol<=1e-3
        @warn("The provided tolerance values are not suitable for the first order algorithm. Consider using (1e-3;1e-3) or higher.")
    end
    if saveat!=Inf error("saveat is not implemented yet") end
     commonQSSdata=createCommonData(prob,Val(Order),finalTime,saveat, initialTime,abstol,reltol,maxErr,maxiters,verbose)
     jac=getClosure(prob.jac)::Function #if in future jac and SD are different datastructures
     SD=getClosure(prob.SD)::Function
    if Solver==:qss
        integrate(al,commonQSSdata,prob,prob.eqs,jac,SD) 
    else
          liqssdata=createLiqssData(Val(PRTYPE),Val(M),Val(T),Val(Order))
          integrate(al,commonQSSdata,liqssdata,prob,prob.eqs,jac,SD,prob.exactJac)  
    end
 end
 function getClosure(jacSD::Vector{Vector{Int}})::Function
    function closureJacSD(i::Int)
         jacSD[i]
    end
    return closureJacSD
  end
#helper methods...extension can be done through creating others via specializing on one PRTYPE or more of the symbols (PRTYPE,T,D,Z,Order) 
#################################################################################################################################################################################
"""
    createCommonData(prob::ODEProblemData{F,PRTYPE,T,D,Z,CS},::Val{Order},finalTime::Float64,saveat::Float64,initialTime::Float64,abstol::Float64,reltol::Float64,maxErr::Float64,maxiters::Int,verbose::Bool) where{F,PRTYPE,T,D,Z,CS,Order}

creates the necessary data for the simulation and stores it in a CommonQSS_Data struct.

# Arguments
- `prob::ODEProblemData{F,PRTYPE,T,D,Z,CS}`: The nonlinear ODE problem to solve.
- `::Val{Order}`: The order of the algorithm.
- `finalTime::Float64`: The final time for the simulation.
- `saveat::Float64`: The time interval at which to save the solution.
- `initialTime::Float64`: The initial time for the simulation.
- `abstol::Float64`: The absolute tolerance for the solver.
- `reltol::Float64`: The relative tolerance for the solver.
- `maxErr::Float64`: The maximum allowable error.
- `maxiters::Int`: The maximum number of iterations.

# Returns
- A data structure containing common data required for the QSS algorithm.
"""
function createCommonData(prob::ODEProblemData{F,PRTYPE,T,D,Z,CS},::Val{Order},finalTime::Float64,saveat::Float64,initialTime::Float64,abstol::Float64,reltol::Float64,maxErr::Float64,maxiters::Int,verbose::Bool) where{F,PRTYPE,T,D,Z,CS,Order}
    quantum =  zeros(T)
    x = Vector{Taylor0}(undef, T)
    q = Vector{Taylor0}(undef, T)
    savedTimes=Vector{Vector{Float64}}(undef, T)
    savedVars = Vector{Vector{Float64}}(undef, T)
    nextStateTime =  zeros(T)
    nextInputTime =   zeros(T)
    tx = zeros(T)
    tq =  zeros(T)
    nextEventTime=@MVector zeros(Z)# only Z number is stored in Staticarray because zcf is usually under 100...
    t = Taylor0(zeros(Order + 1), Order)
    t[1]=1.0
    t[0]=initialTime
    d=zeros(D)
    for i=1:D
        d[i]=prob.discreteVars[i]
    end
    for i = 1:T
        nextInputTime[i]=Inf
        nextStateTime[i]=Inf
        x[i]=Taylor0(zeros(Order + 1), Order) 
        x[i][0]= prob.initConditions[i]      
        q[i]=Taylor0(zeros(Order+1), Order)#q normally 1-order lower than x but since we want f(q) to  be a taylor that holds all info (1,2,3), lets have q of same Order and not update last coeff        
        tx[i] = initialTime
        tq[i] = initialTime
        savedTimes[i]=Vector{Float64}()
        savedVars[i]=Vector{Float64}()
    end
    taylorOpsCache = Vector{Taylor0}(undef, CS)
    for i in 1:CS
        taylorOpsCache[i] = Taylor0(zeros(Order + 1), Order)
    end
    commonQSSdata= CommonQSS_Data(quantum,x,q,tx,tq,d,nextStateTime,nextInputTime ,nextEventTime , t,taylorOpsCache,finalTime, initialTime,abstol,reltol,maxErr,maxiters,verbose,savedTimes,savedVars)
end


"""
    createLiqssData(::Val{:symbolic},::Val{M},::Val{T},::Val{Order}) where{T,Order,M}

Creates LIQSS-specific data required for solving an ODE problem.

# Arguments
- `::Val{M}`: A type parameter indicating which detection mechanism to use.
- `::Val{T}`: The number of continuous variables.
- `::Val{Order}`: The order of the algorithm.

# Returns
- A data structure containing LIQSS-specific data required for liQSS algorithms.

"""
function createLiqssData(::Val{:symbolic},::Val{M},::Val{T},::Val{Order}) where{T,Order,M}
    qaux = Vector{MVector{Order,Float64}}(undef, T)
    dxaux=Vector{MVector{Order,Float64}}(undef, T)
     for i=1:T
        qaux[i]=zeros(MVector{Order,Float64})
        dxaux[i]=zeros(MVector{Order,Float64})
    end
    cacheA=zeros(MVector{1,Float64})
    liqssdata= AexprLiQSS_data(Val(M),1,cacheA,qaux,dxaux,1)    
end

function createLiqssData(::Val{:approximate},::Val{M},::Val{T},::Val{Order})where{T,Order,M}
    a = Vector{Vector{Float64}}(undef, T)
    qaux = Vector{MVector{Order,Float64}}(undef, T)
    dxaux=Vector{MVector{Order,Float64}}(undef, T)
    olddx = Vector{MVector{1,Float64}}(undef, T)
     for i=1:T
        a[i]=zeros(T)
        qaux[i]=zeros(MVector{Order,Float64})
        dxaux[i]=zeros(MVector{Order,Float64})
        olddx[i]=zeros(MVector{1,Float64}) 
    end
    liqssdata= AmanualLiQSS_data(Val(M),a,1,qaux,dxaux,olddx) 
end



function Detection(detectNumber::Int)
    return Val(detectNumber)
end 