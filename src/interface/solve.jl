 # user does not provide solver. default mliqss2 is used
function solve(prob::NLODEProblem{F,PRTYPE,T,D,Z,CS},tspan::Tuple{Float64, Float64};sparsity::Val{Sparsity}=Val(false),saveat=Inf::Float64,abstol=1e-4::Float64,reltol=1e-3::Float64,maxErr=Inf::Float64,maxiters=Int(1e7)::Int) where {F,PRTYPE,T,D,Z,CS,Sparsity}    
   solve(prob,QSSAlgorithm(Val(:nmliqss),Val(2)),tspan;sparsity=sparsity,saveat=saveat,abstol=abstol,reltol=reltol,maxErr=maxErr,maxiters=maxiters)  
end
"""
    solve(prob::NLODEProblem{F,PRTYPE,T,D,Z,CS},al::QSSAlgorithm{SolverType, OrderType},tspan::Tuple{Float64, Float64};sparsity::Val{Sparsity}=Val(false),saveat=Inf::Float64,abstol=1e-3::Float64,reltol=1e-2::Float64,maxErr=Inf::Float64,maxiters=Int(1e7)::Int) where{F,PRTYPE,T,D,Z,CS,SolverType,OrderType,Sparsity}     

dispatches on a specific integrator based on the algorithm provided and send a nonlinear ODE problem to the integrator.

With the exception of the argument prob and tspan, all other arguments are optional and have default values:\n
  - The algorithm defaults to nmliqss2, and it is specified by the QSSAlgorithm type, which is a composite type that has a name and an order. It can be extended independently of the solver.
  - The sparsity argument defaults to false. If true, the integrator will use a sparse representation of the Jacobian matrix (not implemented).
  - The saveat argument defaults to Inf. It specifies the time step at which the integrator will save the solution (not implemented).
  - The abstol argument defaults to 1e-4. It specifies the absolute tolerance of the integrator.
  - The reltol argument defaults to 1e-3. It specifies the relative tolerance of the integrator.
  - The maxErr argument defaults to Inf. It specifies the maximum error allowed by the integrator. This is used as an upper bound for the quantum when a variable goes large.
  - The maxiters argument defaults to 1e7. It specifies the maximum number of steps allowed by the integrator. If the user wants to extend the limit on the maximum number of steps, this argument can be used. 
After the simulation, the solution is returned as a Solution object.
"""
function solve(prob::NLODEProblem{F,PRTYPE,T,D,Z,CS},al::QSSAlgorithm{SolverType, OrderType},tspan::Tuple{Float64, Float64};sparsity::Val{Sparsity}=Val(false),saveat=Inf::Float64,abstol=1e-3::Float64,reltol=1e-2::Float64,maxErr=Inf::Float64,maxiters=Int(1e7)::Int) where{F,PRTYPE,T,D,Z,CS,SolverType,OrderType,Sparsity}    
   custom_Solve(prob,al,Val(Sparsity),tspan[2],saveat,tspan[1],abstol,reltol,maxErr,maxiters)
end
"""
    solve(prob::NLODEProblem{F,PRTYPE,T,D,Z,CS},al::QSSAlgorithm{SolverType, OrderType};sparsity::Val{Sparsity}=Val(false),saveat=Inf::Float64,abstol=1e-3::Float64,reltol=1e-2::Float64,maxErr=Inf::Float64,maxiters=Int(1e7)::Int) where{F,PRTYPE,T,D,Z,CS,SolverType,OrderType,Sparsity}

dispatches on a specific integrator based on the algorithm provided and send a nonlinear ODE problem to the integrator.

# Arguments
- `prob::NLODEProblem{F,PRTYPE,T,D,Z,CS}`: The nonlinear ODE problem to solve.
- `al::QSSAlgorithm{SolverType, OrderType}`: The QSS algorithm to use for solving the problem.
- `sparsity::Val{Sparsity}`: A type parameter indicating whether to use sparsity (default: `Val(false)`).
- `saveat::Float64`: The time interval at which to save the solution (default: `Inf`).
- `abstol::Float64`: The absolute tolerance for the solver (default: `1e-4`).
- `reltol::Float64`: The relative tolerance for the solver (default: `1e-3`).
- `maxErr::Float64`: The maximum allowable error (default: `Inf`).
- `maxiters::Int`: The maximum number of iterations (default: `Int(1e7)`).


"""
function solve(prob::NLODEProblem{F,PRTYPE,T,D,Z,CS},al::QSSAlgorithm{SolverType, OrderType};sparsity::Val{Sparsity}=Val(false),saveat=Inf::Float64,abstol=1e-3::Float64,reltol=1e-2::Float64,maxErr=Inf::Float64,maxiters=Int(1e7)::Int) where{F,PRTYPE,T,D,Z,CS,SolverType,OrderType,Sparsity}    
   tspan=prob.tspan
    custom_Solve(prob,al,Val(Sparsity),tspan[2],saveat,tspan[1],abstol,reltol,maxErr,maxiters)
 end
#default solve method: ...extension or modification is done through creating another custom_solve with different PRTYPE

"""
    custom_Solve(prob::NLODEProblem{F,PRTYPE,T,D,Z,CS},al::QSSAlgorithm{Solver, Order},::Val{Sparsity},finalTime::Float64,saveat::Float64,initialTime::Float64,abstol::Float64,reltol::Float64,maxErr::Float64,maxiters::Int) where{F,PRTYPE,T,D,Z,CS,Solver,Order,Sparsity}

calls the integrator to solve the nonlinear ODE problem.

# Arguments
- `prob::NLODEProblem{F,PRTYPE,T,D,Z,CS}`: The nonlinear ODE problem to solve.
- `al::QSSAlgorithm{Solver, Order}`: The QSS algorithm to use for solving the problem.
- `::Val{Sparsity}`: A type parameter indicating whether to use sparsity.
- `finalTime::Float64`: The final time for the simulation.
- `saveat::Float64`: The time interval at which to save the solution.
- `initialTime::Float64`: The initial time for the simulation.
- `abstol::Float64`: The absolute tolerance for the solver.
- `reltol::Float64`: The relative tolerance for the solver.
- `maxErr::Float64`: The maximum allowable error.
- `maxiters::Int`: The maximum number of iterations.

"""
function custom_Solve(prob::NLODEProblem{F,PRTYPE,T,D,Z,CS},al::QSSAlgorithm{Solver, Order},::Val{Sparsity},finalTime::Float64,saveat::Float64,initialTime::Float64,abstol::Float64,reltol::Float64,maxErr::Float64,maxiters::Int) where{F,PRTYPE,T,D,Z,CS,Solver,Order,Sparsity}
    if saveat!=Inf error("saveat is not implemented") end
     commonQSSdata=createCommonData(prob,Val(Order),finalTime,saveat, initialTime,abstol,reltol,maxErr,maxiters)
     jac=getClosure(prob.jac)::Function #if in future jac and SD are different datastructures
     SD=getClosure(prob.SD)::Function
    if Solver==:qss
        integrate(al,commonQSSdata,prob,prob.eqs,jac,SD)
    else
          liqssdata=createLiqssData(Val(Sparsity),Val(T),Val(Order))
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
    createCommonData(prob::NLODEProblem{F,PRTYPE,T,D,Z,CS},::Val{Order},finalTime::Float64,saveat::Float64,initialTime::Float64,abstol::Float64,reltol::Float64,maxErr::Float64,maxiters::Int) where{F,PRTYPE,T,D,Z,CS,Order}

creates the necessary data for the simulation and stores it in a CommonQSS_Data struct.

# Arguments
- `prob::NLODEProblem{F,PRTYPE,T,D,Z,CS}`: The nonlinear ODE problem to solve.
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
function createCommonData(prob::NLODEProblem{F,PRTYPE,T,D,Z,CS},::Val{Order},finalTime::Float64,saveat::Float64,initialTime::Float64,abstol::Float64,reltol::Float64,maxErr::Float64,maxiters::Int) where{F,PRTYPE,T,D,Z,CS,Order}
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
    taylorOpsCache=Array{Taylor0,1}()# cache= vector of taylor0s of size CS
    for i=1:CS
    push!(taylorOpsCache,Taylor0(zeros(Order+1),Order))
    end
    commonQSSdata= CommonQSS_Data(quantum,x,q,tx,tq,d,nextStateTime,nextInputTime ,nextEventTime , t,taylorOpsCache,finalTime, initialTime,abstol,reltol,maxErr,maxiters,savedTimes,savedVars)
end


"""
    createLiqssData(::Val{false},::Val{T},::Val{Order}) where {T,Order}

Creates LIQSS-specific data required for solving an ODE problem.

# Arguments
- `::Val{Sparsity}`: A type parameter indicating whether to use sparsity.
- `::Val{T}`: The number of continuous variables.
- `::Val{Order}`: The order of the algorithm.

# Returns
- A data structure containing LIQSS-specific data required for liQSS algorithms.

"""
function createLiqssData(::Val{false},::Val{T},::Val{Order}) where {T,Order}
    qaux = Vector{MVector{Order,Float64}}(undef, T)
    dxaux=Vector{MVector{Order,Float64}}(undef, T)
     for i=1:T
        qaux[i]=zeros(MVector{Order,Float64})
        dxaux[i]=zeros(MVector{Order,Float64})
    end
    cacheA=zeros(MVector{1,Float64})
    liqssdata= LiQSS_Data(Val(false),cacheA,qaux,dxaux)
end

