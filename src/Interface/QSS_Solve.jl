 # user does not provide solver. default mliqss2 is used
function solve(prob::NLODEProblem{PRTYPE,T,Z,D,CS},tspan::Tuple{Float64, Float64};sparsity::Val{Sparsity}=Val(false),saveat=1e-9::Float64,abstol=1e-4::Float64,reltol=1e-3::Float64,maxErr=Inf::Float64,maxStepsAllowed=10000000) where {PRTYPE,T,Z,D,CS,Sparsity}    
   solve(prob,QSSAlgorithm(Val(:nmliqss),Val(2)),tspan;sparsity=sparsity,saveat=saveat,abstol=abstol,reltol=reltol,maxErr=maxErr,maxStepsAllowed=maxStepsAllowed)  
end
"""solve(prob::NLODEProblem{PRTYPE,T,Z,D,CS},al::QSSAlgorithm{SolverType, OrderType},tspan::Tuple{Float64, Float64};sparsity::Val{Sparsity}=Val(false),saveat=1e-9::Float64,abstol=1e-4::Float64,reltol=1e-3::Float64,maxErr=Inf::Float64,maxStepsAllowed=10000000) where{PRTYPE,T,Z,D,CS,SolverType,OrderType,Sparsity}  

This function dispatches on a specific integrator based on the algorithm provided.
With the exception of the argument prob and tspan, all other arguments are optional and have default values:\n
-The algorithm defaults to nmliqss2, and it is specified by the QSSAlgorithm type, which is a composite type that has a name and an order. It can be extended independently of the solver.\n
-The sparsity argument defaults to false. If true, the integrator will use a sparse representation of the Jacobian matrix (not implemented).\n
-The saveat argument defaults to 1e-9. It specifies the time step at which the integrator will save the solution (not implemented).\n
-The abstol argument defaults to 1e-4. It specifies the absolute tolerance of the integrator.\n
-The reltol argument defaults to 1e-3. It specifies the relative tolerance of the integrator.\n
-The maxErr argument defaults to Inf. It specifies the maximum error allowed by the integrator. This is used as an upper bound for the quantum when a variable goes large.\n
-The maxStepsAllowed argument defaults to 10000000. It specifies the maximum number of steps allowed by the integrator. If the user wants to extend the limit on the maximum number of steps, this argument can be used.\n
After the simulation, the solution is returned as a Solution object.
"""
function solve(prob::NLODEProblem{PRTYPE,T,Z,D,CS},al::QSSAlgorithm{SolverType, OrderType},tspan::Tuple{Float64, Float64};sparsity::Val{Sparsity}=Val(false),saveat=1e-9::Float64,abstol=1e-4::Float64,reltol=1e-3::Float64,maxErr=Inf::Float64,maxStepsAllowed=10000000) where{PRTYPE,T,Z,D,CS,SolverType,OrderType,Sparsity}    
   custom_Solve(prob,al,Val(Sparsity),tspan[2],saveat,tspan[1],abstol,reltol,maxErr,maxStepsAllowed)
end
#default solve method: this is not to be touched...extension or modification is done through creating another custom_solve with different PRTYPE
function custom_Solve(prob::NLODEProblem{PRTYPE,T,Z,D,CS},al::QSSAlgorithm{Solver, Order},::Val{Sparsity},finalTime::Float64,saveat::Float64,initialTime::Float64,abstol::Float64,reltol::Float64,maxErr::Float64,maxStepsAllowed::Int) where{PRTYPE,T,Z,D,CS,Solver,Order,Sparsity}
    # sizehint=floor(Int64, 1.0+(finalTime/saveat)*0.6)
    commonQSSdata=createCommonData(prob,Val(Order),finalTime,saveat, initialTime,abstol,reltol,maxErr,maxStepsAllowed)
    jac=getClosure(prob.jac)::Function #if in future jac and SD are different datastructures
    SD=getClosure(prob.SD)::Function
    
    if Solver==:qss
        integrate(al,commonQSSdata,prob,prob.eqs,jac,SD)
    else
          liqssdata=createLiqssData(prob,Val(Sparsity),Val(T),Val(Order))
         if VERBOSE println("begin dispatch on liqss algorithm") end
         integrate(al,commonQSSdata,liqssdata,prob,prob.eqs,jac,SD,prob.exactJac)
    
    end
 end
#=  function getClosure(jacSD::Function)::Function # 
   function closureJacSD(i::Int)
        jacSD(i)
   end
   return closureJacSD
 end =#
 function getClosure(jacSD::Vector{Vector{Int}})::Function
    function closureJacSD(i::Int)
         jacSD[i]
    end
    return closureJacSD
  end
#helper methods...not to be touched...extension can be done through creating others via specializing on one PRTYPE or more of the symbols (PRTYPE,T,Z,D,Order) if in the future...
#################################################################################################################################################################################
function createCommonData(prob::NLODEProblem{PRTYPE,T,Z,D,CS},::Val{Order},finalTime::Float64,saveat::Float64,initialTime::Float64,abstol::Float64,reltol::Float64,maxErr::Float64,maxStepsAllowed::Int)where{PRTYPE,T,Z,D,CS,Order}
  if VERBOSE println("begin create common data") end
    quantum =  zeros(T)
    x = Vector{Taylor0}(undef, T)
    q = Vector{Taylor0}(undef, T)
    savedTimes=Vector{Vector{Float64}}(undef, T)
    savedVars = Vector{Vector{Float64}}(undef, T)
    nextStateTime =  zeros(T)
    nextInputTime =   zeros(T)
    tx = zeros(T)
    tq =  zeros(T)
    nextEventTime=@MVector zeros(Z)# only Z number of zcf is usually under 100...so use of SA is ok
    t = Taylor0(zeros(Order + 1), Order)
    t[1]=1.0
    t[0]=initialTime
    integratorCache=Taylor0(zeros(Order+1),Order) #for integratestate only
    d=zeros(D)
    for i=1:D
        d[i]=prob.discreteVars[i]
    end
    for i = 1:T
        nextInputTime[i]=Inf
        nextStateTime[i]=Inf
        x[i]=Taylor0(zeros(Order + 1), Order) 
        x[i][0]= getInitCond(prob,i)        # x[i][0]= prob.initConditions[i] if to remove saving as func
        q[i]=Taylor0(zeros(Order+1), Order)#q normally 1order lower than x but since we want f(q) to  be a taylor that holds all info (1,2,3), lets have q of same Order and not update last coeff        
        tx[i] = initialTime
        tq[i] = initialTime
        savedTimes[i]=Vector{Float64}()
        savedVars[i]=Vector{Float64}()
    end
    taylorOpsCache=Array{Taylor0,1}()# cache= vector of taylor0s of size CS
    for i=1:CS
    push!(taylorOpsCache,Taylor0(zeros(Order+1),Order))
    end
    if VERBOSE println("END create common data") end
    commonQSSdata= CommonQSS_data(quantum,x,q,tx,tq,d,nextStateTime,nextInputTime ,nextEventTime , t, integratorCache,taylorOpsCache,finalTime,saveat, initialTime,abstol,reltol,maxErr,maxStepsAllowed,savedTimes,savedVars)
end
#no sparsity
function createLiqssData(prob::NLODEProblem{PRTYPE,T,Z,D,CS},::Val{false},::Val{T},::Val{Order})where{PRTYPE,T,Z,D,CS,Order}
    qaux = Vector{MVector{Order,Float64}}(undef, T)
    dxaux=Vector{MVector{Order,Float64}}(undef, T)
     for i=1:T
        qaux[i]=zeros(MVector{Order,Float64})
        dxaux[i]=zeros(MVector{Order,Float64})
    end
    cacheA=zeros(MVector{1,Float64})
    liqssdata= LiQSS_data(Val(false),cacheA,qaux,dxaux)
end
# get init conds for normal vect...getinitcond for fun can be found with qssnlsavedprob file
function getInitCond(prob::NLODEContProblem,i::Int)
    return prob.initConditions[i]
end
