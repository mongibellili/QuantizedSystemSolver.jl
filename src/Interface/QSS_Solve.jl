
 # user does not provide solver. default mliqss2
 
 """solve(prob::NLODEProblem{PRTYPE,T,Z,D,CS},tspan::Tuple{Float64, Float64};sparsity::Val{Sparsity}=Val(false)::Float64,saveat=1e-9::Float64::Float64,abstol=1e-4::Float64,reltol=1e-3::Float64,maxErr=Inf::Float64) where {PRTYPE,T,Z,D,CS,Sparsity}    
   
   dispatches on a specific integrator based on the algorithm provided.
   After the simulation, the solution is returned as a Solution object.
   """
function solve(prob::NLODEProblem{PRTYPE,T,Z,D,CS},tspan::Tuple{Float64, Float64};sparsity::Val{Sparsity}=Val(false)::Float64,saveat=1e-9::Float64::Float64,abstol=1e-4::Float64,reltol=1e-3::Float64,maxErr=Inf::Float64) where {PRTYPE,T,Z,D,CS,Sparsity}    
   solve(prob,QSSAlgorithm(Val(:nmliqss),Val(2)),tspan;sparsity=sparsity,saveat=saveat,abstol=abstol,reltol=reltol,maxErr=maxErr)  
end
#main solve interface
function solve(prob::NLODEProblem{PRTYPE,T,Z,D,CS},al::QSSAlgorithm{SolverType, OrderType},tspan::Tuple{Float64, Float64};sparsity::Val{Sparsity}=Val(false),saveat=1e-9::Float64,abstol=1e-4::Float64,reltol=1e-3::Float64,maxErr=Inf::Float64) where{PRTYPE,T,Z,D,CS,SolverType,OrderType,Sparsity}    
   custom_Solve(prob,al,Val(Sparsity),tspan[2],saveat,tspan[1],abstol,reltol,maxErr)
end
#default solve method: this is not to be touched...extension or modification is done through creating another custom_solve with different PRTYPE
function custom_Solve(prob::NLODEProblem{PRTYPE,T,Z,D,CS},al::QSSAlgorithm{Solver, Order},::Val{Sparsity},finalTime::Float64,saveat::Float64,initialTime::Float64,abstol::Float64,reltol::Float64,maxErr::Float64) where{PRTYPE,T,Z,D,CS,Solver,Order,Sparsity}
    # sizehint=floor(Int64, 1.0+(finalTime/saveat)*0.6)
    commonQSSdata=createCommonData(prob,Val(Order),finalTime,saveat, initialTime,abstol,reltol,maxErr)
    jac=getClosure(prob.jac)::Function #if in future jac and SD are different datastructures
    SD=getClosure(prob.SD)::Function
    
    if Solver==:qss
        integrate(al,commonQSSdata,prob,prob.eqs,jac,SD)
    else
          liqssdata=createLiqssData(prob,Val(Sparsity),Val(T),Val(Order))
         specialLiqssData=createSpecialLiqssData(Val(T))
         if VERBOSE println("begin dispatch on liqss algorithm") end
         #integrate()
         integrate(al,commonQSSdata,liqssdata,specialLiqssData,prob)
         integrate(al,commonQSSdata,liqssdata,specialLiqssData,prob,prob.eqs,jac,SD,prob.exactJac)
        #= if Solver==:nmliqss
             nmLiQSS_integrate(commonQSSdata,liqssdata,specialLiqssData,prob,prob.eqs,jac,SD,prob.exactJac)
        elseif Solver==:nliqss
            nLiQSS_integrate(commonQSSdata,liqssdata,specialLiqssData,prob,prob.eqs,jac,SD,prob.exactJac)
        elseif Solver==:mliqss
            mLiQSS_integrate(commonQSSdata,liqssdata,specialLiqssData,prob,prob.eqs,jac,SD,prob.exactJac)
        elseif Solver==:liqss
             LiQSS_integrate(commonQSSdata,liqssdata,specialLiqssData,prob,prob.eqs,jac,SD,prob.exactJac)     
        end =#
    end
 end



 function getClosure(jacSD::Function)::Function # 
   function closureJacSD(i::Int)
        jacSD(i)
   end
   return closureJacSD
 end

 function getClosure(jacSD::Vector{Vector{Int}})::Function
    function closureJacSD(i::Int)
         jacSD[i]
    end
    return closureJacSD
  end




#helper methods...not to be touched...extension can be done through creating others via specializing on one PRTYPE or more of the symbols (PRTYPE,T,Z,D,Order) if in the future...
#################################################################################################################################################################################
function createCommonData(prob::NLODEProblem{PRTYPE,T,Z,D,CS},::Val{Order},finalTime::Float64,saveat::Float64,initialTime::Float64,abstol::Float64,reltol::Float64,maxErr::Float64)where{PRTYPE,T,Z,D,CS,Order}
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
    commonQSSdata= CommonQSS_data(quantum,x,q,tx,tq,d,nextStateTime,nextInputTime ,nextEventTime , t, integratorCache,taylorOpsCache,finalTime,saveat, initialTime,abstol,reltol,maxErr,savedTimes,savedVars)
end





#no sparsity
function createLiqssData(prob::NLODEProblem{PRTYPE,T,Z,D,CS},::Val{false},::Val{T},::Val{Order})where{PRTYPE,T,Z,D,CS,Order}
    if VERBOSE println("begin create Liqss data") end
    a = Vector{Vector{Float64}}(undef, T)
   # u=Vector{Vector{MVector{Order,Float64}}}(undef, T)
    qaux = Vector{MVector{Order,Float64}}(undef, T)
    dxaux=Vector{MVector{Order,Float64}}(undef, T)
    olddx = Vector{MVector{Order,Float64}}(undef, T)
    olddxSpec = Vector{MVector{Order,Float64}}(undef, T)
     for i=1:T
       #=  temparr=Vector{MVector{Order,Float64}}(undef, T)
        for j=1:T
            temparr[j]=zeros(MVector{Order,Float64})
        end
        u[i]=temparr =#
        a[i]=zeros(T)
        qaux[i]=zeros(MVector{Order,Float64})
        olddx[i]=zeros(MVector{Order,Float64})
        dxaux[i]=zeros(MVector{Order,Float64})
        olddxSpec[i]=zeros(MVector{Order,Float64})

    end
    if VERBOSE println("END create Liqss data") end
    liqssdata= LiQSS_data(Val(false),a#= ,u =#,qaux,olddx,dxaux,olddxSpec)
end

#to be removed if sparsity did not help
#= function createLiqssData(prob::NLODEProblem{PRTYPE,T,Z,D,CS},::Val{true},::Val{T},::Val{Order})where{PRTYPE,T,Z,D,CS,Order}
    a = Vector{Vector{Float64}}(undef, T)
    u=Vector{Vector{MVector{Order,Float64}}}(undef, T)
    qaux = Vector{MVector{Order,Float64}}(undef, T)
    olddx = Vector{MVector{Order,Float64}}(undef, T)
    olddxSpec = Vector{MVector{Order,Float64}}(undef, T)
    jacDim=prob.jacDim
    #= @timeit  "liqsssparse" =# for i=1:T
        temparr=Vector{MVector{Order,Float64}}(undef, T)
        for j=1:T
            temparr[j]=zeros(MVector{Order,Float64})
        end
        u[i]=temparr
        a[i]=zeros(jacDim(i))
        qaux[i]=zeros(MVector{Order,Float64})
        olddx[i]=zeros(MVector{Order,Float64})
        olddxSpec[i]=zeros(MVector{Order,Float64})

    end
    liqssdata= LiQSS_data(Val(true),a,u,qaux,olddx,olddxSpec)
end =#


function createSpecialLiqssData(::Val{T})where{T}
    cacheA=zeros(MVector{1,Float64})
    direction= zeros(T)
    qminus= zeros(T)
    buddySimul=zeros(MVector{2,Int})
    prevStepVal = zeros(T)
    specialLiqssData=SpecialLiQSS_data(cacheA,direction,qminus,buddySimul,prevStepVal)
end





# get init conds for normal vect...getinitcond for fun can be found with qssnlsavedprob file
function getInitCond(prob::NLODEContProblem,i::Int)
    return prob.initConditions[i]
end





















