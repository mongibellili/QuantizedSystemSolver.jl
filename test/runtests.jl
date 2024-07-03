using QuantizedSystemSolver
using Test

@testset "QuantizedSystemSolver.jl" begin
    # Write your tests here.
    include("./Taylor0testFunctions.jl")
    include("./Taylor0testArithmetic.jl")
    odeprob = NLodeProblem(quote
        #sys b53
        name=(sysb53,)
        u = [-1.0, -2.0]
        du[1] = -20.0*u[1]-80.0*u[2]+1600.0
        du[2] =1.24*u[1]-0.01*u[2]+0.2
    end)  
    @test typeof(odeprob) <: NLODEProblem{1,2,0,0,3} 
    @test odeprob.prname == :sysb53
    @test odeprob.initConditions == [-1.0, -2.0]
    @test odeprob.jac == [[1,2], [1,2]] || odeprob.jac == [[2,1], [2,1]]
    @test odeprob.SD == [[1,2], [1,2]] || odeprob.SD == [[2,1], [2,1]]
    @test typeof(odeprob.eqs) <: Function
    @test typeof(odeprob.exactJac) <: Function
    @test typeof(odeprob.jacDim) <: Function
    @test odeprob.a == Val(2)
    @test odeprob.b == Val(0)
    @test odeprob.c == Val(0)
    @test odeprob.prtype == Val(1)
   
    
    Order=1
    cache=Array{Taylor0,1}()# cache= vector of taylor0s of size CS
    for i=1:3
    push!(cache,Taylor0(zeros(Order+1),Order))
    end
    q1=Taylor0([1.0,0.0], Order)
    q2=Taylor0([2.0,0.0], Order)
    q=[q1,q2];
    t=Taylor0(zeros(Order + 1), Order)
    
     odeprob.eqs(1, q, t, cache)
     @test cache[1][0]==-20.0*1.0-80.0*2.0+1600.0

     tspan=(0.0,1.0)
     sol=solve(odeprob,nmliqss1(),tspan)
     @test sol.algName == "nmliqss1"
     @test sol(2,0.5) == 19.076750604345154

     odeprob = NLodeProblem(quote   #NLodeProblem(quote ... end);
          name=(buck,)
          C = 1e-4; L = 1e-4; R = 10.0;U = 24.0; T = 1e-4; DC = 0.5; ROn = 1e-5;ROff = 1e5;
          discrete = [1e5,1e-5,1e-4,0.0,0.0];u = [0.0,0.0]
          rd=discrete[1];rs=discrete[2];nextT=discrete[3];lastT=discrete[4];diodeon=discrete[5]
          il=u[1] ;uc=u[2]
          id=(il*rs-U)/(rd+rs) # diode's current
          du[1] =(-id*rd-uc)/L
          du[2]=(il-uc/R)/C
          if t-nextT>0.0 
            lastT=nextT;nextT=nextT+T;rs=ROn
          end
          if t-lastT-DC*T>0.0 
            rs=ROff
          end                          
          if diodeon*(id)+(1.0-diodeon)*(id*rd-0.6)>0
            rd=ROn;diodeon=1.0
          else
            rd=ROff;diodeon=0.0
          end     
    end)
    tspan = (0.0, 0.001)
    sol= solve(odeprob,nmliqss2(),tspan,abstol=1e-4,reltol=1e-3)    
    @test sol(2,0.0005)≈19.209921627620943
    
end
