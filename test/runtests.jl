using QuantizedSystemSolver
using Test
using Plots
using BSON
@testset "QuantizedSystemSolver.jl" begin
    # Write your tests here.
    include("./exampleTest.jl")
    include("./Taylor0testFunctions.jl")
    include("./Taylor0testArithmetic.jl")

    function sysb53(du,u,p,t)
        du[1] = -20.0*u[1]-80.0*u[2]+1600.0
        du[2] =1.24*u[1]-0.01*u[2]+0.2
    end
    u = [-1.0, -2.0]
    tspan=(0.0,1.0)
    odeprob=ODEProblem(sysb53,u,tspan)
    @test typeof(odeprob) <: QuantizedSystemSolver.NLODEProblem{1,2,0,0,3} 
    @test odeprob.prname == :sysb53
    @test odeprob.initConditions == [-1.0, -2.0]
    @test odeprob.jac == [[1,2], [1,2]] || odeprob.jac == [[2,1], [2,1]]
    @test odeprob.SD == [[1,2], [1,2]] || odeprob.SD == [[2,1], [2,1]]
    @test typeof(odeprob.eqs) <: Function
    @test typeof(odeprob.exactJac) <: Function
    @test odeprob.a == Val(2)
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
    d=[0.0]
     odeprob.eqs(1, q, t,d, cache)
     @test cache[1][0]==-20.0*1.0-80.0*2.0+1600.0

     sol=solve(odeprob,nmliqss1())
     @test sol.algName == "nmliqss1"
     @test 18.8<sol(0.5,idxs=2)<19.2
     sol.stats

 
   
end
