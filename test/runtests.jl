using QuantizedSystemSolver
using Test

@testset "QuantizedSystemSolver.jl" begin
    # Write your tests here.
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
    tspan=(0.0,1.0)
    sol=solve(odeprob,nmliqss1(),tspan)
    @test sol.algName == "nmliqss1"
    
end
