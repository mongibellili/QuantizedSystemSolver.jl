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
    @test odeprob.prname == :sysb53
end
