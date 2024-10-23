using QuantizedSystemSolver
function test()
        odeprob = NLodeProblem(quote 
        name=(sysd0,)
        u = [3.0,20.0]
        discrete=[0.0]
        du[1] = u[2]
        du[2] = -9.8+discrete[1]*u[1]
        if -u[1]>0.0
            u[2]=-0.9*u[2]
        end
    end)  
    tspan=(0.0,20.0)
   # sol=solve(odeprob,qss2(),tspan)
    sol= solve(odeprob,nmliqss1(),tspan,abstol=1e-3,reltol=1e-2)    
    @show sol.stats.totalSteps
    save_Sol(sol)
  
    #getAverageErrorByRefs(solRef::Vector{Any},solmliqss::Sol{T,O})
end
test()


