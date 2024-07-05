using QuantizedSystemSolver



function test(solvr,absTol,relTol)
#=     odeprob = NLodeProblem(quote
        name=(tyson,)
        u = [0.0,0.75,0.25,0.0,0.0,0.0]
        du[1] = u[4]-1e6*u[1]+1e3*u[2]
        du[2] =-200.0*u[2]*u[5]+1e6*u[1]-1e3*u[2]
        du[3] = 200.0*u[2]*u[5]-u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)
        du[4] =u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)-u[4]
        du[5] = 0.015-200.0*u[2]*u[5]
        du[6] =u[4]-0.6*u[6]
    end ) 
    println("start tyson solving")
    tspan=(0.0,25.0)
    sol=solve(odeprob,solvr,abstol=absTol,reltol=relTol,tspan)
    println("start saving plot")
    @show sol.totalSteps, sol.simulStepCount
    if sol.totalSteps<1000000
        save_Sol(sol)
    end =#
    odeprob = NLodeProblem(quote
    name=(sysN1,)
    u = [0.0, 1.0]
  #=   du[1] = cos(u[2])
    du[2] = sin(u[1]) =#
    du[1] = (u[2])
    du[2] = -(u[1])
    end)  
    tspan=(0.0,10.0)
    sol=solve(odeprob,solvr,abstol=absTol,reltol=relTol,tspan)
    println("start saving plot")
    @show sol.totalSteps, sol.simulStepCount
    if sol.totalSteps<1000000
        save_Sol(sol)
    end 

end

absTol=1e-5
relTol=1e-2

solvrs=[#= qss1(), =##= qss2() , =##= liqss2() =##= ,liqss2() =##= ,nmliqss1(), =#nmliqss2()]
for solvr in solvrs
    test(solvr,absTol,relTol)
end