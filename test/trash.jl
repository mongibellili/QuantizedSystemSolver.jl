

using QuantizedSystemSolver
#= odeprob = NLodeProblem(quote
    name=(lorenz,)
    u  = [1.0,0.0,0.0]
    du[1] = 10.0(u[2]-u[1])
    du[2] = u[1]*(28.0-u[3]) - u[2]
    du[3] = u[1]*u[2] - (8/3)*u[3]
end)  
tspan=(0.0,10.0)
sol=solve(odeprob,qss2(),tspan)
save_Sol(sol) =#

#= odeprob = NLodeProblem(quote
    name=(vanderpol,)
    u  = [0.0,1.7]
    du[1] = u[2]
    du[2] = (1.0-u[1]*u[1])*u[2]-u[1] 
    
end)  
tspan=(0.0,10.0)
sol=solve(odeprob,nmliqss2(),tspan)
save_Sol(sol)
 =#

#= odeprob = NLodeProblem(quote
    name=(oregonator,)
    u  = [1.0,1.0,0.0]
    du[1] = 100.8*(9.523809523809524e-5*u[2]-u[1]*u[2]+u[1]*(1.0-u[1]))
    du[2] =40320.0*(-9.523809523809524e-5*u[2]-u[1]*u[2]+u[3])
    du[3] = u[1] -u[3]
end)  
tspan=(0.0,10.0)
sol=solve(odeprob,nmliqss2(),tspan)
save_Sol(sol) =#
using Test
odeprob = NLodeProblem(quote
    name=(sysN13,)
    u = [1.0, 0.0,1.0, 0.0,1.0]
    discrete = [0.5]
    du[1] = t
    
    for k in 2:5 
        du[k]=discrete[1]*(u[k]-u[k-1]) ;
    end 
    if t-5.0>0.0
        discrete[1]=0.0
    end
    if t-3.0>0.0
        u[1] = 1.0
        u[2] = 0.0
        u[3] = 1.0
        u[4] = 0.0
        u[5] = 1.0
        discrete[1]=1.0
        
    end
end)  
tspan=(0.0,6.0)
sol=solve(odeprob,qss1(),tspan)
sol=solve(odeprob,liqss1(),tspan)
sol=solve(odeprob,nmliqss2(),tspan)
@test 1.1<sol(1,0.5)<1.3
@test -0.35<sol(2,0.5)<-0.2