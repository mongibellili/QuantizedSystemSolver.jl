using QuantizedSystemSolver

#= odeprob = NLodeProblem(quote
    name=(sysd0,)
u = [10.0]
du[1] =-u[1]
end)  
tspan=(0.0,10.0)
sol=solve(odeprob,qss2(),tspan)
save_Sol(sol) =#
#= odeprob = NLodeProblem(quote 
    name=(sysd0,)
    u = [10.0]
    discrete=[-1e5]
    du[1] =-u[1]/6.0
    if t-4.0>0.0
        discrete[1]=0.0
    end
    if t-4.00000001>0.0
        discrete[1]=-1e5
    end
    if discrete[1]+(4.0-u[1])>0.0
        u[1]=u[1]+10.0
    end
end)  
tspan=(0.0,10.0)
#@show odeprob
sol=solve(odeprob,nmliqss2(),tspan)
save_Sol(sol) =#
using Plots
odeprob = NLodeProblem(quote 
    name=(sysd0,)
    u = [50.0,0.0]
    discrete=[0.0]
    du[1] = u[2]
    du[2] = -9.8+discrete[1]*u[1]
  
    if -u[1]>0.0
        u[2]=-u[2]
    end
end)  
tspan=(0.0,15.0)
#@show odeprob
sol=solve(odeprob,qss2(),tspan)
p=plot_Sol(sol,marker=:none)
savefig(p,"testmarker")
#save_Sol(sol)