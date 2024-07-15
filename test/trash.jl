

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
using Plots
#= odeprob = NLodeProblem(quote
    name=(sysN13,)
    u = [1.0, 0.0]
  
    du[1] = -u[2]
    du[2]=u[1]
    
   
end)  
tspan=(0.0,6.0)

sol=solve(odeprob,nmliqss2(),tspan)
p1=plot(sol,layout=(2,1))
savefig(p1,"testplotlayout") 
 =#


#= p1=plot(sol,idxs=(0,2))
savefig(p1,"testplot(vars0,2)") 
p1=plot(sol,idxs=(1,2))
savefig(p1,"testplot(vars1,2)") 
p1=plot(sol,idxs=(0,1,2))
savefig(p1,"testplot(vars0,1,2)")  =#
#= p1=plot(sol,idxs=(1,2,3))
savefig(p1,"testplot(vars1,2,3)")  =#




p1=plot()
f(x)=1+cos(x)
g(x)=sin(x)
vec=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
fvec=map(f,vec)
gvec=map(g,vec)
p1=plot(vec,gvec)
p2=plot(vec,fvec)
p=plot(p1,p2,layout=2)
savefig(p,"LAYOUTS2")
#= @show sol(0.5,idxs=2)
@show sol(0.5) =#
#= @show sol.stats
print(sol.stats) =#
