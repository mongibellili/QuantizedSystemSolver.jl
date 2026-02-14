
using DifferentialEquations
#using LSODA
#using Sundials
using BenchmarkTools
#using XLSX
using Plots
function adr(du,u,p,t)
    _dx=100.0
    a=1.0
    d=0.1
    r=1000.0
    du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
    for k in 2:999  
        du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
    end 
    du[1000]=-a*_dx*(u[1000]-u[999])+d*_dx*_dx*(2.0*u[999]-2.0*u[1000])+r*u[1000]*u[1000]*(1.0-u[1000]) 
end
tspan = (0.0,10.0)


reltol=1e-3
abstol=1e-3
solvers=[#= QNDF2(),TRBDF2(),Trapezoid(), =#Heun()]
for alg in solvers
  u0=zeros(1000)
  u0[1:333].=1.0
  #Construct the problem
  odeprob = ODEProblem(adr,u0,tspan#= ,jac_mode=:symbolic =#) 
  sol=solve(odeprob,alg,abstol=abstol,reltol=reltol,dense=false)
  @show sol.stats



  #p1=plot(sol,idxs=[1,100,300,400,500,700,800,1000],marker=:circle,markersize=1)
  p1=plot(sol,idxs=[1,2,3,4,5,10,20,30,50,80,100,120,140,160,180,200,220,240,280,300,330,360,400,420,450,500,520,540,560,600,620,640,680,700,750,800,850,900,950,1000],marker=:circle,markersize=1,legend=false)
  savefig(p1, "_adr_d0_1_r1000_$(typeof(alg).name.name)_$(tspan[2])_rel$(reltol).png")

  #@btime solve($odeprob,$alg,abstol=$abstol,reltol=$reltol,dense=false)  


end


