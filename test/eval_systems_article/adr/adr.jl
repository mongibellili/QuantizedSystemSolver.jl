
using QuantizedSystemSolver
using BenchmarkTools
#using XLSX
using Plots

function adr(du,u,p,t)
    _dx=100.0
    a=1.0
    d=0.1
    r=50000.0
    du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
    for k in 2:999  
        du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
    end 
    du[1000]=-a*_dx*(u[1000]-u[999])+d*_dx*_dx*(2.0*u[999]-2.0*u[1000])+r*u[1000]*u[1000]*(1.0-u[1000]) 
end
tspan = (0.0,10.0)


reltol=1e-3
abstol=1e-3
alg=nmliqss2()
for i=2:3
    u0=zeros(1000)
    u0[1:333].=1.0
    #Construct the problem
    odeprob = ODEProblem(adr,u0,tspan#= ,jac_mode=:symbolic =#) 
    sol=solve(odeprob,alg,abstol=abstol,reltol=reltol,maxiters=Int(1e8),detection=Detection(i))
    @show sol.stats.totalSteps
    @show sol.stats.simulSteps

    p1=plot(sol,idxs=[1,100,300,400,500,700,800,1000],marker=:circle,markersize=1)
    savefig(p1, "adrd01_$(alg)_r50000__$(alg)_$(tspan[2])_rel$(reltol)_detect_$(i).png")

   @btime solve($odeprob,$alg,detection=Detection($i),abstol=$abstol,reltol=$reltol) 
end
 













