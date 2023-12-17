
using BenchmarkTools
using Plots;
#gr();



using DifferentialEquations
#using qssv01
using BenchmarkTools
using Plots
function odeDiffEquPackage()
    function f(du,u,p,t)
       
        du[1] = sin(t)
       
    end
   
    u0 = [0.5]
    tspan = (0.0,30.0)
    #p = -9.8
    prob = ODEProblem(f,u0,tspan)
    #sol = solve(prob,Tsit5(),callback=cb)
    sol = solve(prob,BS3())
    p1=plot!(sol,marker=(:circle),markersize=2)
   # p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.0))
  # p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.06e-1))
   savefig(p1, "bs3sin")
 end
#@btime
 odeDiffEquPackage() 