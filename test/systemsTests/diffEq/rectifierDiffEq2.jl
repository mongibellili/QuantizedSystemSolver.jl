


using DifferentialEquations
#using qssv01
using BenchmarkTools
using Plots
function odeDiffEquPackage()
    function f(du,u,p,t)
        
        L=1e3;R=10.0;uu=311.0;ROn = 1e-5; ROff = 1e5;w=314.16
        du[1] = 1000.0*(311.0*sin(314.16*t)-u[1])
        du[2] = L*(u[1]-u[2]*(R))*p
    end
    function condition(u,t,integrator) # Event when event_f(u,t) == 0
        u[1]
    end
    #= function condition2(u,t,integrator) # Event when event_f(u,t) == 0
        u[1]
    end =#
    function affect(integrator)
        if integrator.u[1]>0.0
            integrator.u[2] = 0.0
            integrator.p=0.0
        else
            integrator.p=1.0
        end
            
    end
    
    cb=ContinuousCallback(condition,affect)
    
    u0 = [0.5,0.1]
    tspan = (0.0,0.06)
    p = 1.0
    prob = ODEProblem(f,u0,tspan,p)
    #sol = solve(prob,Tsit5(),callback=cb)
    sol = solve(prob,QBDF2(),callback=cb,reltol=1e-2,abstol=1e-3)
   p1=plot!(sol,marker=(:circle),markersize=2)
   # p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.0))
   #p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.06e-1))
   savefig(p1, "QBDF2_tol-2_rectif")
 end
#@btime
odeDiffEquPackage() 