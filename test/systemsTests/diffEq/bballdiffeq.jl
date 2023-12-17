#using OrdinaryDiffEq # stuck with events
using DifferentialEquations
#using qssv01
using BenchmarkTools
using Plots
function odeDiffEquPackage()
    function f(du,u,p,t)
        du[1] = u[2]
        du[2] = -9.8
    end
    function condition(u,t,integrator) # Event when event_f(u,t) == 0
        u[1]
    end
    function affect!(integrator)
        integrator.u[2] = -integrator.u[2]
    end
    #cb = ContinuousCallback(condition,affect!,save_positions=(true,true),interp_points=100)
    affect_neg! =affect!
    cb=ContinuousCallback(condition,affect!,affect_neg!,
                   initialize = SciMLBase.INITIALIZE_DEFAULT,
                   finalize = SciMLBase.FINALIZE_DEFAULT,
                   idxs = nothing,
                   rootfind=SciMLBase.RightRootFind,
                   save_positions=(true,true),
                   interp_points=1000,
                   abstol=1e-14,reltol=0,repeat_nudge=1//100)
    u0 = [20.0,0.0]
    tspan = (0.0,30.0)
    #p = -9.8
    prob = ODEProblem(f,u0,tspan)
    #sol = solve(prob,Tsit5(),callback=cb)
    sol = solve(prob,BS3(),callback=cb)
   # p1=plot!(sol,marker=(:circle),markersize=2)
   # p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.0))
   p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.06e-1))
   savefig(p1, "tsit5_bball_savePOS_marker_zoom_WebSite_diffeq")
 end
#@btime
 odeDiffEquPackage() 