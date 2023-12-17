



using DifferentialEquations
#using qssv01
using BenchmarkTools
using Plots
function odeDiffEquPackage()
    function f(du, u, d, t)
        C = 1e-4; L = 1e-4; R = 10;U = 24.0; T = 1e-4; DC = 0.5; ROn = 1e-5;ROff = 1e5;
        du[1] =(-((u[1]*d[2]-24.0)*d[1]/(d[1]+d[2]))-u[2])/L
        du[2]=(u[1]-u[2]/R)/C
    end
    function condition(out, u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        out[1] = t-d[3]
        out[2] = t-d[4]-0.5*1e-4
        out[3] = d[5]*((u[1]*d[2]-24.0)/(d[1]+d[2]))+(1.0-d[5])*((u[1]*d[2]-24.0)*d[1]/(d[1]+d[2]))
    end
    
    function affect!(integrator, idx)
        if idx == 1
           @show integrator.t-d[3]
            if integrator.t-d[3]>0.0 
                d[4]=d[3]
                d[3]=d[3]+1e-4
                d[2]=1e-5
            end
        elseif idx == 2
            @show integrator.t-d[4]-0.5*1e-4
           
            if integrator.t-d[4]-0.5*1e-4>0.0 
                d[2]=1e5
            end  
        elseif idx == 3
            @show d[5]*((integrator.u[1]*d[2]-24.0)/(d[1]+d[2]))+(1.0-d[5])*((integrator.u[1]*d[2]-24.0)*d[1]/(d[1]+d[2]))
            if d[5]*((integrator.u[1]*d[2]-24.0)/(d[1]+d[2]))+(1.0-d[5])*((integrator.u[1]*d[2]-24.0)*d[1]/(d[1]+d[2]))>0
                d[1]=1e-5
                d[5]=1.0
              else
                d[1]=1e5
                d[5]=0.0
              end 
        end
    end
   
    cb = VectorContinuousCallback(condition, affect!, 3)

    u0 = [0.0, 0.0]
    tspan = (0.0, 0.0025)
    d = [1e5,1e-5,1e-4,0.0,0.0]
      
    
    prob = ODEProblem(f, u0, tspan, d)
    sol = solve(prob, Tsit5(), callback = cb, reltol=1e-3,abstol=1e-4#= dt = 1e-3, adaptive = false =#)
    
   # p1=plot!(sol, idxs = (1, 3));
   # p1=plot!(sol);

   # p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.0))
   #p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.06e-1))
  # savefig(p1, "tsit5f_buck")
 end
#@btime
 odeDiffEquPackage() 