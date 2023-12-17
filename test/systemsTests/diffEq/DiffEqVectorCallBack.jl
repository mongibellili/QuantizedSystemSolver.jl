using qss
using BenchmarkTools
using Plots;

#= function test()
    odeprob = @NLodeProblem begin
        #= parameter1=3000.0# cache can be dynamic....parameters take this feature
        parameter2=0.00001 =#
        L=1e3;R=10.0;uu=311.0;ROn = 1e-5; ROff = 1e5;w=314.16
        u = [0.0,0.0]
        discrete = [1e5]
        du[1] =1000.0*(uu*sin(w*t)-u[1])
        du[2]=L*(u[1]-u[2]*(R+discrete[1]));


        if -u[2]>0.0 
            discrete[1]=ROff
                                     
        end
        if -u[1]>0.0 
            discrete[1] = ROn;
            
        end
           
    end
   sol= QSS_Solve(odeprob,qss2(),dQmin=1e-3,dQrel=1e-2,finalTime=0.06)
  # @show sol
   #save_Sol(sol)
  # save_Sol(sol,xlims=(0.0,15.0) ,ylims=(-2.04e-1,2.06e-1))
end
#@time 
test() =#

using DifferentialEquations
#using qssv01
using BenchmarkTools
using Plots
function odeDiffEquPackage()
    function f(du, u, p, t)
        du[1] = u[2]
        du[2] = -p
        du[3] = u[4]
        du[4] = 0.0
    end
    function condition(out, u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        out[1] = u[1]
        out[2] = (u[3] - 10.0)u[3]
    end
    
    function affect!(integrator, idx)
        if idx == 1
            integrator.u[2] = -0.9integrator.u[2]
        elseif idx == 2
            integrator.u[4] = -0.9integrator.u[4]
        end
    end
   
    cb = VectorContinuousCallback(condition, affect!, 2)

    u0 = [50.0, 0.0, 0.0, 2.0]
    tspan = (0.0, 15.0)
    p = 9.8
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, Tsit5(), callback = cb, dt = 1e-3, adaptive = false)
    
   # p1=plot!(sol, idxs = (1, 3));
    p1=plot!(sol);

   # p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.0))
   #p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.06e-1))
   savefig(p1, "bs3_rectif")
 end
#@btime
 odeDiffEquPackage() 