#= using qss
using BenchmarkTools =#
#= using Plots;
gr(); =#

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
    function f(du,u,p,t)
        
        L=1e3;R=10.0;uu=311.0;ROn = 1e-5; ROff = 1e5;w=314.16
        du[1] = 1000.0*(311.0*sin(314.16*t)-u[1])
        du[2] = L*(u[1]-u[2]*(R))*p
    end
    function condition1(u,t,integrator) # Event when event_f(u,t) == 0
        u[1]
    end
    #= function condition2(u,t,integrator) # Event when event_f(u,t) == 0
        u[1]
    end =#
    function affect_neg!(integrator)
        #integrator.u[2] = 0.0
        integrator.p=0.0
            
    end
    function affect_pos!(integrator)
        
        integrator.p=1.0
            
    end
    cb=ContinuousCallback(condition1,nothing,affect_neg!)
    cb2=ContinuousCallback(condition1,affect_pos!,nothing)
    cbs = CallbackSet(cb, cb2)
    u0 = [0.5,0.1]
    tspan = (0.0,0.06)
    p = 1.0
    prob = ODEProblem(f,u0,tspan,p)
    #sol = solve(prob,Tsit5(),callback=cb)
    sol = solve(prob,BS3(),callback=cbs)
   p1=plot!(sol,marker=(:circle),markersize=2)
   # p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.0))
   #p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.06e-1))
   savefig(p1, "bs3_rectif")
 end
#@btime
 odeDiffEquPackage() 