



using DifferentialEquations
#using qssv01
using BenchmarkTools
using Plots
function odeDiffEquPackage()
    function f(du, u, p, t)
        C=1e4;L=1e4;R=0.1;uu=24.0;T=1e-4;DC=0.5;Ron = 1e-5; ROff = 1e5;N=4.0
        du[1] =C*(u[2]+u[3]-R*u[1]+u[4]+u[5])
        du[2]=(((uu/discrete[5]) - u[2]) * (discrete[5]*discrete[1]/(discrete[5]+discrete[1])) - u[1])*L
            du[3]=(((uu/discrete[6]) - u[3]) * (discrete[6]*discrete[2]/(discrete[6]+discrete[2])) - u[1])*L
            du[4]=(((uu/discrete[7]) - u[4]) * (discrete[7]*discrete[3]/(discrete[7]+discrete[3])) - u[1])*L
            du[5]=(((uu/discrete[8]) - u[5]) * (discrete[8]*discrete[4]/(discrete[8]+discrete[4])) - u[1])*L
        
    end
    function condition(out, u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        out[1] =  t-discrete[9]
        out[2] = (t-discrete[10]-0.01*1e-4)
        out[3] = (t-discrete[10]-1e-4/4.0-0.01*1e-4)
        out[4] =  (t-discrete[10]-1e-4/2.0-0.01*1e-4)
        out[5] = (t-discrete[10]-1e-4*3.0/4.0-0.01*1e-4)
        out[6] = t-discrete[10]-1e-4*3.0/4.0-0.01*1e-4
        out[7] = t -  discrete[10]-1e-4/4.0-0.5*1e-4/4-0.01*1e-4
        out[8] =  t -  discrete[10]-1e-4/2.0-0.5*1e-4/4-0.01*1e-4
        out[9] = t -  discrete[10]-1e-4*3.0/4-0.5*1e-4/4-0.01*1e-4
      
        
    end
    function condition2(out, u, t, integrator) # Event when condition(out,u,t,integrator) == 0
       
        out[1] = u[2]
        out[2] = u[3]
        out[3] = u[4]
        out[4] = u[5]
        
    end
   


    function affect!(integrator, idx)
        if idx == 1
            discrete[10]=discrete[9] 
            discrete[9]=discrete[9]+1e-4
        elseif idx == 2
            discrete[5] = 1e-5
        discrete[1] = 1e5
        elseif idx == 3
            discrete[6] = 1e-5
            discrete[2] = 1e5
        elseif idx == 4
            discrete[7] = 1e-5
            discrete[3] = 1e5
        elseif idx == 5
            discrete[8] = 1e-5
            discrete[4] = 1e5
        elseif idx == 6
            discrete[5]  = 1e5
            discrete[1]  = 1e-5
        elseif idx == 7
            discrete[6]  = 1e5
            discrete[2]  = 1e-5
        elseif idx == 8
            discrete[7]  = 1e5
            discrete[3] = 1e-5
        elseif idx == 9
            discrete[8] = 1e5
            discrete[4]  = 1e-5
        
        end
        
         
    end
    function affect2!(integrator, idx)
        if idx == 1
              discrete[1] = 1e5
        elseif idx ==2
            discrete[2] = 1e5
        elseif idx == 3
            discrete[3] = 1e5
        elseif idx == 4
            discrete[4] = 1e5
        end
        
         
    end
 
    cbs1 = VectorContinuousCallback(condition, affect!,nothing, 9)
    cbs2 = VectorContinuousCallback(condition2,nothing, affect2!, 4)
    cbs= CallbackSet(cbs1, cbs2)
    tspan = (0.0, 0.0003)
   
    u0 = [0.0,0.0,0.0,0.0,0.0]
    discrete = [1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e-8,0.0]
    
    prob = ODEProblem(f, u0, tspan, discrete)
    sol = solve(prob,Tsit5() #= QBDF2() =#, callback = cbs, reltol=1e-3,abstol=1e-4#= dt = 1e-3, adaptive = false =#)
    
    p1=plot!(sol);
  # @show p[6]
   # p1=plot!(sol,ylims=(-2.04e-1,2.06e-1))

   # p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.0))
   #p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.06e-1))
   savefig(p1, "QBDF2_34_interleaved_vector_ft0025_ogppos")
   #savefig(p1, "zoom1QBDF2_34_buck_ft0025")
  # p2=plot!(sol,ylims=(-2.04e-4,2.06e-4))
   #savefig(p2, "zoom2QBDF2_34_buck_ft0025")
 end
#@btime
 odeDiffEquPackage() 