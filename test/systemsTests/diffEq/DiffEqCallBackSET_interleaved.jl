



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
    function condition11( u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        t-discrete[9]
    end
    function condition21( u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        (t-discrete[10]-0.01*1e-4)
    end
    function condition22( u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        (t-discrete[10]-1e-4/4.0-0.01*1e-4)
    end
    function condition23( u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        (t-discrete[10]-1e-4/2.0-0.01*1e-4)
    end
    function condition24( u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        (t-discrete[10]-1e-4*3.0/4.0-0.01*1e-4)
    end

    function condition31( u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        t-discrete[10]-1e-4*3.0/4.0-0.01*1e-4
    end
    function condition32( u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        t -  discrete[10]-1e-4/4.0-0.5*1e-4/4-0.01*1e-4
    end
    function condition33( u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        t -  discrete[10]-1e-4/2.0-0.5*1e-4/4-0.01*1e-4
    end
    function condition34( u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        t -  discrete[10]-1e-4*3.0/4-0.5*1e-4/4-0.01*1e-4
    end
    function condition41( u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        u[2]
    end
    function condition42( u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        u[3]
   end
   function condition43( u, t, integrator) # Event when condition(out,u,t,integrator) == 0
    u[4]
end
function condition44( u, t, integrator) # Event when condition(out,u,t,integrator) == 0
    u[5]
end


    function affect11!(integrator)
       #=  println("pos at $(integrator.t)")
        @show integrator.t-p[3] =#
        discrete[10]=discrete[9] 
        discrete[9]=discrete[9]+1e-4
         
    end
    function affect21!(integrator)
        #= println("pos at $(integrator.t)")
        @show integrator.t-p[4]-0.5*1e-4 =#
        discrete[5] = 1e-5
        discrete[1] = 1e5
             
    end
    function affect22!(integrator)
      #=   println("pos at $(integrator.t)")
        @show p[5]*((u[1]*p[2]-24.0)/(p[1]+p[2]))+(1.0-p[5])*((u[1]*p[2]-24.0)*p[1]/(p[1]+p[2])) =#
        discrete[6] = 1e-5
        discrete[2] = 1e5
           
end

function affect23!(integrator)
   #=  println("neg at $(integrator.t)")
    @show p[5]*((integrator.u[1]*p[2]-24.0)/(p[1]+p[2]))+(1.0-p[5])*((integrator.u[1]*p[2]-24.0)*p[1]/(p[1]+p[2])) =#
    discrete[7] = 1e-5
    discrete[3] = 1e5

end
function affect24!(integrator)
    #=  println("neg at $(integrator.t)")
     @show p[5]*((integrator.u[1]*p[2]-24.0)/(p[1]+p[2]))+(1.0-p[5])*((integrator.u[1]*p[2]-24.0)*p[1]/(p[1]+p[2])) =#
     discrete[8] = 1e-5
     discrete[4] = 1e5
 
 end
 function affect31!(integrator)
    #=  println("neg at $(integrator.t)")
     @show p[5]*((integrator.u[1]*p[2]-24.0)/(p[1]+p[2]))+(1.0-p[5])*((integrator.u[1]*p[2]-24.0)*p[1]/(p[1]+p[2])) =#
     discrete[5]  = 1e5
                discrete[1]  = 1e-5
 
 end
 function affect32!(integrator)
    #=  println("neg at $(integrator.t)")
     @show p[5]*((integrator.u[1]*p[2]-24.0)/(p[1]+p[2]))+(1.0-p[5])*((integrator.u[1]*p[2]-24.0)*p[1]/(p[1]+p[2])) =#
     discrete[6]  = 1e5
     discrete[2]  = 1e-5
 
 end
 function affect33!(integrator)
    #=  println("neg at $(integrator.t)")
     @show p[5]*((integrator.u[1]*p[2]-24.0)/(p[1]+p[2]))+(1.0-p[5])*((integrator.u[1]*p[2]-24.0)*p[1]/(p[1]+p[2])) =#
     discrete[7]  = 1e5
                discrete[3] = 1e-5
 
 end
 function affect34!(integrator)
    #=  println("neg at $(integrator.t)")
     @show p[5]*((integrator.u[1]*p[2]-24.0)/(p[1]+p[2]))+(1.0-p[5])*((integrator.u[1]*p[2]-24.0)*p[1]/(p[1]+p[2])) =#
     discrete[8] = 1e5
     discrete[4]  = 1e-5
 
 end
 function affect41!(integrator)
    #=  println("neg at $(integrator.t)")
     @show p[5]*((integrator.u[1]*p[2]-24.0)/(p[1]+p[2]))+(1.0-p[5])*((integrator.u[1]*p[2]-24.0)*p[1]/(p[1]+p[2])) =#
     discrete[1] = 1e5
 
 end
 function affect42!(integrator)
    #=  println("neg at $(integrator.t)")
     @show p[5]*((integrator.u[1]*p[2]-24.0)/(p[1]+p[2]))+(1.0-p[5])*((integrator.u[1]*p[2]-24.0)*p[1]/(p[1]+p[2])) =#
     discrete[2] = 1e5
 
 end
 function affect43!(integrator)
    #=  println("neg at $(integrator.t)")
     @show p[5]*((integrator.u[1]*p[2]-24.0)/(p[1]+p[2]))+(1.0-p[5])*((integrator.u[1]*p[2]-24.0)*p[1]/(p[1]+p[2])) =#
     discrete[3] = 1e5
 
 end
 function affect44!(integrator)
    #=  println("neg at $(integrator.t)")
     @show p[5]*((integrator.u[1]*p[2]-24.0)/(p[1]+p[2]))+(1.0-p[5])*((integrator.u[1]*p[2]-24.0)*p[1]/(p[1]+p[2])) =#
     discrete[4] = 1e5
 
 end
    cb1 = ContinuousCallback(condition11, affect11!,nothing; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 )
    cb2 = ContinuousCallback(condition21, affect21!,nothing; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 )
    cb3 = ContinuousCallback(condition22, affect22!,nothing; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 )
    cb4 = ContinuousCallback(condition23, affect23!,nothing; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 )
    cb5 = ContinuousCallback(condition24, affect24!,nothing; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 )
    cb6 = ContinuousCallback(condition31, affect31!,nothing; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 )
    cb7 = ContinuousCallback(condition32, affect32!,nothing; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 )
    cb8 = ContinuousCallback(condition33, affect33!,nothing; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 )
    cb9 = ContinuousCallback(condition34, affect34!,nothing; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 )
    cb10 = ContinuousCallback(condition41,nothing, affect41!; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 )
    cb11= ContinuousCallback(condition42,nothing, affect42!; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 )
    cb12= ContinuousCallback(condition43,nothing, affect43!; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 )
    cb13= ContinuousCallback(condition44,nothing, affect44!; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 )
    cbs= CallbackSet(cb1, cb2,cb3,cb4, cb5,cb6,cb7, cb8,cb9,cb10, cb11,cb12,cb13)
    
    tspan = (0.0, 0.000125)
   
    u0 = [0.0,0.0,0.0,0.0,0.0]
    discrete = [1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e-8,0.0]
    
    prob = ODEProblem(f, u0, tspan, discrete)
    sol = solve(prob, QBDF2(), callback = cbs, reltol=1e-3,abstol=1e-4#= dt = 1e-3, adaptive = false =#)
    
    p1=plot!(sol);
  # @show p[6]
   # p1=plot!(sol,ylims=(-2.04e-1,2.06e-1))

   # p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.0))
   #p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.06e-1))
   savefig(p1, "QBDF2_34_interleaved_ft00005_ogppos")
   #savefig(p1, "zoom1QBDF2_34_buck_ft0025")
  # p2=plot!(sol,ylims=(-2.04e-4,2.06e-4))
   #savefig(p2, "zoom2QBDF2_34_buck_ft0025")
 end
#@btime
 odeDiffEquPackage() 