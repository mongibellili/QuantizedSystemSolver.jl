



using DifferentialEquations
#using qssv01
using BenchmarkTools
using Plots
function odeDiffEquPackage()
    function f(du, u, p, t)
        C = 1e-4; L = 1e-4; R = 10;U = 24.0; T = 1e-4; DC = 0.25; ROn = 1e-5;ROff = 1e5;L1=1e-4;C1=1e-4
        du[1] =(-(((u[1]+u[2])*p[2]-u[4])*p[1]/(p[1]+p[2]))-u[3])/L             #der(iL) =  (-uC-iD*Rd)/L;
        du[2]=(U-u[4]-(((u[1]+u[2])*p[2]-u[4])*p[1]/(p[1]+p[2])))/L1            #der(iL1) =  (U-uC1-iD*Rd)/L1;
        du[3]=(u[1]-u[3]/R)/C                                                   #der(uC) = (iL - uC/R)/C;
        du[4]=(((u[1]+u[2])*p[2]-u[4])/(p[1]+p[2])-u[1])/C1                     #der(uC1) = (iD - iL)/C1;
    end
    function condition1( u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        (t-p[3])
    end
    function condition2( u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        (t-p[4]-0.5*1e-4)
    end
    function condition3( u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        p[5]*(((u[1]+u[2])*p[2]-u[4])/(p[1]+p[2]))+(1.0-p[5])*(((u[1]+u[2])*p[2]-u[4])*p[1]/(p[1]+p[2]))
    end
    function affect1!(integrator)
       #=  println("pos at $(integrator.t)")
        @show integrator.t-p[3] =#
                p[4]=p[3]
                p[3]=p[3]+1e-4
                p[2]=1e-5
                p[6]+=1.0
         
    end
    function affect2!(integrator)
        #= println("pos at $(integrator.t)")
        @show integrator.t-p[4]-0.5*1e-4 =#
                p[2]=1e5
                p[6]+=1.0
             
    end
    function affect3!(integrator)
      #=   println("pos at $(integrator.t)")
        @show p[5]*((u[1]*p[2]-24.0)/(p[1]+p[2]))+(1.0-p[5])*((u[1]*p[2]-24.0)*p[1]/(p[1]+p[2])) =#
            p[1]=1e-5
            p[5]=1.0
            p[6]+=1.0
           
end

function affect6!(integrator)
   #=  println("neg at $(integrator.t)")
    @show p[5]*((integrator.u[1]*p[2]-24.0)/(p[1]+p[2]))+(1.0-p[5])*((integrator.u[1]*p[2]-24.0)*p[1]/(p[1]+p[2])) =#
    p[1]=1e5
    p[5]=0.0
    p[6]+=1.0

end

cb1 = ContinuousCallback(condition1, affect1!,nothing; )
cb2 = ContinuousCallback(condition2, affect2!,nothing;  )
cb3 = ContinuousCallback(condition3, affect3!,affect6!;  )

  #=   cb1 = ContinuousCallback(condition1, affect1!,nothing; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 )
    cb2 = ContinuousCallback(condition2, affect2!,nothing; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 )
    cb3 = ContinuousCallback(condition3, affect3!,affect6!; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 ) =#
    #cb4 = ContinuousCallback(condition3,nothing, affect6!; rootfind=SciMLBase.RightRootFind, save_positions=(true,true),interp_points=100,abstol=10eps(),reltol=0,repeat_nudge=1//100 )
    cbs = CallbackSet(cb1, cb2,cb3)
    u0 = [0.0, 0.0,0.0, 0.0]
    tspan = (0.0, 0.005)
    p = [1e5,1e-5,1e-4,0.0,0.0,0.0]
      
    
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob,Tsit5() #= QBDF2() =#, callback = cbs, reltol=1e-3,abstol=1e-4#= dt = 1e-3, adaptive = false =#)
    
    p1=plot!(sol);
   @show p[6]
   # p1=plot!(sol,ylims=(-2.04e-1,2.06e-1))

   # p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.0))
   #p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.06e-1))
   savefig(p1, "Tsit5()_34_cuk_ft005_")
   #savefig(p1, "zoom1QBDF2_34_buck_ft0025")
  # p2=plot!(sol,ylims=(-2.04e-4,2.06e-4))
   #savefig(p2, "zoom2QBDF2_34_buck_ft0025")
 end
#@btime
 odeDiffEquPackage() 