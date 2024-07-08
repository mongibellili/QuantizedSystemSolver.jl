using OrdinaryDiffEq
#= using BenchmarkTools
using XLSX =#
using Plots



function odeDiffEquPackage()
 
    function funcName(du,u,p,t)# api requires four args
       #=  du[1] = acos(sin(u[2]))
        du[2] = (u[1]) =#
        #= du[1] = t
        du[2] =1.24*u[1]-0.01*u[2]+0.2 =#
       
        du[1] = t
        
        for k in 2:5 
            du[k]=(u[k]-u[k-1]) ;
        end 
    end
    tspan = (0.0,5.0)
    u0= [1.0, 0.0,1.0, 0.0,1.0]
    prob = ODEProblem(funcName,u0,tspan)


   absTol=1e-5
   relTol=1e-2


 solRosenbrock23 = solve(prob,Rosenbrock23(),saveat=0.01,abstol = absTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
 p1=plot!(solRosenbrock23,marker=(:circle),markersize=2)
 savefig(p1, "plot_solRosenbrock23_N8.png")
 
end

odeDiffEquPackage()  



