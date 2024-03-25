using OrdinaryDiffEq
using BSON



function odeDiffEquPackage()
 
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[4]-1e6*u[1]+1e3*u[2]
        du[2] =-200.0*u[2]*u[5]+1e6*u[1]-1e3*u[2]
        du[3] = 200.0*u[2]*u[5]-u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)
        du[4] =u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)-u[4]
        du[5] = 0.015-200.0*u[2]*u[5]
        du[6] =u[4]-0.6*u[6]
    end
    tspan = (0.0,100.0)
    u0= [0.0,0.75,0.25,0.0,0.0,0.0]
    prob = ODEProblem(funcName,u0,tspan)


    solRodas5P = solve(prob,Rodas5P(),saveat=0.01,abstol = 1e-12, reltol = 1e-8)

    solRodas5PVectorTyson=solRodas5P.u

   
  
   BSON.@save "formalqss/ref_bson/solVectAdvection_Tyson_Rodas5Pe-12.bson" solRodas5PVectorTyson


  
end
#@btime 
odeDiffEquPackage()  



