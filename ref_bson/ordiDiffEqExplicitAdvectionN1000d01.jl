using OrdinaryDiffEq

function odeDiffEquPackage()
 
    function funcName(du,u,p,t)# api requires four args
      _dx=100.0#1/dx=N/10=1000/10
      a=1.0
      d=0.1
      r=1000.0
      du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
      for k in 2:999  
          du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
      end 
      du[1000]=-a*_dx*(u[1000]-u[999])+d*_dx*_dx*(2.0*u[999]-2.0*u[1000])+r*u[1000]*u[1000]*(1.0-u[1000]) 
    
    end
    tspan = (0.0,10.0)
    
    u0=zeros(1000)
    u0[1:333].=1.0


    prob = ODEProblem(funcName,u0,tspan)

    solFeagin14 = solve(prob,Feagin14(),saveat=0.01,abstol = 1e-12, reltol = 1e-8)

     solFeagin14VectorN1000d01=solFeagin14.u

     BSON.@save "formalqss/ref_bson/solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01
 
end
#@time 
odeDiffEquPackage()  


