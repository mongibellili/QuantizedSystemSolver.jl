using qss
using BenchmarkTools


function test()
    odeprob = @NLodeProblem begin
      name=(buck,)
         C = 1e-4; L = 4e-4; R = 30;U = 50.0; T = 1e-4; DC = 0.5; ROn = 1e-5;ROff = 1e5;
         #discrete Rd(start=1e5), Rs(start=1e-5), nextT(start=T),lastT,diodeon;
         discrete = [1e5,1e-5,1e-4,0.0,0.0]
       
        u = [0.0,0.0]
      
        du[1] =(-((u[1]*discrete[2]-24.0)*discrete[1]/(discrete[1]+discrete[2]))-u[2])/L
        du[2]=(u[1]-u[2]/R)/C

        if t-discrete[3]>0.0 
            discrete[4]=discrete[3]
            discrete[3]=discrete[3]+T
            discrete[2]=ROn
        end

        if t-discrete[4]-DC*T>0.0 
            discrete[2]=ROff
        end                          
        
      if discrete[5]*((u[1]*discrete[2]-U)/(discrete[1]+discrete[2]))+(1.0-discrete[5])*((u[1]*discrete[2]-U)*discrete[1]/(discrete[1]+discrete[2]))>0
        discrete[1]=ROn
        discrete[5]=1.0
      else
        discrete[1]=ROff
        discrete[5]=0.0
      end 

   
           
    end
    tspan=(0.0,0.01)
   sol= solve(odeprob,nmliqss2(),abstol=1e-4,reltol=1e-3,tspan)
            
  save_Sol(sol)
 # save_Sol(sol,xlims=(0.0,0.0006) ,ylims=(-2.04e-1,40.0))
end
#@btime 
test()
