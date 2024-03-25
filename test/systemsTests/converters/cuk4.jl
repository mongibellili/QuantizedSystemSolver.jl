using QuantizedSystemSolver
#using BenchmarkTools


function test()
    odeprob = @NLodeProblem begin
      name=(cuk4,)
         C = 1e-4; L = 1e-4; R = 10.0;U = 24.0; T = 1e-4; DC = 0.25; ROn = 1e-5;ROff = 1e5;L1=1e-4;C1=1e-4;C2 = 1e-4;L2 = 1e-4;
         #discrete Rd(start=1e5), Rs(start=1e-5), nextT(start=T),lastT,diodeon;
         discrete = [1e5,1e-5,1e-4,0.0,0.0]
       
         u[1:13]=0.0
        for i=1:4
          du[i] =(-(((u[i]+u[i+4])*discrete[2]-u[i+8])*discrete[2]/(discrete[1]+discrete[2]))-u[13])/L2
        end
        for i=5:8
          du[i]=(U-u[i+4]-(((u[i]+u[i-4])*discrete[2]-u[i+4])*discrete[2]/(discrete[1]+discrete[2])))/L1
        end
        for i=9:12
          du[i]=(((u[i-4]+u[i-8])*discrete[2]-u[i])/(discrete[1]+discrete[2])-u[i-8])/C1
        end

        du[13]=(u[1]+u[2]+u[3]+u[4]-u[13]/R)/C2



        if t-discrete[3]>0.0 
            discrete[4]=discrete[3]
            discrete[3]=discrete[3]+T
            discrete[2]=ROn
        end

        if t-discrete[4]-DC*T>0.0 
            discrete[2]=ROff
        end                          
        
      if discrete[5]*(((u[1]+u[5])*discrete[2]-u[9])/(discrete[1]+discrete[2]))+(1.0-discrete[5])*(((u[1]+u[5])*discrete[2]-u[9])*discrete[1]/(discrete[1]+discrete[2]))>0
        discrete[1]=ROn
        discrete[5]=1.0
      else
        discrete[1]=ROff
        discrete[5]=0.0
      end 

      if discrete[5]*(((u[2]+u[6])*discrete[2]-u[10])/(discrete[1]+discrete[2]))+(1.0-discrete[5])*(((u[2]+u[6])*discrete[2]-u[10])*discrete[1]/(discrete[1]+discrete[2]))>0
        discrete[1]=ROn
        discrete[5]=1.0
      else
        discrete[1]=ROff
        discrete[5]=0.0
      end 
           

      if discrete[5]*(((u[3]+u[7])*discrete[2]-u[11])/(discrete[1]+discrete[2]))+(1.0-discrete[5])*(((u[3]+u[7])*discrete[2]-u[11])*discrete[1]/(discrete[1]+discrete[2]))>0
        discrete[1]=ROn
        discrete[5]=1.0
      else
        discrete[1]=ROff
        discrete[5]=0.0
      end 

      if discrete[5]*(((u[4]+u[8])*discrete[2]-u[12])/(discrete[1]+discrete[2]))+(1.0-discrete[5])*(((u[4]+u[8])*discrete[2]-u[12])*discrete[1]/(discrete[1]+discrete[2]))>0
        discrete[1]=ROn
        discrete[5]=1.0
      else
        discrete[1]=ROff
        discrete[5]=0.0
      end 


    end
    tspan=(0.0,0.001)
   sol= solve(odeprob,nmliqss2(),abstol=1e-4,reltol=1e-3,tspan)
            
  save_Sol(sol)
  #save_Sol(sol,xlims=(1.9e-2,2e-2) ,ylims=(-6.0,8))
end
#@btime 
test()
