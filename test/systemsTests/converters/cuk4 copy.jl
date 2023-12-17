using qss
using BenchmarkTools


function test()
    odeprob = @NLodeProblem begin
      name=(cuk4_,)
         C = 1e-4; L = 1e-4; R = 10.0;U = 24.0; T = 1e-4; DC = 0.25; ROn = 1e-5;ROff = 1e5;L1=1e-4;C1=1e-4;C2 = 1e-4;L2 = 1e-4;
         #discrete Rd(start=1e5), Rs(start=1e-5), nextT(start=T),lastT,diodeon;
         discrete = [1e5,1e-5,1e-4,0.0,0.0]
         Rd=discrete[1];Rs=discrete[2];nextT=discrete[3];lastT=discrete[4];diodeon=discrete[5]
         u[1:13]=0.0
        for i=1:4
          du[i] =(-Rs*(((u[i]+u[i+4])*Rs-u[i+8])/(Rd+Rs))-u[13])/L2
        end
        for i=5:8
          du[i]=(U-u[i+4]-Rs*(((u[i]+u[i-4])*Rs-u[i+4])/(Rd+Rs)))/L1
        end
        for i=9:12
          du[i]=(((u[i-4]+u[i-8])*Rs-u[i])/(Rd+Rs)-u[i-8])/C1
        end

        du[13]=(u[1]+u[2]+u[3]+u[4]-u[13]/R)/C2



        if t-nextT>0.0 
            lastT=nextT
            nextT=nextT+T
            Rs=ROn
        end

        if t-lastT-DC*T>0.0 
            Rs=ROff
        end                          
        
      if diodeon*(((u[1]+u[5])*Rs-u[9])/(Rd+Rs))+(1.0-diodeon)*(((u[1]+u[5])*Rs-u[9])*Rd/(Rd+Rs))>0
        Rd=ROn
        diodeon=1.0
      else
        Rd=ROff
        diodeon=0.0
      end 

      if diodeon*(((u[2]+u[6])*Rs-u[10])/(Rd+Rs))+(1.0-diodeon)*(((u[2]+u[6])*Rs-u[10])*Rd/(Rd+Rs))>0
        Rd=ROn
        diodeon=1.0
      else
        Rd=ROff
        diodeon=0.0
      end 
           

      if diodeon*(((u[3]+u[7])*Rs-u[11])/(Rd+Rs))+(1.0-diodeon)*(((u[3]+u[7])*Rs-u[11])*Rd/(Rd+Rs))>0
        Rd=ROn
        diodeon=1.0
      else
        Rd=ROff
        diodeon=0.0
      end 

      if diodeon*(((u[4]+u[8])*Rs-u[12])/(Rd+Rs))+(1.0-diodeon)*(((u[4]+u[8])*Rs-u[12])*Rd/(Rd+Rs))>0
        Rd=ROn
        diodeon=1.0
      else
        Rd=ROff
        diodeon=0.0
      end 


    end
    #@show odeprob
   sol= solve(odeprob,nmliqss2(),abstol=1e-4,reltol=1e-3,finalTime=0.00002501)
            
  save_Sol(sol,1,5,9,13,note="159_13")
 #=  save_Sol(sol,2,6,10,note="26_10")
  save_Sol(sol,3,7,11,note="37_11")
  save_Sol(sol,4,8,12,note="48_12") =#
 # save_Sol(sol,1,ylims=(-5.0,-4.5))
 # save_Sol(sol,1,ylims=(-4.9,-4.85))
 # save_Sol(sol,1,ylims=(-4.9,-4.895))
end
#@btime 
test()
