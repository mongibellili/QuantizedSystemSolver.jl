using qss
using BenchmarkTools


function testCuk4()
    odeprob = @NLodeProblem begin
         name=(cuk4sym,)
         #parameters:
         C = 1e-4; L = 1e-4; R = 10.0;U = 24.0; T = 1e-4; DC = 0.25; ROn = 1e-5;ROff = 1e5;L1=1e-4;C1=1e-4;C2 = 1e-4;L2 = 1e-4;
        #initial conditions:
         u[1:13]=0.0
        #ODEs:
        for i=1:4    #il2
          du[i] =(-uc2-Rs*id1)/L2
        end
        for i=5:8    #il1
          du[i]=(U-uc1_2-id2*Rs)/L1
        end
        for i=9:12    #uc1
          du[i]=(id3-il2_3)/C1
        end
        du[13]=(u[1]+u[2]+u[3]+u[4]-uc2/R)/C2
        #events
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
  
    tspan=(0.0,0.0005)
   sol= solve(odeprob,nmliqss2(),abstol=1e-4,reltol=1e-3,tspan)
   save_Sol(sol)          
 #=  save_Sol(sol,1,note="1")
  save_Sol(sol,5,note="5")
  save_Sol(sol,9,note="9")
  save_Sol(sol,13,note="13") =#
# save_Sol(sol,xlims=(9.226846540089e-5,9.25e-5) ,ylims=(8.118295828650566,8.14))

#  save_Sol(sol,13,xlims=(0.000093,0.000098) ,ylims=(-19.0,-16.0))


 # save_Sol(sol,8,xlims=(0.000093,0.0000935) ,ylims=(8.13,8.15))
end
#@btime 
test()