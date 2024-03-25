using QuantizedSystemSolver
function test()
 
    odeprob = @NLodeProblem begin
          name=(RGElectrical,)
         ROn = 1e-5;ROff = 1e5;
         # Lpr = 4.2*1e-9#0.48*1e-6#0.45 * 1e-6
L1 = 0.6*1e-4 #28.0*1e-6#0.6*1e-6
L2 = 4.0e-6;L3 = 1.1*1e-6

 R1= 4.0e-3; R2 = 0.28*1e-3;R3 = 3.6*1e-3

C = 3.08*1e-3#3.08*1e-3



          discrete = [1e5,1e-5,1.0,1e5,1e-5,1.0,1e5,1e-5,1.0,1e5,1e-5,1.0,0.0,0.0,0.0];u = [0.0,10000.75,0.0,0.0,10000.75,0.0,0.0,10000.75,0.0,0.0,10000.75,0.0]
          rd1=discrete[1];rs1=discrete[2];charge1=discrete[3];
          rd2=discrete[4];rs2=discrete[5];charge2=discrete[6];
          rd3=discrete[7];rs3=discrete[8];charge3=discrete[9];
          rd4=discrete[10];rs4=discrete[11];charge4=discrete[12];
          start2=discrete[13]; start3=discrete[14]
       
          is1=u[1] ;uc1=u[2]; il1=u[3] ;is2=u[4] ;uc2=u[5]; il2=u[6] ;is3=u[7] ;uc3=u[8]; il3=u[9] ;is4=u[10] ;uc4=u[11]; il4=u[12] 
          #α=LR+Lpr;
          RR=4.0e-5
        
          du[1] =((-(R1+rs1+rd1)*is1+rd1*il1+uc1)/L1)*charge1
          du[2]=(-is1/C)*charge1
          du[3]=((-rd1-R2-R3-RR)*il1-RR*(il2+il3+il4)+rd1*is1)/(L2+L3)#((L2+L3+α)*β)/Δ

          du[4] =((-(R1+rs2+rd2)*is2+rd2*il2+uc2)/L1)*charge2
          du[5]=(-is2/C)*charge2
          du[6]=(((-rd2-R2-R3-RR)*il2-RR*(il1+il3+il4)+rd2*is2)/(L2+L3))#((L2+L3+α)*β)/Δ

         
          du[7] =((-(R1+rs3+rd3)*is3+rd3*il3+uc3)/L1)*charge3*start3
          du[8]=(-is3/C)*charge3*start3
          du[9]=start3*(((-rd3-R2-R3-RR)*il3-RR*(il1+il2+il4)+rd3*is3)/(L2+L3))#((L2+L3+α)*β)/Δ


          du[10] =(((-(R1+rs4+rd4)*is4+rd4*il4+uc4)/L1)*charge4)*start2
          du[11]=((-is4/C)*charge4)*start2
          du[12]=(((-rd4-R2-R3-RR)*il4-RR*(il1+il2+il3)+rd4*is4)/(L2+L3))*start2#((L2+L3+α)*β)/Δ

          if t-0.0009>0.0 
            start3=1.0
          end

          if t-0.0016>0.0 
            start2=1.0
           # RR=4.0e-7
          end
        
          if -(uc1)>0.0 
            rs1=ROff;charge1=0.0 # rs off not needed since charge=0
            rd1=ROn;
          #=   uc1=0.0;
            is1=0.0 =#
          end 
          if -(uc2)>0.0 
            rs2=ROff;charge2=0.0
            rd2=ROn;
            @show rd2
            
          end 
         
          if -(uc3)>0.0 
            rs3=ROff;charge3=0.0
            rd3=ROn;
           
          end 
    
          if -(uc4)>0.0 
            rs4=ROff;charge4=0.0
            rd4=ROn;
           
          end 
        
         #=  if -il3>0
            sep_rd3=1e3
          end =#

          
    end
    #@show odeprob
    tspan = (0.0, 0.006)
    sol= solve(odeprob,nmliqss2(),tspan,abstol=1e-4,reltol=1e-3)    
  #=  save_Sol(sol,1)
   save_Sol(sol,2)
   save_Sol(sol,3) 
   save_Sol(sol,4)
  save_Sol(sol,5)
  save_Sol(sol,6) 
   save_Sol(sol,7) 
   save_Sol(sol,8) 
  save_Sol(sol,9)  =#
  save_Sol(sol,10) 
  save_Sol(sol,11) 
 save_Sol(sol,12) 
   save_SolSum(sol,3,6,9,12) #add interpl preference
   
end
test()
