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



          discrete = [1e5,1e-5,1.0,1e5,1e-5,1.0];u = [0.0,7990.75,0.0,0.0,10000.75,0.0]
          rd1=discrete[1];rs1=discrete[2];charge1=discrete[3];rd2=discrete[4];rs2=discrete[5];charge2=discrete[6]#LR = discrete[3];nextT=discrete[4];manualIntg=discrete[5]
          is1=u[1] ;uc1=u[2]; il1=u[3] ;is2=u[4] ;uc2=u[5]; il2=u[6] 
          #α=LR+Lpr;
          RR=4.0e-5
        
          du[1] =((-(R1+rs1+rd1)*is1+rd1*il1+uc1)/L1)*charge1
          du[2]=(-is1/C)*charge1
          du[3]=((-rd1-R2-R3-RR)*il1-RR*il2+rd1*is1)/(L2+L3)#((L2+L3+α)*β)/Δ
          du[4] =((-(R1+rs2+rd2)*is2+rd2*il2+uc2)/L1)*charge2
          du[5]=(-is2/C)*charge2
          du[6]=((-rd2-R2-R3-RR)*il2-RR*il1+rd2*is2)/(L2+L3)#((L2+L3+α)*β)/Δ

          #= if t-0.00025>0.0 
            rs1=ROn
           
          end =#
          if -uc1*is1>0.0 
            rs1=ROff;charge1=0.0
            rd1=ROn;
          end 
          if -uc2*is2>0.0 
            rs2=ROff;charge2=0.0
            rd2=ROn;
          end 
    
       #=    if ((il1-is1)*rd1)-0.7>0
            rd1=ROn;
          else
            rd1=ROff;
          end    =#  

        #=   if diodeon*((il1-is1))+(1.0-diodeon)*((il1-is1)*rd1-0.7)>0
            rd1=ROn;diodeon=1.0
         #=  else
            rd1=ROff;diodeon=0.0 =#
          end      =#



          
    end
    #@show odeprob
    tspan = (0.0, 0.006)
    sol= solve(odeprob,nmliqss2(),tspan,abstol=1e-4,reltol=1e-3)    
    #= save_Sol(sol,1)
    save_Sol(sol,2)
    save_Sol(sol,3)
    save_Sol(sol,4)
    save_Sol(sol,5)
    save_Sol(sol,6) =#
    save_SolSum(sol,3,6)
   
end
test()
