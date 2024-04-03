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



          discrete = [1e5,1e-5,1.0,1e5,1e-5,1.0,1e5,1e-5,1.0,0.0];u = [0.0,9990.75,0.0,0.0,8000.75,0.0,0.0,7000.75,0.0]
          rd1=discrete[1];rs1=discrete[2];charge1=discrete[3];
          rd2=discrete[4];rs2=discrete[5];charge2=discrete[6];
          rd3=discrete[7];rs3=discrete[8];charge3=discrete[9];
          start=discrete[10]#LR = discrete[3];nextT=discrete[4];manualIntg=discrete[5]
          is1=u[1] ;uc1=u[2]; il1=u[3] ;is2=u[4] ;uc2=u[5]; il2=u[6] ;is3=u[7] ;uc3=u[8]; il3=u[9] 
          #α=LR+Lpr;
          RR=4.0e-5
        
      


          rd=discrete[3*i-2];rs=discrete[3*i-1];charge=discrete[3*i]
          is=u[3*i-2] ;uc=u[3*i-1] ;il=u[3*i] 

          I=u[3] +u[6] +u[9] 
          for i in 1:3
            du[3*i-2] =((-(R1+rs+rd)*is+rd*il+uc)/L1)*charge
            du[3*i-1]=(-is/C)*charge
            du[3*i]=((-rd-R2-R3-RR)*il-RR*(I-il)+rd*is)/(L2+L3)#((L2+L3+α)*β)/Δ
          end 




          if t-0.001>0.0 
            start=1.0
           
          end
          if -uc1*is1>0.0 
            rs1=ROff;charge1=0.0
            rd1=ROn;
          end 
          if -uc2*is2>0.0 
            rs2=ROff;charge2=0.0
            rd2=ROn;
          end 
          if -uc3*is3>0.0 
            rs3=ROff;charge3=0.0
            rd3=ROn;
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
   #=  save_Sol(sol,1)
    save_Sol(sol,2)
    save_Sol(sol,3)
    save_Sol(sol,4)
    save_Sol(sol,5)
    save_Sol(sol,6) =#
    save_Sol(sol,3) 
    save_Sol(sol,6) 
    save_Sol(sol,9) 
    save_SolSum(sol,3,6,9) #add interpl preference
   
end
test()
