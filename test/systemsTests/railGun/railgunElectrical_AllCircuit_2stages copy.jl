using QuantizedSystemSolver
using BenchmarkTools
function testRailGun()
    odeprob = @NLodeProblem begin
          name=(RG2Cap,)
          #parameters
          ROn = 1e-5;ROff = 1e0; Lpr = 4.2*1e-9 ;  L1 = 0.6*1e-4 #L2 = 4.0e-6;L3 = 1.1*1e-6 #L23=5.1e-6
          R1= 4.0e-3; R2 = 0.28*1e-3;R3 = 3.6*1e-3;C = 3.08*1e-3;m=0.12
          γ = 50.0*1e6; w = 15.0*1e0; μ = 4.0*3.14*1e-7;Rpr0 = 1.5*1e-6; FNmec = 680.0; α = 0.154; t0=1.0;
          #continuous and discrete vars
          discrete = [1e5,1e-5,1.0,1e5,1e-5,1.0,1.0,0.0,1e-3];u = [0.0,10000.75,0.0,0.0,10000.75,0.0,0.0,0.0]
          rd1=discrete[1];rs1=discrete[2];charge1=discrete[3]; rd2=discrete[4];rs2=discrete[5];charge2=discrete[6];start2=discrete[7];  manualIntg=discrete[8]; nextT=discrete[9]
          is1=u[1] ;uc1=u[2]; il1=u[3] ;is2=u[4] ;uc2=u[5]; il2=u[6] ;x=u[7]; v=u[8] 
          #helper functions
          rr=manualIntg*2.0*sqrt(μ/(3.14*γ))/w
          rpr=Rpr0*(sqrt(t0/(t+1e-4))+(t/t0)^16)/(1.0+(t/t0)^16)
          I=il1+il2
          uf=0.1+0.2*exp(-v/100)
          F=0.5*0.453e-6*I*I*(1-uf*0.124)-uf*FNmec
          #differential equations
          du[1] =((-(R1+rs1+rd1)*is1+rd1*il1+uc1)/L1)#*charge1
          du[2]=(-is1/C)*charge1
          du[3]= ((-I*(rpr + rr + 4.53e-7v) + il2*(-0.00388 - rd2) + is2*rd2)*(-4.2e-9 - 4.53e-7x)) / (2.6052840000000004e-11 + 4.6206e-12x) + ((-I*(rpr + rr + 4.53e-7v) + il1*(-0.00388 - rd1) + is1*rd1)*(5.1042e-6 + 4.53e-7x)) / (2.6052840000000004e-11 + 4.6206e-12x)
          du[4] =((-(R1+rs2+rd2)*is2+rd2*il2+uc2)/L1)#*charge2*start2
          du[5]=(-is2/C)*charge2*start2
          du[6]=((-I*(rpr + rr + 4.53e-7v) + il1*(-0.00388 - rd1) + is1*rd1)*(-4.2e-9 - 4.53e-7x)) / (2.6052840000000004e-11 + 4.6206e-12x) + ((-I*(rpr + rr + 4.53e-7v) + il2*(-0.00388 - rd2) + is2*rd2)*(5.1042e-6 + 4.53e-7x)) / (2.6052840000000004e-11 + 4.6206e-12x)
          du[7]=v
          du[8]=F/m
          #events
          if -(uc1)>0.0 
            rs1=ROff;charge1=0.0 # rs off not needed since charge=0
            rd1=ROn;
          end 
          if -(uc2)>0.0 
            rs2=ROff;charge2=0.0
            rd2=ROn;
          end 
         #=  if -(il1)>0.0 
            
          end 
          if -(il2)>0.0 
           
          end  =#
         #=  if rs1*((il1-is1))+(1.0-rs1)*((il1-is1)*rd1-0.6)>0 
            # rs off not needed since charge=0
            rd1=ROn;rs1=1.0
          else
            rd1=ROff;rs1=0.0
          end 
          if rs2*((il2-is2))+(1.0-rs2)*((il2-is2)*rd2-0.6)>0
            # rs off not needed since charge=0
            rd2=ROn;rs2=1.0
          else
            rd2=ROff;rs2=0.0
          end  =#
          if t-nextT>0
            manualIntg=manualIntg+v/sqrt(t[0])
            nextT=nextT+1e-3
          end
    end
    #solve
    tspan = (0.0, 0.006)
    sol= solve(odeprob,nmliqss1(),tspan,abstol=1e-3,reltol=1e-2)  
    #save  
   save_Sol(sol,1)
   save_Sol(sol,2)
   save_Sol(sol,3) 
   save_Sol(sol,4)
   save_Sol(sol,5)
   save_Sol(sol,6)
   save_Sol(sol,7) 
   save_Sol(sol,8) 
   save_SolSum(sol,3,6)  
end
#@btime 
testRailGun()
