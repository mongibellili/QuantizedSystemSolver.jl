using QuantizedSystemSolver
function test()
 
    odeprob = @NLodeProblem begin
          name=(RGElectrical,)
         ROn = 1e-5;ROff = 1e5;
          Lpr = 4.2*1e-9#0.48*1e-6#0.45 * 1e-6
L1 = 0.6*1e-4 #28.0*1e-6#0.6*1e-6
#L2 = 4.0e-6;L3 = 1.1*1e-6
L23=5.1e-6

 R1= 4.0e-3; R2 = 0.28*1e-3;R3 = 3.6*1e-3

C = 3.08*1e-3#3.08*1e-3

m=0.12
γ = 50.0*1e6; w = 15.0*1e-3#25.0*1e-3#15.0*1e-3 
μ = 4.0*3.14*1e-7
Rpr0 = 15.5*1e-6
FNmec = 680.0; α = 0.154

t0=1.0

          discrete = [1e5,1e-5,1.0,1e5,1e-5,1.0,1e5,1e-5,1.0,1.0,1.0,0.0,1e-3];u = [0.0,10000.75,0.0,0.0,4000.75,0.0,0.0,9000.75,0.0,0.0,0.0]
          rd1=discrete[1];rs1=discrete[2];charge1=discrete[3];
          rd2=discrete[4];rs2=discrete[5];charge2=discrete[6];
          rd3=discrete[7];rs3=discrete[8];charge3=discrete[9];
          start2=discrete[10]; start3=discrete[11]
          manualIntg=discrete[12]; nextT=discrete[13]
          is1=u[1] ;uc1=u[2]; il1=u[3] ;is2=u[4] ;uc2=u[5]; il2=u[6] ;is3=u[7] ;uc3=u[8]; il3=u[9] ;x=u[10]; v=u[11] 
          #α=LR+Lpr;
         # RR=4.0e-5
        rr=manualIntg*2.0*sqrt(μ/(3.14*γ))/w
         rpr=Rpr0*(sqrt(t0/(t+1e-4))+(t/t0)^16)/(1.0+(t/t0)^16)


          I=il1+il2+il3
          uf=0.1+0.2*exp(-v/100)
          F=0.5*0.453e-6*I*I*(1-uf*0.124)-uf*FNmec
          du[1] =((-(R1+rs1+rd1)*is1+rd1*il1+uc1)/L1)*charge1
          du[2]=(-is1/C)*charge1
         # du[3]=((-rd1-R2-R3)*il1-RR*(I)+rd1*is1)/(L2+L3)#((L2+L3+α)*β)/Δ

          du[3]=(-1.0108501920000001e-13il1 + 1.0108501920000001e-13il2 - 2.605284e-11il1*rd1 - 1.7927928e-14il1*x + 2.605284e-11il2*rd2 + 1.7927928e-14il2*x + 2.605284e-11is1*rd1 - 2.605284e-11is2*rd2 - 4.6206e-12il1*rd1*x + 4.6206e-12il2*rd2*x + 4.6206e-12is1*rd1*x - 4.6206e-12is2*rd2*x)/(1.32869484e-16 + 2.356506e-17x - 3.85185988e-34(x^2))

          du[4] =((-(R1+rs2+rd2)*is2+rd2*il2+uc2)/L1)*charge2*start2
          du[5]=(-is2/C)*charge2*start2
          #du[6]=start2*(((-rd2-R2-R3)*il2-RR*(I)+rd2*is2)/(L2+L3))#((L2+L3+α)*β)/Δ
          du[6]=(8.31096e-17(il1 + il3) - 1.0108501920000001e-13il2 - 2.6010000000000002e-11I*(rpr + rr) - 1.1782530000000001e-17I*v + 2.142e-14il1*rd1 + 8.963964e-15il1*x - 2.605284e-11il2*rd2 - 1.7927928e-14il2*x + 2.142e-14il3*rd3 + 8.963964e-15il3*x - 2.142e-14is1*rd1 + 2.605284e-11is2*rd2 - 2.142e-14is3*rd3 + 2.3103e-12il1*rd1*x - 4.6206e-12il2*rd2*x + 2.3103e-12il3*rd3*x - 2.3103e-12is1*rd1*x + 4.6206e-12is2*rd2*x - 2.3103e-12is3*rd3*x)/(1.32869484e-16 + 2.356506e-17x - 3.85185988e-34(x^2))

          du[7] =((-(R1+rs3+rd3)*is3+rd3*il3+uc3)/L1)*charge3*start3
          du[8]=(-is3/C)*charge3*start3
         # du[9]=start3*(((-rd3-R2-R3)*il3-RR*(I)+rd3*is3)/(L2+L3))#((L2+L3+α)*β)/Δ
          du[9]=(8.31096e-17il1 - 1.0100191154e-13il3 - 2.60100005e-11I*(rpr + rr) - 1.17825302265e-17I*v + 2.142e-14il1*rd1 + 8.963964e-15il1*x - 2.60314205e-11il3*rd3 - 8.963964e-15il3*x - 2.142e-14is1*rd1 + 2.60314205e-11is3*rd3 + 2.3103e-12il1*rd1*x - 2.3103e-12il3*rd3*x - 2.3103e-12is1*rd1*x + 2.3103e-12is3*rd3*x)/(1.32869484e-16 + 2.356506e-17x - 3.85185988e-34(x^2))



          du[10]=v
          du[11]=F/m

        #=   if t-0.0009>0.0 
            start3=1.0
          end

          if t-0.0016>0.0 
            start2=1.0
           # RR=4.0e-7
          end =#
        
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
          if -(is1)>0.0 
            rs1=ROff;charge1=0.0 # rs off not needed since charge=0
            rd1=ROn;
          #=   uc1=0.0;
            is1=0.0 =#
          end 
          if -(is2)>0.0 
            rs2=ROff;charge2=0.0
            rd2=ROn;
            @show rd2
            
          end 
         
          if -(is3)>0.0 
            rs3=ROff;charge3=0.0
            rd3=ROn;
           
          end 
    
          if t-nextT>0
            manualIntg=manualIntg+v/sqrt(t[0])
            nextT=nextT+1e-3
          end
        
         #=  if -il3>0
            sep_rd3=1e3
          end =#

          
    end
    #@show odeprob
    tspan = (0.0, 0.001)
    sol= solve(odeprob,nmliqss2(),tspan,abstol=1e-4,reltol=1e-3)    
  #=  save_Sol(sol,1)
   save_Sol(sol,2)
   save_Sol(sol,3) 
   save_Sol(sol,4)
  save_Sol(sol,5)
  save_Sol(sol,6) 
   save_Sol(sol,7) 
   save_Sol(sol,8)  =#
  save_Sol(sol,9,xlims=(0.0,4e-5)) 
#=    save_SolSum(sol,3,6,9) #add interpl preference
   save_Sol(sol,10) 
  save_Sol(sol,11)  =#
   
end
test()
