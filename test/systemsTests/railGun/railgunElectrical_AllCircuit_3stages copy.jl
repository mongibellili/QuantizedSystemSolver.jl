using QuantizedSystemSolver
#using BenchmarkTools
function test()
 
  odeprob = NLodeProblem(
    quote
          name=(RGElectrical3,)
         ROn = 1e-5;ROff = 1e1;
          Lpr = 4.2*1e-9#0.48*1e-6#0.45 * 1e-6
L1 = 0.6*1e-4 #28.0*1e-6#0.6*1e-6
#L2 = 4.0e-6;L3 = 1.1*1e-6
L23=5.1e-6

 R1= 4.0e-3; R2 = 0.28*1e-3;R3 = 3.6*1e-3

C = 3.08*1e-3#3.08*1e-3

m=0.12
γ = 50.0*1e6; w = 15.0*1e0#25.0*1e-3#15.0*1e-3 
μ = 4.0*3.14*1e-7
Rpr0 = 1.5*1e-6
FNmec = 680.0; α = 0.154
rs11=1e-5;rs21=1e-5;rs31=1e-5;
t0=1.0

          discrete = [1e5,1e-5,1.0,1e5,1e-5,1.0,1e5,1e-5,1.0,1.0,1.0,0.0,1e-3,1.0,1.0,1.0];u = [0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,0.0]
          rd1=discrete[1];rs1=discrete[2];charge1=discrete[3];
          rd2=discrete[4];rs2=discrete[5];charge2=discrete[6];
          rd3=discrete[7];rs3=discrete[8];charge3=discrete[9];
          start2=discrete[10]; start3=discrete[11]
          manualIntg=discrete[12]; nextT=discrete[13]
          operate1=discrete[14]; operate2=discrete[15];operate3=discrete[16]; 
          is1=u[1] ;uc1=u[2]; il1=u[3] ;is2=u[4] ;uc2=u[5]; il2=u[6] ;is3=u[7] ;uc3=u[8]; il3=u[9] ;x=u[10]; v=u[11] 
          #α=LR+Lpr;
         # RR=4.0e-5
        rr=manualIntg*2.0*sqrt(μ/(3.14*γ))/w
         rpr=Rpr0*(sqrt(t0/(t+1e-4))+(t/t0)^16)/(1.0+(t/t0)^16)


          I=il1+il2+il3
          uf=0.1+0.2*exp(-v/100)
          F=0.5*0.453e-6*I*I*(1-uf*0.124)-uf*FNmec
          du[1] =((-(R1+rs1+rd1)*is1+rd1*il1+uc1)/L1)*operate1
          du[2]=(-is1/C)*charge1*operate1
         # du[3]=((-rd1-R2-R3)*il1-RR*(I)+rd1*is1)/(L2+L3)#((L2+L3+α)*β)/Δ

          du[3]=operate1*1e6*((-I*(rpr + rr + 4.53e-7v) + il1*(-0.00388 - rd1) + is1*rd1)*((5.1042 + 0.453x)^2 - ((0.0042 + 0.453x)^2))+ ((-I*(rpr + rr + 4.53e-7v) + il2*(-0.00388 - rd2) + is2*rd2)*((0.0042 + 0.453x)^2 - (0.0042 + 0.453x)*(5.1042 + 0.453x))) +((-I*(rpr + rr + 4.53e-7v) + il3*(-0.00388 - rd3) + is3*rd3)*((0.0042 + 0.453x)^2 - (0.0042 + 0.453x)*(5.1042 + 0.453x))) )/(132.97872599999997 + 35.34758999999999x)

          du[4] =((-(R1+rs2+rd2)*is2+rd2*il2+uc2)/L1)*operate2
          du[5]=(-is2/C)*charge2*operate2
          #du[6]=start2*(((-rd2-R2-R3)*il2-RR*(I)+rd2*is2)/(L2+L3))#((L2+L3+α)*β)/Δ
          du[6]=operate2*1e6*(((-I*(rpr + rr + 4.53e-7v) + il2*(-0.00388 - rd2) + is2*rd2)*((5.1042 + 0.453x)^2 - ((0.0042 + 0.453x)^2)))+((-I*(rpr + rr + 4.53e-7v) + il1*(-0.00388 - rd1) + is1*rd1)*((0.0042 + 0.453x)^2 - (0.0042 + 0.453x)*(5.1042 + 0.453x)))+((-I*(rpr + rr + 4.53e-7v) + il3*(-0.00388 - rd3) + is3*rd3)*((0.0042 + 0.453x)^2 - (0.0042 + 0.453x)*(5.1042 + 0.453x))) )/(132.97872599999997 + 35.34758999999999x)

          du[7] =((-(R1+rs3+rd3)*is3+rd3*il3+uc3)/L1)*operate3#*charge3
          du[8]=((-is3/C)*charge3)*operate3
         # du[9]=start3*(((-rd3-R2-R3)*il3-RR*(I)+rd3*is3)/(L2+L3))#((L2+L3+α)*β)/Δ
          du[9]=operate3*1e6*(((-I*(rpr + rr + 4.53e-7v) + il3*(-0.00388 - rd3) + is3*rd3)*((5.1042 + 0.453x)^2 - ((0.0042 + 0.453x)^2)))+((-I*(rpr + rr + 4.53e-7v) + il1*(-0.00388 - rd1) + is1*rd1)*((0.0042 + 0.453x)^2 - (0.0042 + 0.453x)*(5.1042 + 0.453x)))+(-I*(rpr + rr + 4.53e-7v) + il2*(-0.00388 - rd2) + is2*rd2)*((0.0042 + 0.453x)^2 - (0.0042 + 0.453x)*(5.1042 + 0.453x))) /(132.97872599999997 + 35.34758999999999x)
         
       
          du[10]=v
          du[11]=F/m

        #=   if t-0.0009>0.0 
            start3=1.0
          end

          if t-0.0016>0.0 
            start2=1.0
           # RR=4.0e-7
          end =#
     
         #=  if t-0.001>0.0 
            operate3=0.0
            @show is3,il3
          end =#

       #=    if t-0.0006>0.0 
            operate3=1.0
          end =#

 #divT(mulTT(d[16], 1.0e6, addT, cache[2], cache[3]), addt, cache[1])
          if -(uc1)>0.0 
            charge1=0.0 # rs off not needed since charge=0
            rs1=ROff;
            rd1=ROn;
            uc1=0.0
          #=   uc1=0.0;
            is1=0.0 =#
          end 
          if -(uc2)>0.0 
            charge2=0.0
            rs2=ROff;
            rd2=ROn;
            #@show rd2
            uc2=0.0
            
          end 
         
          if -(uc3)>0.0 
            charge3=0.0
            rs3=ROff;
            rd3=ROn;
            uc3=0.0
           
          end 

         


          if -(il1)>0.0 
            operate1=0.0
            il1=0.0
          end 
          if -(il2)>0.0 
            operate2=0.0
            il2=0.0
          end 
          if -(il3)>0.0 
            operate3=0.0
            il3=0.0
          end 

    
          if t-nextT>0
            manualIntg=manualIntg+v/sqrt(t[0])
            nextT=nextT+1e-3
          #=   @show is1,is2,is3
            @show uc1,uc2,uc3
            @show il1,il2,il3 =#
          end
        
      

          
    end)
   # @show odeprob
    tspan = (0.0, 4.0e-3)
    sol= solve(odeprob,nmliqss2(),tspan,abstol=1e-3,reltol=1e-2)    
   save_Sol(sol,1)
   save_Sol(sol,2)
   save_Sol(sol,3) 
   save_Sol(sol,4)
  save_Sol(sol,5)
  save_Sol(sol,6) 
   save_Sol(sol,7) 
   save_Sol(sol,8) 
 
   save_Sol(sol,9) 
 # save_Sol(sol,9,note="z1",xlims=(0.0,0.2e-8),ylims=(-0.005,0.005)) 
   save_SolSum(sol,3,6,9) #add interpl preference
   save_Sol(sol,10) 
  save_Sol(sol,11) 
   
end
#@btime 
test()
