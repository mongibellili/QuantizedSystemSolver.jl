using QuantizedSystemSolver
using BenchmarkTools
function test()
 
  odeprob = NLodeProblem(
    quote
      name=(RGElectrical12,)
      ROn = 1e-5;ROff = 1e1;
      # Lpr = 4.2e-9#0.48*1e-6#0.45 * 1e-6
L1 = 0.6e-4 #28.0*1e-6#0.6*1e-6
#L2 = 4.0e-6;L3 = 1.1*1e-6
#L23=5.1e-6

R1= 4.0e-3; #R2 = 0.28e-3;R3 = 3.6e-3

C = 3.08e-3#3.08*1e-3

m=0.12
#= γ = 50.0e6; w = 7.5#25.0*1e-3#15.0*1e-3 
μ = 4.0*3.14*1e-7 =#
#rr=manualIntg*sqrt(μ/(3.14*γ))/w=manualIntg*coef =manualIntg* 1.1925695879998878e-8
Rpr0 = 2.5e-6
FNmec = 680.0; #α = 0.15
      t0=1.0
       discrete = [1e5,1e-5,1.0,1.0,1e5,1e-5,1.0,1.0,1e5,1e-5,1.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,0.0,1e-3]
       rd1=discrete[1]; rs1=discrete[2]; operate1=discrete[3];  charge1=discrete[4]
       rd2=discrete[5]; rs2=discrete[6]; operate2=discrete[7];  charge2=discrete[8]; 
       rd3=discrete[9]; rs3=discrete[10];operate3=discrete[11]; charge3=discrete[12]; 
       rd4=discrete[13];rs4=discrete[14];operate4=discrete[15]; charge4=discrete[16];
       rd5=discrete[17];rs5=discrete[18];operate5=discrete[19]; charge5=discrete[20];
       rd6=discrete[21];rs6=discrete[22];operate6=discrete[23]; charge6=discrete[24];
       rd7=discrete[25];rs7=discrete[26];operate7=discrete[27]; charge7=discrete[28];
       rd8=discrete[29];rs8=discrete[30];operate8=discrete[31]; charge8=discrete[32];
       rd9=discrete[33];rs9=discrete[34];operate9=discrete[35]; charge9=discrete[36];
       rd10=discrete[37];rs10=discrete[38];operate10=discrete[39]; charge10=discrete[40];
       rd11=discrete[41];rs11=discrete[42];operate11=discrete[43]; charge11=discrete[44];
       rd12=discrete[45];rs12=discrete[46];operate12=discrete[47]; charge12=discrete[48];
       manualIntg=discrete[49]; nextT=discrete[50];
       u = [0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,0.0]    
       is1=u[1] ;uc1=u[2]; il1=u[3] ;is2=u[4] ;uc2=u[5]; il2=u[6] ;is3=u[7] ;uc3=u[8]; il3=u[9] ;is4=u[10] ;uc4=u[11]; il4=u[12];is5=u[13] ;uc5=u[14]; il5=u[15]  ;is6=u[16] ;uc6=u[17]; il6=u[18];is7=u[19] ;uc7=u[20]; il7=u[21] ;is8=u[22] ;uc8=u[23]; il8=u[24] ;is9=u[25] ;uc9=u[26]; il9=u[27]  ;is10=u[28] ;uc10=u[29]; il10=u[30];is11=u[31] ;uc11=u[32]; il11=u[33];is12=u[34] ;uc12=u[35]; il12=u[36] ;x=u[37]; v=u[38] 
       #α=LR+Lpr;
      # RR=4.0e-5
      rr=manualIntg*1.1925695879998878e-8
      rpr=Rpr0*(sqrt(t0/(t+1e-4))+(t)^16)/(1.0+(t)^16)

      rrpp=(rpr + rr + 4.53e-7v)

    #=   rrpp=rrpp 
      x1=(132.97872599999997 + 35.34758999999999x)
      x2=(-0.10924199999999998 - 11.782529999999998x)
      x3=678.7486367999996 + 240.3636119999999x - 8.881784197001252e-16(x^2)
      rd_11=rd_11;rd22=rd22;rd33=rd33;rd44=rd44 =#
         Il=il1+il2-0.0+il3+il4+il5+il6-0.0+il7+il8+il9+il10-0.0+il11+il12
         uf=0.1+0.2*exp(-v/100)
         F=0.5*0.453e-6*Il*Il*(1-uf*0.124)-uf*FNmec
         x4=(5.1(5.150399999999999 + 5.436000000000001x))
         x5=(0.0042 + 0.453x)
         
         rd_11=(-0.00388 - rd1);rd22=(-0.00388 - rd2);rd33=(0.00388 + rd3);rd44=(-0.00388 - rd4); rd55=(-0.00388 - rd5);rd66=(0.00388 + rd6);rd77=(-0.00388 - rd7);rd88=(-0.00388 - rd8); rd99=(0.00388 + rd9);rd100=(-0.00388 - rd10);rd110=(-0.00388 - rd11);rd120=(-0.00388 - rd12)
         isid=il1*rd_11 +is1*rd1+il2*rd22+is2*rd2-il3*rd33+is3*rd3+il4*rd44+is4*rd4+il5*rd55 + is5*rd5-il6*rd66+is6*rd6+il7*rd77+is7*rd7+il8*rd88+is8*rd8-il9*rd99+is9*rd9+il10*rd100+is10*rd10+il11*rd110+is11*rd11+il12*rd120+is12*rd12
         
         du[1] =((-(R1+rs1+rd1)*is1+rd1*il1+uc1)/L1)*operate1
          du[2]=(-is1/C)*charge1*operate1
          du[3]=operate1*1e6*  (-(-12.0*Il*rrpp*(1.0-0.016339869281045753*x4/x5) + (isid))*( x5 / x4) +(il1*rd_11+is1*rd1)*(0.19607843137254904 ))


          du[4] =((-(R1+rs2+rd2)*is2+rd2*il2+uc2)/L1)*operate2
          du[5]=(-is2/C)*charge2*operate2
          du[6]=operate2*1e6* ( -(-12.0*Il*rrpp*(1.0-0.016339869281045753*x4/x5) + (isid))*( x5 / x4) +(il2*rd22+is2*rd2)*(0.19607843137254904 ))

          du[7] =((-(R1+rs3+rd3)*is3+rd3*il3+uc3)/L1)*operate3#*charge3
          du[8]=((-is3/C)*charge3)*operate3
          du[9]=operate3*1e6* (-(-12.0*Il*rrpp*(1.0-0.016339869281045753*x4/x5) + (isid))*( x5 / x4) +(-il3*rd33+is3*rd3)*(0.19607843137254904 ))
         
          du[10] =((-(R1+rs4+rd4)*is4+rd4*il4+uc4)/L1)*operate4#*charge3
          du[11]=((-is4/C)*charge4)*operate4
          du[12]=operate4*1e6* ( -(-12.0*Il*rrpp*(1.0-0.016339869281045753*x4/x5) + (isid))*( x5 / x4) +(il4*rd44+is4*rd4)*(0.19607843137254904 ))
  
          du[13] =((-(R1+rs5+rd5)*is5+rd5*il5+uc5)/L1)*operate5#*charge3
          du[14]=((-is5/C)*charge5)*operate5
          du[15]=operate5*1e6* ( -(-12.0*Il*rrpp*(1.0-0.016339869281045753*x4/x5) + (isid))*( x5 / x4) +(il5*rd55+is5*rd5)*(0.19607843137254904 ))
  
          du[16] =((-(R1+rs6+rd6)*is6+rd6*il6+uc6)/L1)*operate6#*charge3
          du[17]=((-is6/C)*charge6)*operate6
          du[18]=operate6*1e6* ( -(-12.0*Il*rrpp*(1.0-0.016339869281045753*x4/x5) + (isid))*( x5 / x4) +(-il6*rd66+is6*rd6)*(0.19607843137254904 ))
  
          du[19] =((-(R1+rs7+rd7)*is7+rd7*il7+uc7)/L1)*operate7
          du[20]=(-is7/C)*charge7*operate7
          du[21]=operate7*1e6*  ( -(-12.0*Il*rrpp*(1.0-0.016339869281045753*x4/x5) + (isid))*( x5 / x4) +(il7*rd77+is7*rd7)*(0.19607843137254904 ))

          du[22] =((-(R1+rs8+rd8)*is8+rd8*il8+uc8)/L1)*operate8
          du[23]=(-is8/C)*charge8*operate8
          du[24]=operate8*1e6*  ( -(-12.0*Il*rrpp*(1.0-0.016339869281045753*x4/x5) + (isid))*( x5 / x4) +(il8*rd88+is8*rd8)*(0.19607843137254904 ))

          du[25] =((-(R1+rs9+rd9)*is9+rd9*il9+uc9)/L1)*operate9
          du[26]=(-is9/C)*charge9*operate9
          du[27]=operate9*1e6*  ( -(-12.0*Il*rrpp*(1.0-0.016339869281045753*x4/x5) + (isid))*( x5 / x4) +(-il9*rd99+is9*rd9)*(0.19607843137254904 ))

          du[28] =((-(R1+rs10+rd10)*is10+rd10*il10+uc10)/L1)*operate10
          du[29]=(-is10/C)*charge10*operate10
          du[30]=operate10*1e6*  ( -(-12.0*Il*rrpp*(1.0-0.016339869281045753*x4/x5) + (isid))*( x5 / x4) +(il10*rd100+is10*rd10)*(0.19607843137254904 ))

          du[31] =((-(R1+rs11+rd11)*is11+rd11*il11+uc11)/L1)*operate11
          du[32]=(-is11/C)*charge11*operate11
          du[33]=operate11*1e6*  ( -(-12.0*Il*rrpp*(1.0-0.016339869281045753*x4/x5) + (isid))*( x5 / x4) +(il11*rd110+is11*rd11)*(0.19607843137254904 ))

          du[34] =((-(R1+rs12+rd12)*is12+rd12*il12+uc12)/L1)*operate12
          du[35]=(-is12/C)*charge12*operate12
          du[36]=operate12*1e6*  ( -(-12.0*Il*rrpp*(1.0-0.016339869281045753*x4/x5) + (isid))*( x5 / x4) +(il12*rd120+is12*rd12)*(0.19607843137254904 ))


          du[37]=v #    v=u[26]
          du[38]=F/m

          if t-0.00019>0.0 
            operate4=1.0
          end 

          if t-0.00038>0.0 
            operate5=1.0
          end 
          if t-0.00057>0.0 
            operate6=1.0
          end 
          if t-0.00095>0.0 
            operate7=1.0
          end 
          if t-0.00115>0.0 
            operate8=1.0
          end 
          if t-0.00142>0.0 
            operate9=1.0
          end 
          if t-0.00175>0.0 
            operate10=1.0
          end 
          if t-0.00195>0.0 
            operate11=1.0
          end 
          if t-0.00210>0.0 
            operate12=1.0
          end 

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
          if -(uc4)>0.0 
            charge4=0.0
            rs4=ROff;
            rd4=ROn;
            uc4=0.0
           
          end 
         
          if -(uc5)>0.0 
            charge5=0.0
            rs5=ROff;
            rd5=ROn;
            uc5=0.0
           
          end 
          if -(uc6)>0.0 
            charge6=0.0
            rs6=ROff;
            rd6=ROn;
            uc6=0.0
           
          end 
          if -(uc7)>0.0 
            charge7=0.0
            rs7=ROff;
            rd7=ROn;
            uc7=0.0
           
          end 
          if -(uc8)>0.0 
            charge8=0.0
            rs8=ROff;
            rd8=ROn;
            uc8=0.0
           
          end 
          if -(uc9)>0.0 
            charge9=0.0
            rs9=ROff;
            rd9=ROn;
            uc9=0.0
          end 
          if -(uc10)>0.0 
            charge10=0.0
            rs10=ROff;
            rd10=ROn;
            uc10=0.0
          end 
          if -(uc11)>0.0 
            charge11=0.0
            rs11=ROff;
            rd11=ROn;
            uc11=0.0
          end 
          if -(uc12)>0.0 
            charge12=0.0
            rs12=ROff;
            rd12=ROn;
            uc12=0.0
          end 
    #=       if -(il1)>0.0 
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
          if -(il4)>0.0 
            operate4=0.0
            il4=0.0
          end 
    
          if -(il5)>0.0 
            operate5=0.0
            il5=0.0
          end 
          if -(il6)>0.0 
            operate6=0.0
            il6=0.0
          end 

          if -(il7)>0.0 
            operate7=0.0
            il7=0.0
          end 
          if -(il8)>0.0 
            operate8=0.0
            il8=0.0
          end 
          if -(il9)>0.0 
            operate9=0.0
            il9=0.0
          end 
          if -(il10)>0.0 
            operate10=0.0
            il10=0.0
          end 
          if -(il11)>0.0 
            operate11=0.0
            il11=0.0
          end 
          if -(il12)>0.0 
            operate12=0.0
            il12=0.0
          end  =#
          if t-nextT>0
            manualIntg=manualIntg+v/sqrt(t[0])
            nextT=nextT+1e-3
       
          end
    
        end
          
    )
    #@show odeprob.eqs
  println("start solving")
    tspan = (0.0, 9.0e-3)
   sol= solve(odeprob,nmliqss2(),tspan,abstol=1e-3,reltol=1e-2) 
  #  sol= solve(odeprob,nmliqss1(),tspan,abstol=1e-3,reltol=1e-2)    
#=    save_Sol(sol,1)
   save_Sol(sol,2)
   save_Sol(sol,3) 

   save_Sol(sol,4)
  save_Sol(sol,5)
  save_Sol(sol,6) 

   save_Sol(sol,7) 
   save_Sol(sol,8) 
   save_Sol(sol,9) 

   save_Sol(sol,10) 
  save_Sol(sol,11) 
  save_Sol(sol,12) 

  save_Sol(sol,13) 
  save_Sol(sol,14) 
  save_Sol(sol,15)  =#
 # save_Sol(sol,16) 
 # save_Sol(sol,17) 
#=   save_Sol(sol,18) 
 # save_Sol(sol,19) 
 # save_Sol(sol,20) 
  save_Sol(sol,21) 
 # save_Sol(sol,22) 
 # save_Sol(sol,23) 
  save_Sol(sol,24) 
 # save_Sol(sol,25) 
 # save_Sol(sol,26) 
  save_Sol(sol,27) 
 # save_Sol(sol,28) 
 # save_Sol(sol,29) 
  save_Sol(sol,30)  =#

 # save_Sol(sol,9,note="z1",xlims=(0.0,0.2e-8),ylims=(-0.005,0.005)) 
   save_SolSum(sol,3,6,9,12,15,18,21,24,27,30,33,36,interp=0.00001) #add interpl preference
   

 #= save_Sol(sol,31) 
  save_Sol(sol,32) 
  save_Sol(sol,33) 
 # save_Sol(sol,34) 
 # save_Sol(sol,35) 
  save_Sol(sol,36) 
   =#
  save_Sol(sol,37) 
  save_Sol(sol,38) 
   
end
#@time 
#test()
#@btime 
test()
