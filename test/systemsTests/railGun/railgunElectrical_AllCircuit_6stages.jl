using QuantizedSystemSolver
using BenchmarkTools
function test()
 
    odeprob = @NLodeProblem begin
      name=(RGElectrical6,)
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
rs11=1e-5;rs21=1e-5;rs31=1e-5;rs41=1e-5;;rs51=1e-5;rs61=1e-5;
t0=1.0

       discrete = [1e5,1.0,1.0,1e5,1.0,1.0,1e5,1.0,1.0,1e5,0.0,1.0,1e5,0.0,1.0,1e5,0.0,1.0,0.0,1e-3,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5];u = [0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,0.0]
     
       rd1=discrete[1]; operate1=discrete[2];charge1=discrete[3];
       rd2=discrete[4];operate2=discrete[5];charge2=discrete[6]; 
       rd3=discrete[7];operate3=discrete[8]; charge3=discrete[9]; 
       rd4=discrete[10];operate4=discrete[11]; charge4=discrete[12];
       rd5=discrete[13];operate5=discrete[14]; charge5=discrete[15];
       rd6=discrete[16];operate6=discrete[17]; charge6=discrete[18];
       manualIntg=discrete[19]; nextT=discrete[20];rs1=discrete[21];rs2=discrete[22];rs3=discrete[23];rs4=discrete[24];rs5=discrete[25];rs6=discrete[26];
       is1=u[1] ;uc1=u[2]; il1=u[3] ;is2=u[4] ;uc2=u[5]; il2=u[6] ;is3=u[7] ;uc3=u[8]; il3=u[9] ;is4=u[10] ;uc4=u[11]; il4=u[12];is5=u[13] ;uc5=u[14]; il5=u[15]  ;is6=u[16] ;uc6=u[17]; il6=u[18]  ;x=u[19]; v=u[20] 
       #α=LR+Lpr;
      # RR=4.0e-5
     rr=manualIntg*2.0*sqrt(μ/(3.14*γ))/w
      rpr=Rpr0*(sqrt(t0/(t+1e-4))+(t/t0)^16)/(1.0+(t/t0)^16)

    #=   rrpp=(rpr + rr + 4.53e-7v) 
      x1=(132.97872599999997 + 35.34758999999999x)
      x2=(-0.10924199999999998 - 11.782529999999998x)
      x3=678.7486367999996 + 240.3636119999999x - 8.881784197001252e-16(x^2)
      rd11=(-0.00388 - rd1);rd22=(-0.00388 - rd2);rd33=(-0.00388 - rd3);rd44=(-0.00388 - rd4) =#
       Il=il1+il2+il3+il4+il5+il6
       uf=0.1+0.2*exp(-v/100)
       F=0.5*0.453e-6*Il*Il*(1-uf*0.124)-uf*FNmec
          du[1] =((-(R1+rs1+rd1)*is1+rd1*il1+uc1)/L1)#*operate1
          du[2]=(-is1/C)*charge1#*operate1
          du[3]=operate1*1e6*  (  -(-Il*(rpr + rr + 4.53e-7v) + il3*(-0.00388 - rd3) + is3*rd3)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il4*(-0.00388 - rd4) + is4*rd4)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il5*(-0.00388 - rd5) + is5*rd5)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il2*(-0.00388 - rd2) + is2*rd2)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il6*(-0.00388 - rd6) + is6*rd6)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) + (-Il*(rpr + rr + 4.53e-7v) + il1*(-0.00388 - rd1) + is1*rd1)*(0.19607843137254904 + (-0.0042 - 0.453x) / (5.1(5.1251999999999995 + 2.718x))))


          du[4] =((-(R1+rs2+rd2)*is2+rd2*il2+uc2)/L1)*operate2
          du[5]=(-is2/C)*charge2*operate2
          du[6]=operate2*1e6* ( -(-Il*(rpr + rr + 4.53e-7v) + il3*(-0.00388 - rd3) + is3*rd3)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il4*(-0.00388 - rd4) + is4*rd4)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il5*(-0.00388 - rd5) + is5*rd5)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il1*(-0.00388 - rd1) + is1*rd1)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il6*(-0.00388 - rd6) + is6*rd6)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) + (-Il*(rpr + rr + 4.53e-7v) + il2*(-0.00388 - rd2) + is2*rd2)*(0.19607843137254904 + (-0.0042 - 0.453x) / (5.1(5.1251999999999995 + 2.718x))) )

          du[7] =((-(R1+rs3+rd3)*is3+rd3*il3+uc3)/L1)*operate3#*charge3
          du[8]=((-is3/C)*charge3)*operate3
          du[9]=operate3*1e6* (-(-Il*(rpr + rr + 4.53e-7v) + il4*(-0.00388 - rd4) + is4*rd4)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il5*(-0.00388 - rd5) + is5*rd5)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il1*(-0.00388 - rd1) + is1*rd1)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il2*(-0.00388 - rd2) + is2*rd2)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il6*(-0.00388 - rd6) + is6*rd6)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) + (-Il*(rpr + rr + 4.53e-7v) + il3*(-0.00388 - rd3) + is3*rd3)*(0.19607843137254904 + (-0.0042 - 0.453x) / (5.1(5.1251999999999995 + 2.718x))) )
         
          du[10] =((-(R1+rs4+rd4)*is4+rd4*il4+uc4)/L1)*operate4#*charge3
          du[11]=((-is4/C)*charge4)*operate4
          du[12]=operate4*1e6* (-(-Il*(rpr + rr + 4.53e-7v) + il3*(-0.00388 - rd3) + is3*rd3)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il5*(-0.00388 - rd5) + is5*rd5)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il1*(-0.00388 - rd1) + is1*rd1)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il2*(-0.00388 - rd2) + is2*rd2)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il6*(-0.00388 - rd6) + is6*rd6)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) + (-Il*(rpr + rr + 4.53e-7v) + il4*(-0.00388 - rd4) + is4*rd4)*(0.19607843137254904 + (-0.0042 - 0.453x) / (5.1(5.1251999999999995 + 2.718x))))
  
          du[13] =((-(R1+rs5+rd5)*is5+rd5*il5+uc5)/L1)*operate5#*charge3
          du[14]=((-is5/C)*charge5)*operate5
          du[15]=operate5*1e6* (-(-Il*(rpr + rr + 4.53e-7v) + il3*(-0.00388 - rd3) + is3*rd3)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il4*(-0.00388 - rd4) + is4*rd4)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il1*(-0.00388 - rd1) + is1*rd1)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il2*(-0.00388 - rd2) + is2*rd2)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il6*(-0.00388 - rd6) + is6*rd6)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) + (-Il*(rpr + rr + 4.53e-7v) + il5*(-0.00388 - rd5) + is5*rd5)*(0.19607843137254904 + (-0.0042 - 0.453x) / (5.1(5.1251999999999995 + 2.718x))))
  
          du[16] =((-(R1+rs6+rd6)*is6+rd6*il6+uc6)/L1)*operate6#*charge3
          du[17]=((-is6/C)*charge6)*operate6
          du[18]=operate6*1e6* (-(-Il*(rpr + rr + 4.53e-7v) + il3*(-0.00388 - rd3) + is3*rd3)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il4*(-0.00388 - rd4) + is4*rd4)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il5*(-0.00388 - rd5) + is5*rd5)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il1*(-0.00388 - rd1) + is1*rd1)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) - (-Il*(rpr + rr + 4.53e-7v) + il2*(-0.00388 - rd2) + is2*rd2)*((0.0042 + 0.453x) / (5.1(5.1251999999999995 + 2.718x))) + (-Il*(rpr + rr + 4.53e-7v) + il6*(-0.00388 - rd6) + is6*rd6)*(0.19607843137254904 + (-0.0042 - 0.453x) / (5.1(5.1251999999999995 + 2.718x))))
  

          du[19]=v #    =u[20]
          du[20]=F/m

          if t-0.0003>0.0 
            operate4=1.0
          end 

          if t-0.0007>0.0 
            operate5=1.0
          end 
          if t-0.001>0.0 
            operate6=1.0
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




          if t-nextT>0
            manualIntg=manualIntg+v/sqrt(t[0])
            nextT=nextT+1e-3
       
          end
    

          
    end
  #  @show odeprob
  println("start solving")
    tspan = (0.0, 5.0e-3)
    sol= solve(odeprob,nmliqss1(),tspan,abstol=1e-3,reltol=1e-2) 
  #  sol= solve(odeprob,nmliqss1(),tspan,abstol=1e-3,reltol=1e-2)    
   save_Sol(sol,1)
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
  save_Sol(sol,15) 
  save_Sol(sol,16) 
  save_Sol(sol,17) 
  save_Sol(sol,18) 
 # save_Sol(sol,9,note="z1",xlims=(0.0,0.2e-8),ylims=(-0.005,0.005)) 
   save_SolSum(sol,3,6,9,12,15,18,interp=0.00001) #add interpl preference
   
  save_Sol(sol,19) 
  save_Sol(sol,20) 
   
end
#@btime 
test()

function convert(par)
  exp_as_str = " a = $par"

  exp = Meta.parse(exp_as_str)
   
  return exp
end

exp = convert(156)

function test(exp)
  println(exp)
end

test(convert(6))