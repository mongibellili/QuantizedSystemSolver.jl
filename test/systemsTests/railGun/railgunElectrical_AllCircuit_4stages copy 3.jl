using QuantizedSystemSolver
#using BenchmarkTools
function test()
 
  odeprob = NLodeProblem(
    quote
      name=(RGElectrical4,)
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
Rpr0 = 1.5e-6
FNmec = 680.0; #α = 0.154

t0=1.0

discrete = [1e5,1e-5,1.0,1.0,1e5,1e-5,1.0,1.0,1e5,1e-5,1.0,1.0,1e5,1e-5,1.0,1.0,0.0,1e-3]
rd1=discrete[1]; rs1=discrete[2]; operate1=discrete[3];  charge1=discrete[4]
rd2=discrete[5]; rs2=discrete[6]; operate2=discrete[7];  charge2=discrete[8]; 
rd3=discrete[9]; rs3=discrete[10];operate3=discrete[11]; charge3=discrete[12]; 
rd4=discrete[13];rs4=discrete[14];operate4=discrete[15]; charge4=discrete[16];
       manualIntg=discrete[17]; nextT=discrete[18]
       u = [0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,0.0]
       is1=u[1] ;uc1=u[2]; il1=u[3] ;is2=u[4] ;uc2=u[5]; il2=u[6] ;is3=u[7] ;uc3=u[8]; il3=u[9] ;is4=u[10] ;uc4=u[11]; il4=u[12] ;x=u[13]; v=u[14] 
       #α=LR+Lpr;
      # RR=4.0e-5
     rr=manualIntg*1.1925695879998878e-8
      rpr=Rpr0*(sqrt(t0/(t+1e-4))+(t)^16)/(1.0+(t)^16)

      rrpp=(rpr + rr + 4.53e-7v) 
    #=   x1=(132.97872599999997 + 35.34758999999999x)
      x2=(-0.10924199999999998 - 11.782529999999998x)
      x3=678.7486367999996 + 240.3636119999999x - 8.881784197001252e-16(x^2) =#
      x4=(5.1(5.1168 + 1.812x))
      x5=(0.0042 + 0.453x)
      rd11=(-0.00388 - rd1);rd22=(-0.00388 - rd2);rd33=(-0.00388 - rd3);rd44=(-0.00388 - rd4)
       Il=il1+il2+il3+il4
       uf=0.1+0.2*exp(-v/100.0)
       F=2.265e-7*Il*Il*(1.0-uf*0.124)-uf*FNmec
          du[1] =((-(R1+rs1+rd1)*is1+rd1*il1+uc1)/L1)*operate1#*charge1
          du[2]=(-is1/C)*charge1*operate1
          du[3]=operate1*1e6*  (-(-4*Il*(rrpp)*(1.0-0.04901960784313726 *x4/x5)+il2*rd22+is2*rd2+il3*rd33+is3*rd3+il1*rd11+is1*rd1+il4*rd44+is4*rd4)*(x5 / x4)+ (il1*rd11+is1*rd1)*(0.19607843137254904 ))

          du[4] =((-(R1+rs2+rd2)*is2+rd2*il2+uc2)/L1)*operate2#*charge2
          du[5]=(-is2/C)*charge2*operate2
          du[6]=operate2*1e6* ( -(-4*Il*(rrpp)*(1.0-0.04901960784313726 *x4/x5) + il2*rd22+is2*rd2+il3*rd33+is3*rd3+il1*rd11+is1*rd1+il4*rd44+is4*rd4)*(x5 / x4)+( il2*rd22 + is2*rd2)*(0.19607843137254904 ))

          du[7] =((-(R1+rs3+rd3)*is3+rd3*il3+uc3)/L1)*operate3#*charge3
          du[8]=((-is3/C)*charge3)*operate3
          du[9]=operate3*1e6* ( -(-4*Il*(rrpp)*(1.0-0.04901960784313726 *x4/x5) + il2*rd22+is2*rd2+il3*rd33+is3*rd3+il1*rd11+is1*rd1+il4*rd44+is4*rd4)*(x5 / x4)+( il3*rd33 + is3*rd3)*(0.19607843137254904 ))
         
          du[10] =((-(R1+rs4+rd4)*is4+rd4*il4+uc4)/L1)*operate4#*charge3
          du[11]=((-is4/C)*charge4)*operate4
          du[12]=operate4*1e6* ( -(-4*Il*(rrpp)*(1.0-0.04901960784313726 *x4/x5) + il2*rd22+is2*rd2+il3*rd33+is3*rd3+il1*rd11+is1*rd1+il4*rd44+is4*rd4)*(x5 / x4)+ ( il4*rd44 + is4*rd4)*(0.19607843137254904 ))
  
          du[13]=v
          du[14]=F/m

         #=  if t-0.00025>0.0 
            operate3=1.0
          end 

          if t-0.0003>0.0 
            operate4=1.0
          end  =#
      

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
    
   

          if t-nextT>0
            manualIntg=manualIntg+v/sqrt(t[0])
            nextT=nextT+1e-3
       
          end
    

          
    end)
   # @show odeprob.eqs
    tspan = (0.0, 4.0e-3)
   #println("start solving")
    sol= solve(odeprob,nmliqss2(),tspan,abstol=1e-3,reltol=1e-2)    
 #  save_Sol(sol,1)
 #  save_Sol(sol,2)
#=   save_Sol(sol,3)
 msg="_"
  #     save_Sol(sol,4,note="$msg")
 # save_Sol(sol,5)
  save_Sol(sol,6,note="$msg") 
 #  save_Sol(sol,7,note="$msg") 
  # save_Sol(sol,8) 
 
   save_Sol(sol,9,note="$msg") 
 #  save_Sol(sol,10,note="$msg") 
 # save_Sol(sol,11) 
  save_Sol(sol,12,note="$msg") 
 # save_Sol(sol,12,note="q_",xlims=(1.326e-6,1.33e-6),ylims=(221.1,221.75)) 
   save_SolSum(sol,3,6,9,12,interp=0.00001) #add interpl preference 
   
  save_Sol(sol,13) 
  save_Sol(sol,14) 
    =#
end
@time test()
