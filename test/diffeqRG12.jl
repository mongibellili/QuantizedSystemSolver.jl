using DifferentialEquations
using Plots
function odeDiffEquPackage()
    function f(du, u, p, t)
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
         p = [1e5,1e-5,1.0,1.0,1e5,1e-5,1.0,1.0,1e5,1e-5,1.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,0.0,1e-3]
         rd1=p[1]; rs1=p[2]; operate1=p[3];  charge1=p[4]
         rd2=p[5]; rs2=p[6]; operate2=p[7];  charge2=p[8]; 
         rd3=p[9]; rs3=p[10];operate3=p[11]; charge3=p[12]; 
         rd4=p[13];rs4=p[14];operate4=p[15]; charge4=p[16];
         rd5=p[17];rs5=p[18];operate5=p[19]; charge5=p[20];
         rd6=p[21];rs6=p[22];operate6=p[23]; charge6=p[24];
         rd7=p[25];rs7=p[26];operate7=p[27]; charge7=p[28];
         rd8=p[29];rs8=p[30];operate8=p[31]; charge8=p[32];
         rd9=p[33];rs9=p[34];operate9=p[35]; charge9=p[36];
         rd10=p[37];rs10=p[38];operate10=p[39]; charge10=p[40];
         rd11=p[41];rs11=p[42];operate11=p[43]; charge11=p[44];
         rd12=p[45];rs12=p[46];operate12=p[47]; charge12=p[48];
         manualIntg=p[49]; nextT=p[50];
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
    end
    function condition(out, u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        out[1] = t-0.00019
        out[2] =  t-0.00038
        out[3] = t-0.00057
        out[4] = t-0.00095
        out[5] = t-0.00115
        out[6] = t-0.00142
        out[7] = t-0.00175
        out[8] =  t-0.00195
        out[9] = t-0.00210

        out[10] = -u[2]
        out[11] = -u[5]
        out[12] = -u[8]
        out[13] = -u[11]
        out[14] =  -u[14]
        out[15] = -u[17]
        out[16] = -u[20]
        out[17] = -u[23]
        out[18] = -u[26]
        out[19] = -u[29]
        out[20] =  -u[32]
        out[21] = -u[35]
        out[22] = t-p[50]
       
    end
  
    function affect!(integrator, idx)
        if idx == 1
            p[15]=1.0
        elseif idx == 2
            p[19]=1.0
        elseif idx == 3
            p[23]=1.0
        elseif idx == 4
            p[27]=1.0
        elseif idx == 5
            p[31]=1.0
        elseif idx == 6
            p[35]=1.0
        elseif idx == 7
            p[39]=1.0
        elseif idx == 8
            p[43]=1.0
        elseif idx == 9
            p[47]=1.0
        elseif idx == 10
            p[4]=0.0
            p[2]=10.0
            p[1]=1e-5
            integrator.u[2]=0.0
        elseif idx == 11
            p[8]=0.0
            p[6]=10.0
            p[5]=1e-5
            integrator.u[5]=0.0
        elseif idx == 12
            p[12]=0.0
            p[10]=10.0
            p[9]=1e-5
            integrator.u[8]=0.0
        elseif idx == 13
            p[16]=0.0
            p[14]=10.0
            p[13]=1e-5
            integrator.u[11]=0.0
        elseif idx == 14
            p[20]=0.0
            p[18]=10.0
            p[17]=1e-5
            integrator.u[14]=0.0
        elseif idx == 15
            p[24]=0.0
            p[22]=10.0
            p[21]=1e-5
            integrator.u[17]=0.0
        elseif idx == 16
            p[28]=0.0
            p[26]=10.0
            p[25]=1e-5
            integrator.u[20]=0.0
        elseif idx == 17
            p[32]=0.0
            p[30]=10.0
            p[29]=1e-5
            integrator.u[23]=0.0
        elseif idx == 18
            p[36]=0.0
            p[34]=10.0
            p[33]=1e-5
            integrator.u[26]=0.0
        elseif idx == 19
            p[40]=0.0
            p[38]=10.0
            p[37]=1e-5
            integrator.u[29]=0.0
        elseif idx == 20
            p[44]=0.0
            p[42]=10.0
            p[41]=1e-5
            integrator.u[32]=0.0
        elseif idx == 21
            p[48]=0.0
            p[46]=10.0
            p[45]=1e-5
            integrator.u[35]=0.0
        elseif idx == 22
            p[49]=p[49]+integrator.u[38]/sqrt(integrator.t)
            p[50]=p[50]+1e-3
      
        end
    end

  
    cbs = VectorContinuousCallback(condition, affect!, 22)
    u0 = [0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,10750.0,0.0,0.0,0.0]    
    tspan = (0.0, 9.0e-3)
    p =  [1e5,1e-5,1.0,1.0,1e5,1e-5,1.0,1.0,1e5,1e-5,1.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,1e5,1e-5,0.0,1.0,0.0,1e-3]
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, Rodas5P(), callback = cbs, reltol=1e-3,abstol=1e-4#= ,dtmax=1e-12 =#,maxiters=Int(1e12)#= dt = 1e-3, adaptive = false =#)
    p1=plot(sol,idxs=1);
    savefig(p1, "Rodas5P()_vector_rg12_1")
    p2=plot(sol,idxs=2);
    savefig(p2, "Rodas5P()_vector_rg12_2")
    p3=plot(sol,idxs=3);
    savefig(p3, "Rodas5P()_vector_rg12_3")
 end
 odeDiffEquPackage() 