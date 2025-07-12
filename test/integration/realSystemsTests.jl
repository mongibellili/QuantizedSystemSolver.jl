

function buck(dy,y,p,t)# api requires four args
  C = 1e-4; L = 1e-4; R = 10.0;U = 24.0; T = 1e-4; DC = 0.5; ROn = 1e-5;ROff = 1e5;RR=0.1;LL=1e4;CC=1e4
  #rd=p[1];rs=p[2];nextT=p[3];lastT=p[4];#diodeon=p[5]
  rd,rs,nextT,lastT=p
  il=y[1] ;uc=y[2]
  id=(il*rs-U)/(rd+rs) # diode's current
  dy[1] =(-id*rd-uc)*LL
  dy[2]=(il-uc*RR)*CC
  if t-nextT>0.0 
    lastT=nextT;nextT=nextT+T;rs=ROn
  end
  if t-lastT-DC*T>0.0 
    rs=ROff
  end                          
  #if diodeon*(id)+(1.0-diodeon)*(id)>0
  if (id)>0
    rd=ROn;#diodeon=1.0
  else
    rd=ROff;#diodeon=0.0
  end     
end
tspan = (0.0,0.001)
p = [1e5,1e-5,1e-4,0.0];u0 = [0.0,0.0]
odeprob = ODEProblem(buck,u0,tspan,p)
sol= solve(odeprob,liqss2(),abstol=1e-4,reltol=1e-2) 
@test  sol.stats.totalSteps<500
#save_Sol(sol)
@test 19.2<sol(0.0005,idxs=2)<19.9
sol= solve(odeprob,nmliqss2(),abstol=1e-3,reltol=1e-2,detection=Detection(3))   
@test  sol.stats.totalSteps<500
@test 19.2<sol(0.0005,idxs=2)<19.9

function adr(du,u,p,t)
   _dx=100.0#1/dx=N/10=1000/10
  a=1.0;d=0.1;r=1000.0
  #p=[0.0]
  du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
  for k in 2:999  
    du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
  end 
  du[1000]=-a*_dx*(u[1000]-u[999])+d*_dx*_dx*(2.0*u[999]-2.0*u[1000])+r*u[1000]*u[1000]*(1.0-u[1000]) 
  #=  if u[1]-10.0>0.0 #fake to test discreteintgrator & loop
  p[1]=1.0
  end =#
end
tspan=(0.0,5.0)
u=zeros(1000)
u[1:333].=1.0
odeprob=ODEProblem(adr,u,tspan)
sol=solve(odeprob,nmliqss1(),abstol=1e-5,reltol=1e-2)#
@test  sol.stats.totalSteps<600000
sol=solve(odeprob,abstol=1e-5,reltol=1e-2,tspan)#
@test  sol.stats.totalSteps<160000
BSON.@load "solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01
solnmliqssInterp=solInterpolated(sol,0.01)
getErrorByRefs(solnmliqssInterp,1,solFeagin14VectorN1000d01)
err4=getAverageErrorByRefs(solnmliqssInterp,solFeagin14VectorN1000d01)
@test err4<0.05
@test 0.33<sol(1.5,idxs=1)<0.39
@test 0.62<sol(1.5,idxs=2)<0.67
@test 0.92<sol(1.5,idxs=400)<1.0
@test 0.92<sol(1.5,idxs=600)<1.0
@test 0.92<sol(1.5,idxs=1000)<1.0


function tyson(du,u,p,t)
  x=[1e6,1e3]
du[1] = u[4]-x[1]*u[1]+x[2]*u[2]
du[2] =-200.0*u[2]*u[5]+1e6*u[1]-1e3*u[2]
du[3] = 200.0*u[2]*u[5]-u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)
du[4] =u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)-u[4]
du[5] = 0.015-200.0*u[2]*u[5]
du[6] =u[4]-0.6*u[6]
end 
tspan=(0.0,25.0)
u = [0.0,0.75,0.25,0.0,0.0,0.0]
odeprob=ODEProblem(tyson,u,tspan)
sol=solve(odeprob,nmliqss2(),abstol=1e-5,reltol=1e-4,detection=Detection(3))
@test  sol.stats.totalSteps<2000
@test 0.00001<sol(20.0,idxs=1)<0.01
@test 0.8<sol(20.0,idxs=2)<0.99
@test 0.03<sol(20.0,idxs=3)<0.09
@test 0.00001<sol(20.0,idxs=4)<0.01
@test 0.0<sol(20.0,idxs=5)<8.0e-4
@test 0.01<sol(20.0,idxs=6)<0.03
solInterpolated(sol,1,0.01)
BSON.@load "solRodas5PVectorTyson.bson" solRodas5PVectorTyson
solnmliqssInterp=solInterpolated(sol,0.01)
err=getAverageErrorByRefs(solnmliqssInterp,solRodas5PVectorTyson) 
@test err<1.8

function cuk4(du,u,p,t)
   C = 1e-4; L = 1e-4; R = 10.0;U = 24.0; T = 1e-4; DC = 0.25; ROn = 1e-5;ROff = 1e5;L1=1e-4;C1=1e-4;C2 = 1e-4;L2 = 1e-4;
   #p Rd(start=1e5), Rs(start=1e-5), nextT(start=T),lastT,diodeon;
   Rd=p[1];Rs=p[2];nextT=p[3];lastT=p[4];diodeon=p[5]
   uc2=u[13]
   il2_1=u[i] ;il2_2=u[i-4] ;il2_3=u[i-8] ;il1_1=u[i+4] ;il1_2=u[i] ;il1_3=u[i-4] ;uc1_1=u[i+8];uc1_2=u[i+4] ;uc1_3=u[i] ;
   id1=(((il2_1+il1_1)*Rs-uc1_1)/(Rd+Rs))
   id2=(((il2_2+il1_2)*Rs-uc1_2)/(Rd+Rs))
   id3=(((il2_3+il1_3)*Rs-uc1_3)/(Rd+Rs))
  for i=1:4    #il2
    du[i] =(-uc2-Rs*id1)/L2
  end
  for i=5:8    #il1
    du[i]=(U-uc1_2-id2*Rs)/L1
  end
  for i=9:12    #uc1
    du[i]=(id3-il2_3)/C1
  end
  du[13]=(u[1]+u[2]+u[3]+u[4]-uc2/R)/C2
  if t-nextT>0.0 
      lastT=nextT
      nextT=nextT+T
      Rs=ROn
  end
  if t-lastT-DC*T>0.0 
      Rs=ROff 
  end                          
  if diodeon*(((u[1]+u[5])*Rs-u[9])/(Rd+Rs))+(1.0-diodeon)*(((u[1]+u[5])*Rs-u[9])*Rd/(Rd+Rs))>0
    Rd=ROn
    diodeon=1.0
  else
    Rd=ROff
    diodeon=0.0
  end 
  if diodeon*(((u[2]+u[6])*Rs-u[10])/(Rd+Rs))+(1.0-diodeon)*(((u[2]+u[6])*Rs-u[10])*Rd/(Rd+Rs))>0
    Rd=ROn
    diodeon=1.0
  else
    Rd=ROff
    diodeon=0.0
  end 
  if diodeon*(((u[3]+u[7])*Rs-u[11])/(Rd+Rs))+(1.0-diodeon)*(((u[3]+u[7])*Rs-u[11])*Rd/(Rd+Rs))>0
    Rd=ROn
    diodeon=1.0
  else
    Rd=ROff
    diodeon=0.0
  end 
  if diodeon*(((u[4]+u[8])*Rs-u[12])/(Rd+Rs))+(1.0-diodeon)*(((u[4]+u[8])*Rs-u[12])*Rd/(Rd+Rs))>0
    Rd=ROn
    diodeon=1.0
  else
    Rd=ROff
    diodeon=0.0
  end 
end

tspan=(0.0,0.0005)
u = zeros(13)
p = [1e5,1e-5,1e-4,0.0,0.0]
odeprob=ODEProblem(cuk4,u,tspan,p)
sol= solve(odeprob,nmliqss2(),abstol=1e-4,reltol=1e-3)
@test  sol.stats.totalSteps<20000
@test -0.46<sol(0.0004,idxs=1)<-0.4
@test 0.4<sol(0.0004,idxs=5)<0.45
@test 2.55<sol(0.0004,idxs=9)<2.57
@test -2.98<sol(0.0004,idxs=13)<-2.9
