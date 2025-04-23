

function sysb1(du,u,p,t)
    du[1] = -2.0
    du[2] =1.24*u[1]-0.01*u[2]+0.2
end
tspan=(0.0,1.0)
u = [-1.0, -2.0]
odeprob=ODEProblem(sysb1,u,tspan)
sol=solve(odeprob,qss1())
sol=solve(odeprob,qss2())
sol=solve(odeprob,liqss1())
sol=solve(odeprob,liqss2())
sol=solve(odeprob,nmliqss1())
sol=solve(odeprob,nmliqss2())
@test -2.6<sol(0.7,idxs=1)<-2.2
@test -3.517<sol(0.7,idxs=2)<-3.117

sol(0.5)
function sysb2(du,u,p,t)
    name=("sysb2",)
    u = [-1.0, -2.0]
    du[1] = -u[2]
    du[2] =1.24*u[1]+0.01*u[2]-0.2
end  
tspan=(0.0,1.0)
u = [-1.0, -2.0]
odeprob=ODEProblem(sysb2,u,tspan)
sol=solve(odeprob,qss1())
sol=solve(odeprob,qss2())
sol=solve(odeprob,liqss1())
sol=solve(odeprob,liqss2())
sol=solve(odeprob,nmliqss1())
sol=solve(odeprob,nmliqss2())
@test 0.5<sol(0.7,idxs=1)<0.7
@test -2.4<sol(0.7,idxs=2)<-2.2

function sysb53(du,u,p,t)
  du[1] = -20.0*u[1]-80.0*u[2]+1600.0
  du[2] =1.24*u[1]-0.01*u[2]+0.2
end
tspan=(0.0,1.0)
u = [-1.0, -2.0]
odeprob=ODEProblem(sysb53,u,tspan)
for N=0:6
  solve(odeprob,nmliqss2(),detection=Detection(N))
end
for N=0:13
  solve(odeprob,nmliqss1(),detection=Detection(N))
end


sol=solve(odeprob,nmliqss2(),reltol=1e-3,abstol=1e-5)

u1, u2 = -8.73522174738572, -7.385745994549763
λ1, λ2 = -10.841674966758294, -9.168325033241706
c1, c2 = 121.14809142478035, -143.14809142478035
xp1, xp2 = 0.0, 20.0
x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2


solnmliqssInterp=solInterpolated(sol,0.01)
er1=getError(solnmliqssInterp,1,x1)  
er2=getError(solnmliqssInterp,2,x2) 
@test 0.0<er1<0.01
@test 0.0<er2<0.01
avgErr=getAverageError(solnmliqssInterp,[x1,x2])
@test 0.0<avgErr<0.01

function buck(dy,y,p,t)# api requires four args
      
  C = 1e-4; L = 1e-4; R = 10.0;U = 24.0; T = 1e-4; DC = 0.5; ROn = 1e-5;ROff = 1e5;RR=0.1;LL=1e4;CC=1e4
  
  rd=p[1];rs=p[2];nextT=p[3];lastT=p[4];#diodeon=p[5]
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
p = [1e5,1e-5,1e-4,0.0,0.0];u0 = [0.0,0.0]
odeprob = ODEProblem(buck,u0,tspan,p)
sol= solve(odeprob,liqss2(),abstol=1e-4,reltol=1e-2) 
#save_Sol(sol)
@test 19.2<sol(0.0005,idxs=2)<19.9
sol= solve(odeprob,nmliqss2(),abstol=1e-3,reltol=1e-2)   
@test 19.2<sol(0.0005,idxs=2)<19.9

BSON.@load "solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01
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
sol=solve(odeprob,abstol=1e-5,reltol=1e-2,tspan)#
solnmliqssInterp=solInterpolated(sol,0.01)
getErrorByRefs(solnmliqssInterp,1,solFeagin14VectorN1000d01)
err4=getAverageErrorByRefs(solnmliqssInterp,solFeagin14VectorN1000d01)
@test err4<0.05
@test 0.35<sol(1.5,idxs=1)<0.39
@test 0.62<sol(1.5,idxs=2)<0.67
@test 0.92<sol(1.5,idxs=400)<1.0
@test 0.92<sol(1.5,idxs=600)<1.0
@test 0.92<sol(1.5,idxs=1000)<1.0


function tyson(du,u,p,t)
du[1] = u[4]-1e6*u[1]+1e3*u[2]
du[2] =-200.0*u[2]*u[5]+1e6*u[1]-1e3*u[2]
du[3] = 200.0*u[2]*u[5]-u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)
du[4] =u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)-u[4]
du[5] = 0.015-200.0*u[2]*u[5]
du[6] =u[4]-0.6*u[6]
end 
tspan=(0.0,25.0)
u = [0.0,0.75,0.25,0.0,0.0,0.0]
odeprob=ODEProblem(tyson,u,tspan)
sol=solve(odeprob,nmliqss2())
@test 0.00001<sol(20.0,idxs=1)<0.01
@test 0.8<sol(20.0,idxs=2)<0.99
@test 0.03<sol(20.0,idxs=3)<0.09
@test 0.00001<sol(20.0,idxs=4)<0.01
@test 0.0<sol(20.0,idxs=5)<8.0e-4
@test 0.01<sol(20.0,idxs=6)<0.03
solInterpolated(sol,1,0.01)
solnmliqssInterp=solInterpolated(sol,0.01)
BSON.@load "solRodas5PVectorTyson.bson" solRodas5PVectorTyson
err=getAverageErrorByRefs(solnmliqssInterp,solRodas5PVectorTyson) 
@test err<1.8
function sico(du,u,p,t)
    du[1] = -sin(t)
    du[2] = u[1]
end
tspan=(0.0,1.0)
u = [1.0, 0.0]
odeprob=ODEProblem(sico,u,tspan)
sol=solve(odeprob,qss2())
@test 0.5<sol(0.7,idxs=1)<0.9
@test 0.4<sol(0.7,idxs=2)<0.8
function expo(du,u,p,t)
    du[1] = -exp(u[2])
    du[2] = (u[1])
end
tspan=(0.0,1.0)
u = [1.0, 0.0]
odeprob=ODEProblem(expo,u,tspan)
sol=solve(odeprob,qss2())
@test 0.07<sol(0.7,idxs=1)<0.1
@test 0.2<sol(0.7,idxs=2)<0.6
function abso(du,u,p,t)
    du[1] = -abs(u[2])
    du[2] = (u[1])
end
tspan=(0.0,1.0)
u = [1.0, 0.0]
odeprob=ODEProblem(abso,u,tspan)
sol=solve(odeprob,qss1())
sol=solve(odeprob,liqss1())
@test 0.6<sol(0.7,idxs=1)<0.8
@test 0.5<sol(0.7,idxs=2)<0.7
function acossin(du,u,p,t)
    u = [1.0, 0.0]
    du[1] = acos(sin(u[2]))
    du[2] = (u[1])
end 
tspan=(0.0,1.0)
u = [1.0, 0.0]
odeprob=ODEProblem(acossin,u,tspan)
sol=solve(odeprob,qss2())
@test 1.3<sol(0.5,idxs=1)<1.8
@test 0.4<sol(0.5,idxs=2)<0.8
function squar(du,u,p,t)
    name=("sysbN7",)
    u = [1.0, 0.0]
    du[1] = -(u[2]+t^2)
    du[2] = u[1]
end
tspan=(0.0,1.0)
u = [1.0, 0.0]
odeprob=ODEProblem(squar,u,tspan)
sol=solve(odeprob,qss2())
@test 0.7<sol(0.5,idxs=1)<0.9
 
function onet(du,u,p,t)
    du[1] = t
    du[2] =1.24*u[1]-0.01*u[2]+0.2
end 
tspan=(0.0,5.0)
u = [1.0, 0.0]
odeprob=ODEProblem(onet,u,tspan)
sol=solve(odeprob,qss1())
sol=solve(odeprob,qss2())
sol=solve(odeprob,liqss1())
sol=solve(odeprob,liqss2())
sol=solve(odeprob,nmliqss1())
sol=solve(odeprob,nmliqss2())
@test 0.6<sol(0.5,idxs=2)<0.8
 

 function oneteven(du,u,p,t)
du[1] = t
du[2] =1.24*u[1]-0.01*u[2]+0.2
if t-2.0>0.0
    u[1] = 0.0
else
    u[2] = 0.0
end
end 
tspan=(0.0,5.0)
u = [1.0, 0.0]
p=[]
odeprob=ODEProblem(oneteven,u,tspan,p)
sol=solve(odeprob,qss2())
sol=solve(odeprob,liqss2())
sol=solve(odeprob,nmliqss2())
@test 0.6<sol(0.5,idxs=2)<0.8

function onetevloop(du,u,p,t)
    
    du[1] = t
    
    for k in 2:5 
        du[k]=(u[k]-u[k-1]) ;
    end 
    if t-3.0>0.0
        u[1] = 1.0
        u[2] = 0.0
        u[3] = 1.0
        u[4] = 0.0
        u[5] = 1.0
        
    end
end  
tspan=(0.0,6.0)
u = [1.0, 0.0,1.0, 0.0,1.0]
p=[]
odeprob=ODEProblem(onetevloop,u,tspan,p)
sol=solve(odeprob,qss2())
sol=solve(odeprob,liqss2())
sol=solve(odeprob,nmliqss2())
@test -0.75<sol(0.5,idxs=2)<-0.6
plot_SolSum(sol)
plot_SolSum(sol,1,2,xlims=(0.0,1.0),ylims=(0.0,2.0))
plot_SolSum(sol,1,2,xlims=(0.0,1.0),ylims=(0.0,0.0))
plot_SolSum(sol,1,2,xlims=(0.0,0.0),ylims=(0.0,1.0))
plot_SolSum(sol,1,2)
#= save_Sol(sol,1,2,3,4,5)
save_SolSum(sol,1,2)  =#

plot(sol,idxs=[1],xlims=(0.0,1.0),ylims=(0.0,2.0))
plot(sol,idxs=[1],xlims=(0.0,1.0),ylims=(0.0,0.0))
plot(sol,idxs=[1],xlims=(0.0,0.0),ylims=(0.0,1.0))
plot(sol)
plot(sol,idxs=[1,2,3])
plot(sol,idxs=(1,2))
plot(sol,idxs=(1,2,3))

function onetevsloop(du,u,p,t)
  
    du[1] = t
    
    for k in 2:5 
        du[k]=p[1]*(u[k]-u[k-1]) ;
    end 
    if t-5.0>0.0
        p[1]=0.0
    end
    if t-3.0>0.0
        u[1] = 1.0
        u[2] = 0.0
        u[3] = 1.0
        u[4] = 0.0
        u[5] = 1.0
        p[1]=1.0
        
    end
end 
tspan=(0.0,6.0)
u = [1.0, 0.0,1.0, 0.0,1.0]
p=[0.5]
odeprob=ODEProblem(onetevsloop,u,tspan,p)
sol=solve(odeprob,qss1())
sol=solve(odeprob,liqss1())
sol=solve(odeprob,nmliqss2())
@test 1.1<sol(0.5,idxs=1)<1.3
@test -0.35<sol(0.5,idxs=2)<-0.2


function onetevsloop2(du,u,p,t)
    name=("sysN14",)
    u = [1.0, 0.0,1.0, 0.0,1.0]
    p = [0.5]
    du[1] = t
    
    for k in 2:5 
        du[k]=p[1]*(u[k]-u[k-1]) ;
    end 
    if t>5.0
        p[1]=0.0
    end
    if u[1]-2.0>0.0
        u[1] = 1.0
        u[2] = 0.0
        u[3] = 1.0
        u[4] = 0.0
        u[5] = 1.0
        p[1]=1.0
        
    end
end
tspan=(0.0,6.0)
u = [1.0, 0.0,1.0, 0.0,1.0]
p=[0.5]
odeprob=ODEProblem(onetevsloop2,u,tspan,p)
sol=solve(odeprob,nmliqss2())
@test -0.35<sol(0.5,idxs=2)<-0.2
sol=solve(odeprob,liqss2())



function onetevsloop3(du,u,p,t)
    du[1] = t+u[2]
    
    for k in 2:5 
        du[k]=p[k]*(u[k*1]-u[k-1]-u[k+1])+(p[k+1]+p[k-1])*p[k*1] ;
    end 
    if t-5.0>0.0
        p[2]=0.0
    end
    if u[1]+u[2]>3.0
        u[1] = 1.0
        u[3] = 1.0
        u[4] = 0.0
        p[1]=1.0
    end
    if p[1]-u[4]+u[5]>0.0
        u[2]=0.0
    end
end
tspan=(0.0,6.0)
u = [1.0, 1.0,1.0, 0.0,1.0,1.0]
p=[0.5,1.0,1.0,1.0,1.0,1.0]
odeprob=ODEProblem(onetevsloop3,u,tspan,p)
sol=solve(odeprob,qss2())
sol=solve(odeprob,liqss2())
sol=solve(odeprob,nmliqss2())
@test 1.7<sol(0.7,idxs=1)<1.95
@test 0.2<sol(0.7,idxs=2)<0.35

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
@test -0.46<sol(0.0004,idxs=1)<-0.4
@test 0.4<sol(0.0004,idxs=5)<0.45
@test 2.55<sol(0.0004,idxs=9)<2.57
@test -2.98<sol(0.0004,idxs=13)<-2.9



function Power_Grid(du, u, p, t)
  # System size and parameters
  N = 5  
  C=1
  M = [1.0, 1.2, 0.9, 1.1, 1.0] 
  D = [0.2, 0.3, 0.25, 0.22, 0.2] 
  α, β= 12.0, 10.001
  a=u[2]
    # Differential equations   
  for k in C:N-1
      du[k] = β*u[k+N] / M[k] 
  end
  du[N-1+C] = β*u[N+N] / M[N] 
  for k in 6:10
      du[k] = - a*α*u[k] - D[k-5] * u[k] 
  end
end

# Initial conditions
u0 = [0.1, 0.2, 0.15, 0.05, 0.1,5.0, 1.0, 0.5, 2.0, 0.0]

# Time span
tspan = (0.0, 5.0)

# odeprob
odeprob = ODEProblem(Power_Grid, u0, tspan)

