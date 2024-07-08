odeprob = NLodeProblem(quote
    name=(sysb1,)
    u = [-1.0, -2.0]
    du[1] = -2.0
    du[2] =1.24*u[1]-0.01*u[2]+0.2
end)  
tspan=(0.0,1.0)
sol=solve(odeprob,qss1(),tspan)
sol=solve(odeprob,qss2(),tspan)
sol=solve(odeprob,liqss1(),tspan)
sol=solve(odeprob,liqss2(),tspan)
sol=solve(odeprob,nmliqss1(),tspan)
sol=solve(odeprob,nmliqss2(),tspan)

odeprob = NLodeProblem(quote
    name=(sysb2,)
    u = [-1.0, -2.0]
    du[1] = -u[2]
    du[2] =1.24*u[1]+0.01*u[2]-0.2
end)  
tspan=(0.0,1.0)
sol=solve(odeprob,qss1(),tspan)
sol=solve(odeprob,qss2(),tspan)
sol=solve(odeprob,liqss1(),tspan)
sol=solve(odeprob,liqss2(),tspan)
sol=solve(odeprob,nmliqss1(),tspan)
sol=solve(odeprob,nmliqss2(),tspan)

odeprob = NLodeProblem(quote
    name=(sysb53,)
    u = [-1.0, -2.0]
    du[1] = -20.0*u[1]-80.0*u[2]+1600.0
    du[2] =1.24*u[1]-0.01*u[2]+0.2
end)  
tspan=(0.0,10.0)
sol=solve(odeprob,qss1(),tspan)
sol=solve(odeprob,qss2(),tspan)
sol=solve(odeprob,liqss1(),tspan)
sol=solve(odeprob,liqss2(),tspan)
sol=solve(odeprob,nmliqss1(),tspan)
sol=solve(odeprob,nmliqss2(),tspan)
plot_Sol(sol)
plot_Sol(sol,1)
getPlot(sol)

u1, u2 = -8.73522174738572, -7.385745994549763
λ1, λ2 = -10.841674966758294, -9.168325033241706
c1, c2 = 121.14809142478035, -143.14809142478035
xp1, xp2 = 0.0, 20.0
x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2


solnmliqssInterp=solInterpolated(sol,0.01)
er1=getError(solnmliqssInterp,1,x1)  
er2=getError(solnmliqssInterp,2,x2) 
@test 0.0009<er1<0.00099
@test 7.0e-5<er2<7.4e-5

odeprob = NLodeProblem(quote   #NLodeProblem(quote ... end);
name=(buck,)
C = 1e-4; L = 1e-4; R = 10.0;U = 24.0; T = 1e-4; DC = 0.5; ROn = 1e-5;ROff = 1e5;
discrete = [1e5,1e-5,1e-4,0.0,0.0];u = [0.0,0.0]
rd=discrete[1];rs=discrete[2];nextT=discrete[3];lastT=discrete[4];diodeon=discrete[5]
il=u[1] ;uc=u[2]
id=(il*rs-U)/(rd+rs) # diode's current
du[1] =(-id*rd-uc)/L
du[2]=(il-uc/R)/C
if t-nextT>0.0 
  lastT=nextT;nextT=nextT+T;rs=ROn
end
if t-lastT-DC*T>0.0 
  rs=ROff
end                          
if diodeon*(id)+(1.0-diodeon)*(id*rd-0.6)>0
  rd=ROn;diodeon=1.0
else
  rd=ROff;diodeon=0.0
end     
end)
tspan = (0.0, 0.001)
sol= solve(odeprob,nmliqss2(),tspan,abstol=1e-4,reltol=1e-3)    
@test 19.0<sol(2,0.0005)<19.4
sol= solve(odeprob,nmliqss1(),tspan,abstol=1e-4,reltol=1e-3)    
@test 19.0<sol(2,0.0005)<19.4
sol= solve(odeprob,qss2(),tspan,abstol=1e-4,reltol=1e-3)  


BSON.@load "./solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01
prob=NLodeProblem(quote
  name=(adrN1000d01,)
  u[1:333]=1.0
  u[334:1000]=0.0
  _dx=100.0#1/dx=N/10=1000/10
  a=1.0;d=0.1;r=1000.0
  #discrete=[0.0]
  du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
  for k in 2:999  
    du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
  end 
  du[1000]=-a*_dx*(u[1000]-u[999])+d*_dx*_dx*(2.0*u[999]-2.0*u[1000])+r*u[1000]*u[1000]*(1.0-u[1000]) 
  #=  if u[1]-10.0>0.0 #fake to test discreteintgrator & loop
  discrete[1]=1.0
  end =#
end)



tspan=(0.0,5.0)
solnmliqss=solve(prob,abstol=1e-5,reltol=1e-2,tspan)#
solnmliqssInterp=solInterpolated(solnmliqss,0.01)
getErrorByRefs(solnmliqssInterp,1,solFeagin14VectorN1000d01)
err4=getAverageErrorByRefs(solnmliqssInterp,solFeagin14VectorN1000d01)

odeprob = NLodeProblem(quote
    name=(sysbN5,)
    u = [1.0, 0.0]
    du[1] = -sin(t)
    du[2] = (u[1])
end)  
tspan=(0.0,1.0)
sol=solve(odeprob,qss2(),tspan)
odeprob = NLodeProblem(quote
    name=(sysbN6,)
    u = [1.0, 0.0]
    du[1] = -exp(u[2])
    du[2] = (u[1])
end)  
tspan=(0.0,1.0)
sol=solve(odeprob,qss2(),tspan)
odeprob = NLodeProblem(quote
    name=(sysbN6,)
    u = [1.0, 0.0]
    du[1] = -abs(u[2])
    du[2] = (u[1])
end)  
tspan=(0.0,1.0)
sol=solve(odeprob,qss1(),tspan)
odeprob = NLodeProblem(quote
    name=(sysbN8,)
    u = [1.0, 0.0]
    du[1] = acos(sin(u[2]))
    du[2] = (u[1])
end)  
tspan=(0.0,1.0)
sol=solve(odeprob,qss2(),tspan)
odeprob = NLodeProblem(quote
    name=(sysbN7,)
    u = [1.0, 0.0]
    du[1] = -(u[2]+t^2)
    du[2] = u[1]
end)  
tspan=(0.0,1.0)
sol=solve(odeprob,qss2(),tspan)
@test 0.7<sol(1,0.5)<0.9

odeprob = NLodeProblem(quote
    name=(sysN10,)
    u = [1.0, 0.0]
    du[1] = t
    du[2] =1.24*u[1]-0.01*u[2]+0.2
end)  
tspan=(0.0,5.0)
sol=solve(odeprob,qss2(),tspan)
sol=solve(odeprob,liqss2(),tspan)
sol=solve(odeprob,nmliqss2(),tspan)
@test 0.6<sol(2,0.5)<0.8


odeprob = NLodeProblem(quote
    name=(sysN11,)
u = [1.0, 0.0]
du[1] = t
du[2] =1.24*u[1]-0.01*u[2]+0.2
if t-2.0>0.0
    u[1] = 0.0
    u[2] = 0.0
end
end)  
tspan=(0.0,5.0)
sol=solve(odeprob,qss2(),tspan)
sol=solve(odeprob,liqss2(),tspan)
sol=solve(odeprob,nmliqss2(),tspan)
@test 0.6<sol(2,0.5)<0.8

odeprob = NLodeProblem(quote
    name=(sysN12,)
    u = [1.0, 0.0,1.0, 0.0,1.0]
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
end)  
tspan=(0.0,6.0)
sol=solve(odeprob,qss2(),tspan)
sol=solve(odeprob,liqss2(),tspan)
sol=solve(odeprob,nmliqss2(),tspan)
getPlot(sol,1)
plot_SolSum(sol,1,2)
@test -0.75<sol(2,0.5)<-0.6


odeprob = NLodeProblem(quote
    name=(sysN13,)
    u = [1.0, 0.0,1.0, 0.0,1.0]
    discrete = [0.5]
    du[1] = t
    
    for k in 2:5 
        du[k]=discrete[1]*(u[k]-u[k-1]) ;
    end 
    if t-5.0>0.0
        discrete[1]=0.0
    end
    if t-3.0>0.0
        u[1] = 1.0
        u[2] = 0.0
        u[3] = 1.0
        u[4] = 0.0
        u[5] = 1.0
        discrete[1]=1.0
        
    end
end)  
tspan=(0.0,6.0)
sol=solve(odeprob,qss1(),tspan)
sol=solve(odeprob,liqss1(),tspan)
sol=solve(odeprob,nmliqss2(),tspan)
#@test -0.35<sol(2,0.5)<-0.2


odeprob = NLodeProblem(quote
    name=(sysN14,)
    u = [1.0, 0.0,1.0, 0.0,1.0]
    discrete = [0.5]
    du[1] = t
    
    for k in 2:5 
        du[k]=discrete[1]*(u[k]-u[k-1]) ;
    end 
    if t-5.0>0.0
        discrete[1]=0.0
    end
    if u[1]-3.0>0.0
        u[1] = 1.0
        u[2] = 0.0
        u[3] = 1.0
        u[4] = 0.0
        u[5] = 1.0
        discrete[1]=1.0
        
    end
end)  
tspan=(0.0,6.0)
sol=solve(odeprob,nmliqss2(),tspan)
@test -0.35<sol(2,0.5)<-0.2