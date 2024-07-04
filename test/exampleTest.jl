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
@test sol(2,0.0005)≈19.209921627620943
sol= solve(odeprob,nmliqss1(),tspan,abstol=1e-4,reltol=1e-3)    
@test sol(2,0.0005)≈19.20307419699851
sol= solve(odeprob,qss2(),tspan,abstol=1e-4,reltol=1e-3)  


BSON.@load "./ref_bson/solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01
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
getErrorByRefs(solFeagin14VectorN1000d01,solnmliqssInterp,1)
err4=getAverageErrorByRefs(solFeagin14VectorN1000d01,solnmliqssInterp)
@test err4≈0.00811021