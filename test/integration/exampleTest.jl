
# test du= constant; test all algs; test evaluations
function system1(du,u,p,t)
    du[1] = -2.0
    du[2] =1.24*u[1]-0.01*u[2]+0.2
end
tspan=(0.0,1.0)
u = [-1.0, -2.0]
odeprob=ODEProblem(system1,u,tspan)
sol=solve(odeprob,qss1())
sol=solve(odeprob,qss2())
sol=solve(odeprob,liqss1())
sol=solve(odeprob,liqss2())
sol=solve(odeprob,nmliqss1())
sol=solve(odeprob,nmliqss2())
@show sol(0.7)
@test -2.6<sol(0.7,idxs=1)<-2.2
@test -3.517<sol(0.7,idxs=2)<-3.117


# test du= one var; test all algs; test manual jacobian computation during simulation run
function system2(du,u,p,t)
    u = [-1.0, -2.0]
    du[1] = u[2]
    du[2] =1.24*u[1]+0.01*u[2]-0.2
end  
tspan=(0.0,1.0)
u = [-1.0, -2.0]
odeprob=ODEProblem(system2,u,tspan,jac_mode= :approximate)
sol=solve(odeprob,qss1())
sol=solve(odeprob,qss2())
sol=solve(odeprob,liqss1())
sol=solve(odeprob,liqss2())
sol=solve(odeprob,nmliqss1())
sol=solve(odeprob,nmliqss2())
@test -3.0<sol(0.7,idxs=1)<-2.8

# test cycle detections; test interpolation; test getError
function system3(du,u,p,t)
  du[1] = -20.0*u[1]-80.0*u[2]+1600.0
  du[2] =1.24*u[1]-0.01*u[2]+0.2
end
tspan=(0.0,1.0)
u = [-1.0, -2.0]
odeprob=ODEProblem(system3,u,tspan)
for N=0:6
  solve(odeprob,nmliqss2(),detection=Detection(N))
end
for N=0:12
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

# test du = t
function one_t(du,u,p,t)
    du[1] = t 
    du[2] =1.24*u[1]-0.01*u[2]+0.2
end 
tspan=(0.0,5.0)
u = [1.0, 0.0]
odeprob=ODEProblem(one_t,u,tspan)
sol=solve(odeprob,qss1())
sol=solve(odeprob,qss2())
sol=solve(odeprob,liqss1())
sol=solve(odeprob,liqss2())
sol=solve(odeprob,nmliqss1())
sol=solve(odeprob,nmliqss2())
@test 0.6<sol(0.5,idxs=2)<0.8
 
# test du = t for a discrete problem; test p=[]
function one_t_one_event(du,u,p,t)
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
odeprob=ODEProblem(one_t_one_event,u,tspan,p)
sol=solve(odeprob,qss2())
sol=solve(odeprob,liqss2())
sol=solve(odeprob,nmliqss2())
@test 0.6<sol(0.5,idxs=2)<0.8

# test du = t and a loop in a discrete problem; test plot
function one_t_discrete_loop(du,u,p,t)
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
odeprob=ODEProblem(one_t_discrete_loop,u,tspan,p)
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

# previous test + another event; order1 algs
function one_t_two_events(du,u,p,t)
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
odeprob=ODEProblem(one_t_two_events,u,tspan,p)
sol=solve(odeprob,qss1())
sol=solve(odeprob,liqss1())
sol=solve(odeprob,nmliqss2())
@test 1.1<sol(0.5,idxs=1)<1.3
@test -0.35<sol(0.5,idxs=2)<-0.2

# previous test + another event; order2 algs
function one_t_three_events(du,u,p,t)
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
odeprob=ODEProblem(one_t_three_events,u,tspan,p)
sol=solve(odeprob,qss2())
sol=solve(odeprob,liqss2())
sol=solve(odeprob,nmliqss2())
@test 1.7<sol(0.7,idxs=1)<1.95
@test 0.2<sol(0.7,idxs=2)<0.35


function testParams(du, u, p, t)
  # System size and parameters
  N = 5  
  @show N
  C=1
  M = [1.0, 1.2, 0.9, 1.1, 1.0] 
  D = [0.2, 0.3, 0.25, 0.22, 0.2] 
  α, β= 12.0, -1.001
  a=u[2] 
  for k in C:N-1
      du[k] = β
  end
  du[N-1+C] = β*u[N+N] / M[N] 
  for k in 6:10
      du[k] = - a*α*u[k] - D[k-5] * u[k] 
  end
  if t>2.0
    u[2]=0.0
  end
end
u0 = [0.1, 0.2, 0.15, 0.05, 0.1,5.0, 1.0, 0.5, 2.0, 0.0]
tspan = (0.0, 5.0)
odeprob = ODEProblem(testParams, u0, tspan)



function simple_loop(du,u,p,t)
    du[1] = u[3]
    for k in 2:3
        du[k]=u[k-1] ;
    end 
end  

tspan=(0.0,2.0)
u = [-0.5, 0.0,1.0] 
odeprob=ODEProblem(simple_loop,u,tspan,jac_mode= :approximate)

sol=solve(odeprob,nmliqss2())
@test 0.1<sol(0.7,idxs=1)<0.3

function simple_loop_discrete(du,u,p,t)
    du[1] = u[3]
    for k in 2:3
        du[k]=u[k-1] ;
    end 
    if u[1] > 0.0
        @show t
        m=1.0+u[2]
        u[1] = m
    end
end  

tspan=(0.0,2.0)
u = [-0.5, 0.0,1.0] 
odeprob=ODEProblem(simple_loop_discrete,u,tspan,jac_mode= :approximate)

sol=solve(odeprob,nmliqss2())
@show sol.stats
@test 0.9<sol(0.7,idxs=1)<1.1