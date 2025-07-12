


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
 