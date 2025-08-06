

function expExtern(x)
  exp(x)
end
function testExternal(du,u,p,t)
    du[1] = -expExtern(u[2])
    du[2] = (u[1])
end
tspan=(0.0,1.0)
u = [1.0, 0.0]
odeprob=ODEProblem(testExternal,u,tspan)
sol=solve(odeprob,qss2())
@test 0.07<sol(0.7,idxs=1)<0.1


function test()
  function testInnerScopeExternal(du,u,p,t)
    du[1] = -expExtern(u[2])
    du[2] = (u[1])
  end
  tspan=(0.0,1.0)
  u = [1.0, 0.0]
  odeprob=ODEProblem(testInnerScopeExternal,u,tspan)
  sol=solve(odeprob,qss2())
  @test 0.07<sol(0.7,idxs=1)<0.1
end
test()


function testExpoInternal(du,u,p,t)
  function expInter(x)
    exp(x)
  end
    du[1] = -expInter(u[2])
    du[2] = (u[1])
end
tspan=(0.0,1.0)
u = [1.0, 0.0]
odeprob=ODEProblem(testExpoInternal,u,tspan)
sol=solve(odeprob,qss2())
@test 0.07<sol(0.7,idxs=1)<0.1

# same as above for a discrete problem

function testExternalDiscrete(du,u,p,t)
    du[1] = -expExtern(u[2])
    du[2] = (u[1])
    if t>2.0
      u[1]=0.0
    end
end
tspan=(0.0,1.0)
u = [1.0, 0.0]
odeprob=ODEProblem(testExternalDiscrete,u,tspan)
sol=solve(odeprob,qss2())
@test 0.07<sol(0.7,idxs=1)<0.1


function test1()
  function testInnerScopeExternal(du,u,p,t)
    du[1] = -expExtern(u[2])
    du[2] = (u[1])
    if t>2.0
      u[1]=0.0
    end
  end
  tspan=(0.0,1.0)
  u = [1.0, 0.0]
  odeprob=ODEProblem(testInnerScopeExternal,u,tspan)
  sol=solve(odeprob,qss2())
  @test 0.07<sol(0.7,idxs=1)<0.1
end
test1()


function testExpoInternalDiscrete(du,u,p,t)
  function expInter(x)
    exp(x)
  end
    du[1] = -expInter(u[2])
    du[2] = (u[1])
     if t>2.0
      u[1]=0.0
    end
end
tspan=(0.0,1.0)
u = [1.0, 0.0]
odeprob=ODEProblem(testExpoInternalDiscrete,u,tspan)
sol=solve(odeprob,qss2())
@test 0.07<sol(0.7,idxs=1)<0.1

println("End external function tests.")

