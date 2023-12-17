using qss
using BenchmarkTools


function test()
  odeprob = @NLodeProblem begin
    name=(cuk,)
    C = 1e-4; L = 1e-4; R = 10.0;U = 24.0; T = 1e-4; DC = 0.25; ROn = 1e-5;ROff = 1e5;L1=1e-4;C1=1e-4
    discrete = [1e5,1e-5,1e-4,0.0,0.0]
    u = [0.0,0.0,0.0,0.0]
    rd=discrete[1];rs=discrete[2];nextT=discrete[3];lastT=discrete[4];diodeon=discrete[5]
    il=u[1] ;uc=u[3];il1=u[2] ;uc1=u[4]
    id=((il+il1)*rs-uc1)/(rd+rs) # diode's current
    du[1] =(-id*rd-uc)/L #il
    du[2]=(U-uc1-id*rd)/L1 #il1
    du[3]=(il-uc/R)/C  #uc
    du[4]=(id-il)/C1   #uc1
    if t-nextT>0.0 
        lastT=nextT
        nextT=nextT+T
        rs=ROn
    end
    if t-lastT-DC*T>0.0 
        rs=ROff
    end                          
    if diodeon*(id)+(1.0-diodeon)*(id*rd)>0
      rd=ROn
      diodeon=1.0
    else
      rd=ROff
      diodeon=0.0
    end  
  end
  tspan=(0.0,0.001)
  sol= solve(odeprob,nmliqss2(),abstol=1e-4,reltol=1e-3,tspan)
  save_Sol(sol,note="selectivecancel")
end
test()
