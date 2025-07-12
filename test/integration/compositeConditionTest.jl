

function one_t_three_events(du,u,p,t) 
    du[1] = t-u[2]   
    du[2]=u[1]-u[2]*p[1] ;
    function constantNumber()
        return 0.2
    end

    if 5.0 < t || u[1]+u[2]>3.0
        u[2]=-0.5
    end

    if 4.0 < t || u[1]-u[2]>0.5
        u[2]=u[2]
    else
        u[1]=-0.4
    end

    if u[1]>2.0 && u[2]>2.0
        u[1] = 0.1
        u[2] = 0.1
    end

    if u[1]>2.0 && t>2.0
        u[1] = constantNumber()
        u[2] = constantNumber()
    else
        u[1] = 0.3
        u[2] = 0.3

    end 

     if u[1]>2.0 && u[2]>2.0 || u[2]<0.0
        u[1] = 1.0
        u[2] = 1.0
    end
end
tspan=(0.0,6.0)
u = [1.0, 1.0]
p=[0.5]
odeprob=ODEProblem(one_t_three_events,u,tspan,p)
#@show odeprob.eqs
sol=solve(odeprob,qss2())
@show sol.stats


@test 0.2<sol(0.7,idxs=1)<0.7
@test 0.9<sol(0.7,idxs=2)<1.3