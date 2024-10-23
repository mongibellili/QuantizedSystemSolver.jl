using DifferentialEquations
using Plots
function odeDiffEquPackage()
    function f(du, u, p, t)
        C = 1e-4; L = 1e-4; R = 10;U = 24.0; T = 1e-4; DC = 0.5; ROn = 1e-5;ROff = 1e5;
        rd=p[1];rs=p[2];
        il=u[1] ;uc=u[2]
        id=(il*rs-U)/(rd+rs) # diode's current
        du[1] =(-id*rd-uc)/L # inductor's current
        du[2]=(il-uc/R)/C    # capacitor's voltage
    end
    function condition1( u, t, integrator) 
        (t-p[3])
    end
    function condition2( u, t, integrator) 
        (t-p[4]-0.5*1e-4)
    end
    function condition3( u, t, integrator) 
         (p[5]*((u[1]*p[2]-24.0)/(p[1]+p[2]))+(1.0-p[5])*((u[1]*p[2]-24.0)*p[1]/(p[1]+p[2])))
    end
    function affect1!(integrator)
                p[4]=p[3]
                p[3]=p[3]+1e-4
                p[2]=1e-5
    end
    function affect2!(integrator)
                p[2]=1e5
    end
    function affect3!(integrator)
            p[1]=1e-5
            p[5]=1.0
    end
    function affect33!(integrator)
        p[1]=1e5
        p[5]=0.0
    end
    cb1 = ContinuousCallback(condition1, affect1!,nothing;  )
    cb2 = ContinuousCallback(condition2, affect2!,nothing; )
    cb3 = ContinuousCallback(condition3, affect3!,affect33!;  )
    cbs = CallbackSet(cb1, cb2,cb3)
    u0 = [0.0, 0.0]
    tspan = (0.0, 0.0002)
    p = [1e5,1e-5,1e-4,0.0,0.0,0.0]
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, ImplicitEuler(), callback = cbs, reltol=1e-2,abstol=1e-3#= ,dtmax=1e-12 =#,maxiters=Int(1e12)#= dt = 1e-3, adaptive = false =#)
    p1=plot!(sol);
    savefig(p1, "Rodas5P()_34_buck_ft0002_")
 end
 odeDiffEquPackage() 