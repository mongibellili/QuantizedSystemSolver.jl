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
    function condition(out, u, t, integrator) # Event when condition(out,u,t,integrator) == 0
        out[1] = (t-p[3])
        out[2] =  (t-p[4]-0.5*1e-4)
        out[3] = (p[5]*((u[1]*p[2]-24.0)/(p[1]+p[2]))+(1.0-p[5])*((u[1]*p[2]-24.0)*p[1]/(p[1]+p[2])))
        out[4] = -(p[5]*((u[1]*p[2]-24.0)/(p[1]+p[2]))+(1.0-p[5])*((u[1]*p[2]-24.0)*p[1]/(p[1]+p[2])))
    end
  
    function affect!(integrator, idx)
        if idx == 1
            p[4]=p[3]
            p[3]=p[3]+1e-4
            p[2]=1e-5
        elseif idx == 2
            p[2]=1e5
        elseif idx == 3
            p[1]=1e-5
            p[5]=1.0
        elseif idx == 4
            p[1]=1e5
            p[5]=0.0
        end
    end

  
    cbs = VectorContinuousCallback(condition, affect!, 4)
    u0 = [0.0, 0.0]
    tspan = (0.0, 0.0002)
    p = [1e5,1e-5,1e-4,0.0,0.0,0.0]
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, Rodas5P(), callback = cbs, reltol=1e-6,abstol=1e-9,dtmax=1e-12,maxiters=Int(1e12)#= dt = 1e-3, adaptive = false =#)
    p1=plot!(sol);
    savefig(p1, "Rodas5P()_vector_buck_ft0002_")
 end
 odeDiffEquPackage() 