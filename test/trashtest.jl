#using QuantizedSystemSolver


#= 
    # Write your tests here.
    odeprob = NLodeProblem(quote
        #sys b53
        name=(sysb53,)
        u = [-1.0, -2.0]
        du[1] = -20.0*u[1]-80.0*u[2]+1600.0
        du[2] =1.24*u[1]-0.01*u[2]+0.2
    end)  
    tspan=(0.0,1.0)
    sol=solve(odeprob,nmliqss1(),tspan)
    #save_Sol(sol)
    xp=sol(2,0.5)
    @show xp =#
 #=    cache1=Taylor0([1.0,1.0,1.0],2)
    addT(4.1,5.8,cache1)
    @show cache1[0] =#
    @show 0.1+1.8