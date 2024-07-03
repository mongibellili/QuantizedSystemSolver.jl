using QuantizedSystemSolver



 #=    # Write your tests here.
    odeprob = NLodeProblem(quote
        #sys b53
        name=(sysb53,)
        u = [-1.0, -2.0]
        du[1] = -20.0*u[1]-80.0*u[2]+1600.0
        du[2] =1.24*u[1]-0.01*u[2]+0.2
    end)  
    tspan=(0.0,1.0)
    sol=solve(odeprob,nmliqss1(),tspan)
    save_Sol(sol)
    xp=sol(2,0.5)
    @test   xp  =#
    t1=Taylor0([1.0,1.0,0.0],2)
    t2=Taylor0([0.5,2.0,3.0],2)
    t3=Taylor0([2.0,1.0,0.0],2)
    cache1=Taylor0([0.0,0.0,0.0],2)
    cache2=Taylor0([0.0,0.0,0.0],2)
    cache3=Taylor0([0.0,0.0,0.0],2)
    #= @test   exp(t2,cache1)[0]≈1.6487212707001282
    @test   log(t2,cache1)[0]≈-0.6931471805599453
    @test   sin(t2,cache1,cache2)[0]≈0.479425538604203
    @test   cos(t2,cache1,cache2)[0]≈0.8775825618903728
    @test   tan(t2,cache1,cache2)[0]≈0.5463024898437905
    @test   asin(t2,cache1,cache2,cache3)[0]≈0.5235987755982989
    @test   acos(t2,cache1,cache2,cache3)[0]≈1.0471975511965979 =#
    @test   (t2^2)[0]≈0.25
    @test   (t2^3.0)[0]≈0.125
    @test   sqrt(t2)[0]≈0.7071067811865476
    @test   powerT(t2,2,cache1)[0]≈0.25
    @test   powerT(t2,3.0,cache1)[0]≈0.125
    @test   sqrt(t2,cache1)[0]≈0.7071067811865476