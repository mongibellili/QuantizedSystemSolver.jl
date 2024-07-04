using QuantizedSystemSolver

odeprob = NLodeProblem(quote
    name=(sysb1,)
    u = [-1.0, -2.0]
    du[1] = -2.0
    du[2] =1.24*u[1]-0.01*u[2]+0.2
end)  
tspan=(0.0,1.0)
sol=solve(odeprob,nmliqss2(),tspan)
save_Sol(sol)
xp=sol(2,0.5)
@show   xp  
odeprob = NLodeProblem(quote
    name=(sysb2,)
    u = [-1.0, -2.0]
    du[1] = -u[2]
    du[2] =1.24*u[1]-0.01*u[2]+0.2
end)  
tspan=(0.0,1.0)
sol=solve(odeprob,nmliqss2(),tspan)
save_Sol(sol)
xp=sol(2,0.5)
@show   xp  #xp = -2.1809630246611893 liqss1
    # Write your tests here.
    odeprob = NLodeProblem(quote
        #sys b53
        name=(sysb53,)
        u = [-1.0, -2.0]
        du[1] = -20.0*u[1]-80.0*u[2]+1600.0
        du[2] =1.24*u[1]-0.01*u[2]+0.2
    end)  
    tspan=(0.0,10.0)
    sol=solve(odeprob,qss1(),tspan)
   # sol=solve(odeprob,qss2(),tspan,abstol=1e-6,reltol=1e-3)
    save_Sol(sol)
    odeprob = NLodeProblem(quote
    #sys b53
    name=(sysb53,)
    u = [-1.0, -2.0]
    du[1] = -20.0*u[1]-80.0*u[2]+1600.0
    du[2] =1.24*u[1]-0.01*u[2]+0.2
end)  
tspan=(0.0,10.0)
sol=solve(odeprob,qss2(),tspan)
# sol=solve(odeprob,qss2(),tspan,abstol=1e-6,reltol=1e-3)
save_Sol(sol)
   #=  xp=sol(2,0.5)
    @show   xp  =#
   #=  t1=Taylor0([1.0,1.0,0.0],2)
    t2=Taylor0([0.5,2.0,3.0],2)
    t3=Taylor0([2.0,1.0,0.0],2)
    cache1=Taylor0([0.0,0.0,0.0],2)
    cache2=Taylor0([0.0,0.0,0.0],2)
    cache3=Taylor0([0.0,0.0,0.0],2) =#
    #= @test   exp(t2,cache1)[0]≈1.6487212707001282
    @test   log(t2,cache1)[0]≈-0.6931471805599453
    @test   sin(t2,cache1,cache2)[0]≈0.479425538604203
    @test   cos(t2,cache1,cache2)[0]≈0.8775825618903728
    @test   tan(t2,cache1,cache2)[0]≈0.5463024898437905
    @test   asin(t2,cache1,cache2,cache3)[0]≈0.5235987755982989
    @test   acos(t2,cache1,cache2,cache3)[0]≈1.0471975511965979 =#
  #=   @test   (t2^2)[0]≈0.25
    @test   (t2^3.0)[0]≈0.125
    @test   sqrt(t2)[0]≈0.7071067811865476
    @test   powerT(t2,2,cache1)[0]≈0.25
    @test   powerT(t2,3.0,cache1)[0]≈0.125
    @test   sqrt(t2,cache1)[0]≈0.7071067811865476 =#



#=     using Coverage
# process '*.cov' files
coverage = process_folder() # defaults to src/; alternatively, supply the folder name as argument
coverage = append!(coverage, process_folder("deps"))  # useful if you want to analyze more than just src/
# process '*.info' files, if you collected them
coverage = merge_coverage_counts(coverage, filter!(
    let prefixes = (joinpath(pwd(), "src", ""),
                    joinpath(pwd(), "deps", ""))
        c -> any(p -> startswith(c.filename, p), prefixes)
    end,
    LCOV.readfolder("test")))
# Get total coverage for all Julia files
covered_lines, total_lines = get_summary(coverage)
# Or process a single file
@show get_summary(process_file(joinpath("src", "QuantizedSystemSolver"))) =#

#= using Coverage
coverage = process_folder()
open("lcov.info", "w") do io
    LCOV.write(io, coverage)
end; =#