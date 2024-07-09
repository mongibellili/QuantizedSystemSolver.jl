using QuantizedSystemSolver



#= odeprob = NLodeProblem(quote
    name=(sysb53,)
    u = [-1.0, -2.0]
    du[1] = -20.0*u[1]-80.0*u[2]+1600.0
    du[2] =1.24*u[1]-0.01*u[2]+0.2
end)  
tspan=(0.0,10.0)

sol=solve(odeprob,nmliqss2(),tspan)

u1, u2 = -8.73522174738572, -7.385745994549763
λ1, λ2 = -10.841674966758294, -9.168325033241706
c1, c2 = 121.14809142478035, -143.14809142478035
xp1, xp2 = 0.0, 20.0
x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2


solnmliqssInterp=solInterpolated(sol,0.01)
er1=getError(solnmliqssInterp,1,x1)  
er2=getError(solnmliqssInterp,2,x2) 
@show er1,er2 =#

#= import Base:  *, /


/(a::Int,b::Int)=Int(Float64(a)/b)
a=5;b=1
@show a*b
@show a/b =#

#= using Plots
odeprob = NLodeProblem(quote
    name=(sysN14,)
    u = [1.0, 0.0,1.0, 0.0,1.0,1.0]
    discrete = [0.5,1.0,1.0,1.0,1.0,1.0]
    du[1] = t+u[2]
    
    for k in 2:5 
        du[k]=discrete[k]*(u[k*1]-u[k-1]-u[k+1])+(discrete[k+1]+discrete[k-1])*discrete[k*1] ;
    end 
    if t-5.0>0.0
        discrete[2]=0.0
    end
    if u[1]-3.0>0.0
        u[1] = 1.0
        u[2] = 0.0
        u[3] = 1.0
        u[4] = 0.0
        u[5] = 1.0
        discrete[1]=1.0
    end
    if discrete[1]-5.0>0.0
        u[2]=0.0
    end
end)  
tspan=(0.0,6.0)
sol=solve(odeprob,nmliqss2(),tspan)
p1=plot_SolSum(sol,1,2,xlims=(0.0,5.0),ylims=(-5.0,20.0))
 savefig(p1, "plot1")
p2=getPlot(sol,1,xlims=(0.0,5.0),ylims=(-6.0,20.0))
savefig(p2, "plot2") =#
#= save_Sol(sol)
xp=sol(2,0.5)
@show   xp   =#
#= 
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
@test   xp  #xp = -2.1809630246611893 liqss1
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
save_Sol(sol) =#
   #=  xp=sol(2,0.5)
    @test   xp  =#
   #=  t1=Taylor0([1.0,1.0,0.0],2)
    t2=Taylor0([0.5,2.0,3.0],2)
    t3=Taylor0([2.0,1.0,0.0],2)
    cache1=Taylor0([0.0,0.0,0.0],2)
    cache2=Taylor0([0.0,0.0,0.0],2)
    cache3=Taylor0([0.0,0.0,0.0],2) =#
    #= @test    exp(t2,cache1)[0]≈1.6487212707001282
    @test    log(t2,cache1)[0]#≈-0.6931471805599453
    @test    sin(t2,cache1,cache2)[0]≈0.479425538604203
    @test    cos(t2,cache1,cache2)[0]≈0.8775825618903728
    @test    tan(t2,cache1,cache2)[0]≈0.5463024898437905
    @test    asin(t2,cache1,cache2,cache3)[0]≈0.5235987755982989
    @test    acos(t2,cache1,cache2,cache3)[0]≈1.0471975511965979 =#
  #=   @test    (t2^2)[0]≈0.25
    @test    (t2^3.0)[0]≈0.125
    @test    sqrt(t2)[0]≈0.7071067811865476
    @test    powerT(t2,2,cache1)[0]≈0.25
    @test    powerT(t2,3.0,cache1)[0]≈0.125
    @test    sqrt(t2,cache1)[0]≈0.7071067811865476 =#



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
@test get_summary(process_file(joinpath("src", "QuantizedSystemSolver"))) =#
#= using Pkg
Pkg.test("QuantizedSystemSolver"; coverage=true) # run this then comment and run next code =#
#= using Coverage
coverage = process_folder()
open("lcov.info", "w") do io
    LCOV.write(io, coverage)
end;
 =#



#= t2=1.0

cache1=Taylor0([1.0,1.0,1.0],2)
cache2=Taylor0([1.0,1.0,1.0],2)
cache3=Taylor0([0.0,0.0,0.0],2)
@test    exp(t2,cache1)[0]≈2.718281828459045
@test    log(t2,cache1)[0]≈0.0
@test    sin(t2,cache1,cache2)[0]≈0.8414709848078965
@test    cos(t2,cache1,cache2)[0]≈0.5403023058681398
@test    tan(t2,cache1,cache2)[0]≈1.5574077246549023
 
@test    abs(t2,cache1)[0]≈1.0


@test    (t2^2)≈1.0
@test    (t2^3.0)≈1.0
@test    sqrt(t2)≈1.0
@test    powerT(t2,2,cache1)[0]≈1.0
@test    powerT(t2,3.0,cache1)[0]≈1.0
@test    sqrt(t2,cache1)[0]≈1.0 =#


#= using Test
x=3.0*6.0+1.5
@test 19.0<x<19.4 =#
#= jac=[Int64[], [1]]
for i = 1:2
    @test !isempty(jac[i])
end =#


#= t2=1.0
    cache1=Taylor0([0.0,0.0,0.0],2)
    cache2=Taylor0([0.0,0.0,0.0],2)
@show one(cache1) =#
#= using Test
acceptedi=Vector{Vector{Float64}}(undef,3)
for i =1:3
acceptedi[i]=[0.0,0.0]#zeros(2)
end

constructIntrval(acceptedi,-1.0,-2.0,-3.0,4.0)
@test acceptedi[1]==[0.0, 4.0]

constructIntrval(acceptedi,-1.0,2.0,3.0,4.0)
@test acceptedi[1]==[0.0, 2.0]
@test acceptedi[2]==[3.0, 4.0]
constructIntrval(acceptedi,1.0,2.0,3.0,4.0)
@test acceptedi[1]  ==[0.0, 1.0]
@test acceptedi[2]==[2.0, 3.0]
@test acceptedi[3] ==[4.0, Inf] =#

@show 1==1.0