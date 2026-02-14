
#using BSON
using QuantizedSystemSolver
using Plots
#using BenchmarkTools

include("QSSsystems.jl")


alg= liqss2()
 abstol, reltol =  1e-4,1e-3

#= tspan=(0.0, 40.0) 
sol=get_Sol_LTI(alg,tspan ,abstol, reltol)
BSON.@load joinpath(@__DIR__,  "BSON_data", "ref_sol_LTI_Feagin12.bson") ref_sol_LTI_Feagin12
solInterp=solInterpolated(sol,0.01)
LTI_err=getAverageErrorByRefs(solInterp,ref_sol_LTI_Feagin12)
@show LTI_err =#


tspan=(0.0, 40.0)
sol=get_Sol_Bball(alg,tspan, abstol, reltol)
BSON.@load joinpath(@__DIR__,  "BSON_data", "ref_sol_Bball.bson") ref_sol_Bball
solInterp=solInterpolated(sol,0.01)
p= plot(solInterp, marker=:circle, title="Bball_system", markersize=1)
    savefig(p, "$(alg)_interp_Bball_system_tol_$(abstol)_$(tspan[2]).png")

#= for k in 1:1
    #p1=plot!(sol.t, [sol.u[i][k] for i in 1:length(sol.t)], label="solution",marker=:circle, markersize=1)
    p= plot(solInterp,idxs=[1], marker=:circle, title="Bball_system", markersize=1)
    p1=plot!(solInterp[1][1], [ref_sol_Bball[i][k] for i in 1:length(solInterp[1][1])], label="reference",marker=:star, markersize=1,
        title="Bball_system variable $k")
    savefig(p1, "$(alg)_interp_variable_$(k)_interleaved_tol_$(abstol)_$(tspan[2]).png")
end   =#

Bball_err=getAverageErrorByRefs(solInterp,ref_sol_Bball)
@show Bball_err


tspan = (0.0, 300)
sol=get_Sol_TsodyksMarkram(alg,tspan, abstol, reltol)
BSON.@load joinpath(@__DIR__,  "BSON_data", "ref_sol_TsodyksMarkram.bson") ref_sol_TsodyksMarkram
solInterp=solInterpolated(sol,0.1)  
TsodyksMarkram_err=getAverageErrorByRefs(solInterp,ref_sol_TsodyksMarkram)
@show TsodyksMarkram_err

tspan= (0.0, 0.0005)
sol=get_Sol_Interleaved(alg,tspan, abstol, reltol)  
BSON.@load joinpath(@__DIR__,  "BSON_data", "ref_sol_InterleavedRodas.bson") ref_sol_InterleavedRodas
solInterp=solInterpolated(sol,1e-6) 
Interleaved_err=getAverageErrorByRefs(solInterp,ref_sol_InterleavedRodas)
@show Interleaved_err


tspan=(0.0,300.0)
sol=get_Sol_IZH(alg,tspan, abstol, reltol)
BSON.@load joinpath(@__DIR__,  "BSON_data", "ref_sol_IZH.bson") ref_sol_IZH
solInterp=solInterpolated(sol,0.1)
IZH_err=getAverageErrorByRefs(solInterp,ref_sol_IZH)
@show IZH_err


tspan=(0.0,0.005)
sol=get_Sol_Railgun(alg,tspan, abstol, reltol) 
BSON.@load joinpath(@__DIR__,  "BSON_data", "ref_sol_Railgun.bson") ref_sol_Railgun
solInterp=solInterpolated(sol,1e-5)
Railgun_err=getAverageErrorByRefs(solInterp,ref_sol_Railgun)
@show Railgun_err




#= liqss2
LTI_err = 0.11170333156750877
Bball_err = 0.08852397591895964 

qss2
LTI_err = 0.11202475475322707
Bball_err = 0.08220989766101719
=#
