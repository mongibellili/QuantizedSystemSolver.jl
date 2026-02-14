
using BSON
using DifferentialEquations
using Plots

include("DEsystems.jl")



 

#alg=  Tsit5()
alg=Feagin12()
abstol, reltol =  1e-12,1e-8



#= tspan=(0.0, 40.0)
sol=get_Sol_LTI(alg, tspan, abstol, reltol)
ref_sol_LTI_Feagin12=sol.u
BSON.@save joinpath(@__DIR__,  "BSON_data", "ref_sol_LTI_Feagin12.bson") ref_sol_LTI_Feagin12=#


#= tspan=(0.0, 40.0)
sol=get_Sol_Bball(alg, tspan, abstol, reltol)
ref_sol_Bball=sol.u
BSON.@save joinpath(@__DIR__,  "BSON_data", "ref_sol_Bball.bson") ref_sol_Bball  =#

#= tspan = (0.0, 300)
sol=get_Sol_TsodyksMarkram(alg, tspan, abstol, reltol)
ref_sol_TsodyksMarkram=sol.u        
@show length(ref_sol_TsodyksMarkram)
BSON.@save joinpath(@__DIR__,  "BSON_data", "ref_sol_TsodyksMarkram.bson") ref_sol_TsodyksMarkram
 
 =#

#=   alg=Rodas5()
  abstol, reltol =  1e-7,1e-5
  tspan=(0.0,0.0005)
sol=get_Sol_Interleaved(alg, tspan, abstol, reltol)
ref_sol_InterleavedRodas=sol.u  
@show length(ref_sol_InterleavedRodas) 
BSON.@save joinpath(@__DIR__,  "BSON_data", "ref_sol_InterleavedRodas.bson") ref_sol_InterleavedRodas
 =#

#=  alg=Tsit5()
 abstol, reltol =  1e-7,1e-5 =#
#=     tspan=(0.0,300.0)
sol=get_Sol_IZH(alg,tspan, abstol, reltol)
ref_sol_IZH=sol.u   
@show length(ref_sol_IZH)
BSON.@save joinpath(@__DIR__,  "BSON_data", "ref_sol_IZH.bson") ref_sol_IZH  =#
#= alg=Rodas5()
 abstol, reltol =  1e-7,1e-5
    tspan=(0.0,0.005)
sol=get_Sol_Railgun(alg,tspan, abstol, reltol)
ref_sol_Railgun_compact=sol.u   
@show length(ref_sol_Railgun_compact)
BSON.@save joinpath(@__DIR__,  "BSON_data", "ref_sol_Railgun_compact.bson") ref_sol_Railgun_compact

 =#