
#using BSON
using DifferentialEquations
#using BenchmarkTools
using Plots
include("DEsystems.jl")



function getAverageErrorByRefs(sol::Vector{Vector{Float64}},solRef::Vector{Any},T::Int,numPoints::Int)
    #numPoints=length(solmliqss.savedTimes[1])
    allErrors=0.0
    for index=1:T
        sumTrueSqr=0.0
        sumDiffSqr=0.0
        relerror=0.0
        for i = 1:numPoints #
            ts=solRef[i][index]
            Ns=sol[i][index]
             if abs(Ns-ts)<1.0 # to avoid big errors due to discontinuities  
                if ts*ts>1e-12
                    sumDiffSqr+=(Ns-ts)*(Ns-ts)
                    sumTrueSqr+=ts*ts
                end
            end
        end
         if  abs(sumTrueSqr)>1e-6
            relerror=sqrt(sumDiffSqr/sumTrueSqr)
        else
            relerror=0.0
        end
        allErrors+= relerror
    end
    return allErrors/T
end

 
 abstol, reltol =  1e-4,1e-3

for alg in [ TRBDF2() ,QNDF2(), Trapezoid()   ,Heun()  ]

@show typeof(alg).name.name
#= tspan=(0.0, 40.0) 
sol=get_Sol_LTI(alg,tspan, abstol, reltol)
BSON.@load joinpath(@__DIR__,  "BSON_data", "ref_sol_LTI_Feagin12.bson") ref_sol_LTI_Feagin12
LTI_err=getAverageErrorByRefs(sol.u,ref_sol_LTI_Feagin12,2,4000)
@show LTI_err =#


tspan=(0.0, 40.0)
sol=get_Sol_Bball(alg,tspan, abstol, reltol)
BSON.@load joinpath(@__DIR__,  "BSON_data", "ref_sol_Bball.bson") ref_sol_Bball
Bball_err=getAverageErrorByRefs(sol.u,ref_sol_Bball,2,4000)
@show Bball_err


tspan = (0.0, 300)
sol=get_Sol_TsodyksMarkram(alg,tspan, abstol, reltol)
BSON.@load joinpath(@__DIR__,  "BSON_data", "ref_sol_TsodyksMarkram.bson") ref_sol_TsodyksMarkram
TsodyksMarkram_err=getAverageErrorByRefs(sol.u,ref_sol_TsodyksMarkram,7,3000)
@show TsodyksMarkram_err

tspan = (0.0, 0.0005)
sol=get_Sol_Interleaved(alg,tspan, abstol, reltol)  
 BSON.@load joinpath(@__DIR__,  "BSON_data", "ref_sol_InterleavedRodas.bson") ref_sol_InterleavedRodas
Interleaved_err=getAverageErrorByRefs(sol.u,ref_sol_InterleavedRodas,5,500)
@show Interleaved_err


tspan= (0.0, 300.0)
sol=get_Sol_IZH(alg,tspan, abstol, reltol)  
BSON.@load joinpath(@__DIR__,  "BSON_data", "ref_sol_IZH.bson") ref_sol_IZH
IZH_err=getAverageErrorByRefs(sol.u,ref_sol_IZH,2,3000)
@show IZH_err

tspan=(0.0,0.005)
sol=get_Sol_Railgun(alg,tspan, abstol, reltol)
@show sol.u[50][1]
BSON.@load joinpath(@__DIR__,  "BSON_data", "ref_sol_Railgun_compact.bson") ref_sol_Railgun_compact
Railgun_err=getAverageErrorByRefs(sol.u,ref_sol_Railgun_compact,18,500)
@show Railgun_err 

end

#= TRBDF2
LTI_err = 0.01935000107940335
Bball_err = 0.06920133785815154 =#