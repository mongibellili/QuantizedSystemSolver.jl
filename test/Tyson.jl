
using QuantizedSystemSolver
#using XLSX
using BenchmarkTools
#using BSON
#using TimerOutputs
#using Plots
function test(case,solvr)
  absTol=1e-5
     relTol=1e-2
   

    
     
 # BSON.@load "formalA2/ref_bson/solVect_Tyson_Rodas5Pe-12.bson" solRodas5PVectorTyson
     odeprob = NLodeProblem(quote
         name=(tyson,)
         u = [0.0,0.75,0.25,0.0,0.0,0.0]
         du[1] = u[4]-1e6*u[1]+1e3*u[2]
         du[2] =-200.0*u[2]*u[5]+1e6*u[1]-1e3*u[2]
         du[3] = 200.0*u[2]*u[5]-u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)
         du[4] =u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)-u[4]
         du[5] = 0.015-200.0*u[2]*u[5]
         du[6] =u[4]-0.6*u[6]
     end  ) 
     println("start tyson solving")
     timenmliqss=0.0
     tspan=(0.0,25.0)
     solnmliqss=solve(odeprob,solvr,abstol=absTol,saveat=0.01,reltol=relTol,tspan,maxiters=10000#= ,maxErr=100*relTol =#)
     save_Sol(solnmliqss)
    # @show solnmliqss.totalSteps
    #save_Sol(solnmliqss,note="cancelcriteria-3_xi10dxi"#= xlims=(10.3778695,14.5789) =#)  


    # save_Sol(solnmliqss,1,note="x1 intrval13  ",xlims=(4.0,4.38),ylims=(0.0007,0.000723))
    # solnmliqssInterp=solInterpolated(solnmliqss,0.01)
     err3=0.0
  #  err3=getAverageErrorByRefs(solRodas5PVectorTyson,solnmliqssInterp) 

  #timenmliqss=@belapsed solve($odeprob,$solvr,abstol=$absTol,saveat=0.01,reltol=$relTol,$tspan#= ,maxErr=1000*$relTol =#)
    resnmliqss11E_2= ("$(solnmliqss.algName)",relTol,err3,solnmliqss.stats.totalSteps,solnmliqss.stats.simulStepCount,timenmliqss)
    @show resnmliqss11E_2 

end
case="order1_"

println("compareBounds")
test(case,nmliqss1())  #compareBounds
