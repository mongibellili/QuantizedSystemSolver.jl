
using QuantizedSystemSolver
#using XLSX
#using BenchmarkTools
using BSON
#using TimerOutputs
#using Plots
function test(case,solvr)
  absTol=1e-5
     relTol=1e-2
     

      BSON.@load "./ref_bson/solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01
     prob=NLodeProblem(quote
        name=(adrN1000d01,)
        u[1:333]=1.0
        u[334:1000]=0.0
        _dx=100.0#1/dx=N/10=1000/10
        a=1.0;d=0.1;r=1000.0
        #discrete=[0.0]
        du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
        for k in 2:999  
            du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
        end 
        du[1000]=-a*_dx*(u[1000]-u[999])+d*_dx*_dx*(2.0*u[999]-2.0*u[1000])+r*u[1000]*u[1000]*(1.0-u[1000]) 
       #=  if u[1]-10.0>0.0 #fake to test discreteintgrator & loop
          discrete[1]=1.0
        end =#
    end)
    println("start adr solving")
    ttnmliqss=0.0

    tspan=(0.0,5.0)
    solnmliqss=solve(prob,solvr,abstol=absTol,reltol=relTol,tspan)#
   # @show solnmliqss.totalSteps,solnmliqss.simulStepCount  
    #save_Sol(solnmliqss,1,2,600,1000,note="analy1") 
    #@show solnmliqss.savedVars
    #save_Sol(solnmliqss) 
      solnmliqssInterp=solInterpolated(solnmliqss,0.01)
      #@show solnmliqssInterp.savedVars
    err4=getAverageErrorByRefs(solFeagin14VectorN1000d01,solnmliqssInterp)
  #  @show err4,solnmliqss.totalSteps

   #ttnmliqss=@belapsed solve($prob,$solvr,abstol=$absTol,saveat=0.01,reltol=$relTol,$tspan)
    resnmliqss12E_2= ("$(solnmliqss.algName)",relTol,err4,solnmliqss.totalSteps,solnmliqss.simulStepCount,ttnmliqss)
    @show resnmliqss12E_2   



    #=  XLSX.openxlsx("3sys $(solvr)_$(case)_$(relTol).xlsx", mode="w") do xf
        sheet = xf[1]
        sheet["A1"] = "3sys __$case)"
        sheet["A4"] = collect(("solver","Tolerance","error","totalSteps","simul_steps","time"))
        sheet["A5"] = collect(resnmliqss1E_2)
        sheet["A6"] = collect(resnmliqss11E_2)
        sheet["A7"] = collect(resnmliqss12E_2)

       #=  sheet["A8"] = collect(resnmliqss2E_2)
        sheet["A9"] = collect(resnmliqss2E_3)
        sheet["A10"] = collect(resnmliqss2E_4) =#
     end =#
end

case="order1_"
println("golden search")
#test(case,mliqss1())   #goldenSearch
println("compareBounds")
test(case,nmliqss2())  #compareBounds
#test(case,mliqssBounds1())
println("iterations")
#test(case,nliqss1())  #iterations