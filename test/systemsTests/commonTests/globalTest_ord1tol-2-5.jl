
using QuantizedSystemSolver
#using XLSX
#using BenchmarkTools
using BSON
#using TimerOutputs
#using Plots
function test(case,solvr)
  absTol=1e-5
     relTol=1e-2
      odeprob = @NLodeProblem begin
         #sys b53
         name=(sysb53,)
         u = [-1.0, -2.0]
         du[1] = -20.0*u[1]-80.0*u[2]+1600.0
         du[2] =1.24*u[1]-0.01*u[2]+0.2
     end  
    # @show odeprob.prname
     u1, u2 = -8.73522174738572, -7.385745994549763
     λ1, λ2 = -10.841674966758294, -9.168325033241706
     c1, c2 = 121.14809142478035, -143.14809142478035
     xp1, xp2 = 0.0, 20.0
     x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
     x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2

     timenmliqss=0.0;er1=0.0;er2=0.0
     
     println("start LTI solving")
     tspan=(0.0,100.0)
     solnmliqss=solve(odeprob,solvr,abstol=absTol,saveat=0.01,reltol=relTol,tspan#= ,maxErr=10000*relTol =#)
     @show solnmliqss.totalSteps
   # save_Sol(solnmliqss,note="xi10dxi")

     solnmliqssInterp=solInterpolated(solnmliqss,0.01)
     er1=getError(solnmliqssInterp,1,x1)  
     er2=getError(solnmliqssInterp,2,x2) 
   # timenmliqss=@belapsed solve($odeprob,$solvr,abstol=$absTol,saveat=0.01,reltol=$relTol,tspan#= ,maxErr=1000*$relTol =#)
     resnmliqss1E_2= ("$(solnmliqss.algName)",relTol,(er1+er2)/2,solnmliqss.totalSteps,solnmliqss.simulStepCount,timenmliqss)
     @show resnmliqss1E_2

    
    
  BSON.@load "ref_bson/solVect_Tyson_Rodas5Pe-12.bson" solRodas5PVectorTyson
     odeprob = @NLodeProblem begin
         name=(tyson,)
         u = [0.0,0.75,0.25,0.0,0.0,0.0]
         du[1] = u[4]-1e6*u[1]+1e3*u[2]
         du[2] =-200.0*u[2]*u[5]+1e6*u[1]-1e3*u[2]
         du[3] = 200.0*u[2]*u[5]-u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)
         du[4] =u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)-u[4]
         du[5] = 0.015-200.0*u[2]*u[5]
         du[6] =u[4]-0.6*u[6]
     end   
     println("start tyson solving")
     timenmliqss=0.0
     tspan=(0.0,25.0)
     solnmliqss=solve(odeprob,solvr,abstol=absTol,saveat=0.01,reltol=relTol,tspan#= ,maxErr=100*relTol =#)
     @show solnmliqss.totalSteps
    #save_Sol(solnmliqss,note="cancelcriteria-3_xi10dxi"#= xlims=(10.3778695,14.5789) =#)  


    # save_Sol(solnmliqss,1,note="x1 intrval13  ",xlims=(4.0,4.38),ylims=(0.0007,0.000723))
     solnmliqssInterp=solInterpolated(solnmliqss,0.01)
    err3=getAverageErrorByRefs(solRodas5PVectorTyson,solnmliqssInterp) 

  # timenmliqss=@belapsed solve($odeprob,$solvr,abstol=$absTol,saveat=0.01,reltol=$relTol,tspan#= ,maxErr=1000*$relTol =#)
    resnmliqss11E_2= ("$(solnmliqss.algName)",relTol,err3,solnmliqss.totalSteps,solnmliqss.simulStepCount,timenmliqss)
    @show resnmliqss11E_2 
     


      BSON.@load "ref_bson/solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01
     prob=@NLodeProblem begin
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
    end
    println("start adr solving")
    ttnmliqss=0.0

    tspan=(0.0,10.0)
    solnmliqss=solve(prob,solvr,abstol=absTol,saveat=0.01,reltol=relTol,tspan)#
    @show solnmliqss.totalSteps,solnmliqss.simulStepCount  
    #save_Sol(solnmliqss,1,2,600,1000,note="analy1") 
    #@show solnmliqss.savedVars
    #save_Sol(solnmliqss) 
      solnmliqssInterp=solInterpolated(solnmliqss,0.01)
      #@show solnmliqssInterp.savedVars
    err4=getAverageErrorByRefs(solFeagin14VectorN1000d01,solnmliqssInterp)
    @show err4,solnmliqss.totalSteps

  # ttnmliqss=@belapsed solve($prob,$solvr,abstol=$absTol,saveat=0.01,reltol=$relTol,tspan)
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
test(case,nmliqss1())
#test(case,nmliqss2())