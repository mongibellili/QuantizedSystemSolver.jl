
using QuantizedSystemSolver
#using XLSX
#using BenchmarkTools
#using BSON

 function Oregon(solRodas5PVectorOregon,note,solvr,cycleDetection)
  
    
       odeprob = @NLodeProblem begin
            name=(linear,)
            u  = [1.0,1.0,10.0,10.0]
            # u  = [-1.0,-2.0]
            #= du[1] = u[2]-u[1]*u[2]-u[1]
            du[2] = u[2]+u[1]*u[2]-u[3]-2.0
            du[3] = u[2]*u[1] -u[3]*u[1]-0.5 =#
            du[1] = -80.0*u[2]-20.0*u[1]+20.0*u[4]
            du[2] = -0.01*u[2]+1.24*u[1]-20.0*u[3]
            du[3] =  -0.05*u[3]
            du[4] =  -0.05*u[4]
       end  
     
    
       println("start solving")
       timenmliqss=0.0
     
       absTol=1e-5
       relTol=1e-3
       tspan=(0.0,200.0)
       sol=solve(odeprob,solvr,abstol=absTol,saveat=0.01,reltol=relTol,tspan,maxiters=100000#= ,maxErr=100*relTol =#)
      save_Sol(sol)
      #solInterp=solInterpolated(sol,0.01)
      err3=0.0
      #err3=getAverageErrorByRefs(solRodas5PVectorOregon,solInterp)
     # timenmliqss=@belapsed QSS_Solve($odeprob,nmliqss1(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0,maxErr=10*$relTol)
      resnmliqss1E_2= ("$(sol.algName)",err3,sol.stats.totalSteps,sol.stats.simulStepCount,timenmliqss)
      @show resnmliqss1E_2

      return resnmliqss1E_2
    
end 
solRodas5PVectorOregon=0.0
note="abs-3_1Delta" 
resmliqss21=Oregon(solRodas5PVectorOregon,note,nmliqss2(),5)
#resmliqss22=Oregon(solRodas5PVectorOregon,note,mliqss2_SimulIter(),2)
#= resmliqss23=Oregon(solRodas5PVectorOregon,note,mliqss2_SimulIter(),3)
resmliqss24=Oregon(solRodas5PVectorOregon,note,mliqss2_SimulIter(),4) =#
