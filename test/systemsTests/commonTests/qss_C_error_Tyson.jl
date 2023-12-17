
using CSV
using DataFrames
#= using Plots
using BSON =#
using OrdinaryDiffEq

function getAverageErrorByRefs(solT::Vector{Float64},sol::Vector{Float64},solRef,index::Int,numPoints::Int)
  #numPoints=length(solmliqss.savedTimes[1])
 
      sumTrueSqr=0.0
      sumDiffSqr=0.0
      relerror=0.0
      for i = 1:numPoints #
          ts=solRef(solT[i],idxs=index) #diffEq interpolates
          Ns=sol[i]
          sumDiffSqr+=(Ns-ts)*(Ns-ts)
          sumTrueSqr+=ts*ts
      end
      if  abs(sumTrueSqr)>1e-12
        relerror=sqrt(sumDiffSqr/sumTrueSqr)
        else
          relerror=0.0
        end
      
      relerror
  
end
function Tyson()
  function funcName(du,u,p,t)# api requires four args
    du[1] = u[4]-1e6*u[1]+1e3*u[2]
    du[2] =-200.0*u[2]*u[5]+1e6*u[1]-1e3*u[2]
    du[3] = 200.0*u[2]*u[5]-u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)
    du[4] =u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)-u[4]
    du[5] = 0.015-200.0*u[2]*u[5]
    du[6] =u[4]-0.6*u[6]
end
tspan = (0.0,25.0)
u0= [0.0,0.75,0.25,0.0,0.0,0.0]
prob = ODEProblem(funcName,u0,tspan)

 solRodas5P = solve(prob,Rodas5P(),saveat=0.0001,abstol = 1e-12, reltol = 1e-8)

  #BSON.@load "qss/ref_bson/solVect_Tyson_Rodas5Pe-12.bson" solRodas5PVectorTyson
sumError=0.0
filenames=["C2","CP","pM","M","Y","yP"]
for i=1:6
  df=DataFrame(CSV.File("qss/Csolver/mliqss1_TYSON_e-2-5/$(filenames[i]).dat"#= , delim=" " =#))
 #=  a,b=names(df)  # output of c solver in .dat files does not have a title
  t=parse(Float64,a)  # first row seen as title: store values...sometimes converting first row to title change type...no longer can be parsed to float
  u=parse(Float64,b)
  insert!.(eachcol(df), 1, [t, u]) #store them in first row =#
  rename!(df,[:A,:B])    #rename title
  err1=getAverageErrorByRefs(df.A,df.B,solRodas5P,i,length(df.A))
  sumError+=err1
end

error=sumError/6
@show error


end

Tyson()