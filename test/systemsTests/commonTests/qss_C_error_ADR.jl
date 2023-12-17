
using CSV
using DataFrames

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
function ADR()
  function funcName(du,u,p,t)# api requires four args
    _dx=100.0#1/dx=N/10=1000/10
    a=1.0
    d=0.1
    r=1000.0
    du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
    for k in 2:999  
        du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
    end 
    du[1000]=-a*_dx*(u[1000]-u[999])+d*_dx*_dx*(2.0*u[999]-2.0*u[1000])+r*u[1000]*u[1000]*(1.0-u[1000]) 
  
  end
  tspan = (0.0,10.0)
  
  u0=zeros(1000)
  u0[1:333].=1.0


  prob = ODEProblem(funcName,u0,tspan)

   solFeagin14 = solve(prob,Feagin14(),saveat=0.0001,abstol = 1e-12, reltol = 1e-8)


sumError=0.0 

for i=1:1000

  df=DataFrame(CSV.File("qss/Csolver/advectionN1000d01r1000mliqss2e-2/u[$(i)].dat"#= , delim=" " =#))  #order2:error = 0.02960787395245143
 # df=DataFrame(CSV.File("qss/Csolver/order1_mliqss_adr1000d01_e-2-5/mliqss_adr/u[$(i)].dat"#= , delim=" " =#))#order1:error = 0.14204469389384278
  #=a,b=names(df) 
   a,b=names(df)  # output of c solver in .dat files does not have a title
  t=parse(Float64,a)  # first row seen as title: store values...sometimes converting first row to title change type...no longer can be parsed to float
  u=parse(Float64,b) =#
  nr=nrow(df)
  if nr>0  # 1 step tables do not add to error
      #insert!.(eachcol(df), 1, [t, u]) #store them in first row
        rename!(df,[:A,:B])    #rename title
      err1=getAverageErrorByRefs(df.A,df.B,solFeagin14,i,length(df.A))
      sumError+=err1
      #= else
        push!(df,[t, u]) =#
  end
  
  
  

end

error=sumError/1000
@show error


end

ADR()