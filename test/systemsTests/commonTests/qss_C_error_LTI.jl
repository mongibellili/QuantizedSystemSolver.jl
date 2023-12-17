
using CSV
using DataFrames
#using Plots


function getAverageErrorByRefs(solT::Vector{Float64},sol::Vector{Float64},solRef::Function,numPoints::Int)
    #numPoints=length(solmliqss.savedTimes[1])
    
        sumTrueSqr=0.0
        sumDiffSqr=0.0
        relerror=0.0
        for i = 1:numPoints #
            ts=solRef(solT[i])
            Ns=sol[i]
            sumDiffSqr+=(Ns-ts)*(Ns-ts)
            sumTrueSqr+=ts*ts
        end
        if  abs(sumTrueSqr)>1e-12
          relerror=sqrt(sumDiffSqr/sumTrueSqr)
          else
            relerror=0.0
          end
        
       
    
    return relerror
  end
function LTI1()
  u1, u2 = -8.73522174738572, -7.385745994549763
  λ1, λ2 = -10.841674966758294, -9.168325033241706
  c1, c2 = 121.14809142478035, -143.14809142478035
  xp1, xp2 = 0.0, 20.0
  x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
  x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
sumError=0.0
x=[x1,x2]
for i=1:2
  df=DataFrame(CSV.File("qss/Csolver/mliqss_LTI1_e-2-5/x$(i).dat"#= , delim=" " =#))
  a,b=names(df)
  t=parse(Float64,a)
  u=parse(Float64,b)
  insert!.(eachcol(df), 1, [t, u])
  rename!(df,[:A,:B])
  err1=getAverageErrorByRefs(df.A,df.B,x[i],length(df.A))
  sumError+=err1
end

error=sumError/2
@show error


end

LTI1()