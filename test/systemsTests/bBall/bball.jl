

using qss
using BenchmarkTools
#= using Plots;
gr(); =#
#include("D:/models/bball.jl")
function test()
    odeprob = @NLodeProblem begin
        name=(bball,)
        #= parameter1=3000.0# cache can be dynamic....parameters take this feature
        parameter2=0.00001 =#
        u = [20.0,0.0]
        discrete = [0.0]
        du[1] =u[2]
        du[2] =-9.8
        if -u[1]>0   #5*discrte gave error
            u[2]=-u[2]  #discrete=0.0-->type Symbol has no field args...find to personalize error msg  
            #discrete[1]=5.0          
       #=  else
            discrete[1]=-1.0   =#                                  
        end
       #=  if u[2]>0
            discrete[1]=discrete[1]-5.0 #the zc will happen once going up and once going down...so step is -10
        end =#
    end
    tspan=(0.0,20.0)
   sol= solve(odeprob,nmliqss2()#= ,abstol=absTol,saveat=0.01,reltol=relTol =#,tspan)
   #save_Sol(sol,xlims=(0.0,15.0) ,ylims=(-2.04e-1,2.06e-1))
   save_Sol(sol)
   save_SolDer(sol)
end
#@time 
test()

