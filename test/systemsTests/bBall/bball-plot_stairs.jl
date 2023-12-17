

using qss
using BenchmarkTools
#= using Plots;
gr(); =#
#include("D:/models/bball.jl")
function test()
    odeprob = @NLodeProblem begin
        name=(bballStairs,)
        #= parameter1=3000.0# cache can be dynamic....parameters take this feature
        parameter2=0.00001 =#
        u = [20.0,-0.001]
        discrete = [15.0]
        du[1] =u[2]
        du[2] =-9.8
        if discrete[1]-u[1]>0   #5*discrte gave error
            u[2]=-u[2]*0.5   #discrete=0.0-->type Symbol has no field args...find to personalize error msg  
           # discrete[1]=discrete[1]-5.0          
       #=  else
            discrete[1]=-1.0   =#                                  
        end

        if u[2]>0
            discrete[1]=discrete[1]-5.0 #the zc will happen once going up and once going down...so step is -10
        end
                                        #=    if u[2]>0
                                                u[1]=5.0

                                            end =#
    end
    #sol= solve(odeprob,2.3,qss2(),saveat(0.01),0.0,1e-6,1e-3)
   # save_prob_to_model(odeprob,"D:/models/bball.jl","bball") #any location you want

   tspan=(0.0,10.0)
   sol= solve(odeprob,qss2()#= ,abstol=absTol,saveat=0.01,reltol=relTol =#,tspan)
    save_Sol(sol)


   #=  sol=QSS_Solve_from_model(bball,odeprob,10.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(sol) =#
   #=  save_Sol(sol," ";xlims=(1.40,1.433),ylims=(-2.04e-1,2.06e-1)) =#
    #= sol= solve(odeprob,1.0,qss3(),saveat(0.01),0.0,1e-6,1e-3)
      save_Sol(sol) =#
end
#@time 
test()
#= function test2()
    odeprob = @NLodeProblem begin
        parameter2=30.0# why 300 bad
        parameter1=1e5
        u = [10.0,0.0]
        discrete = [0.0]
        du[1] =u[2]
        du[2] =-9.8-(discrete[1])*(parameter1*u[1]+parameter2*u[2])
        if -u[1]>0   #5*discrte gave error
            discrete[1]=1.0   #discrete=0.0-->type Symbol has no field args...find to personalize error msg            
        else
            discrete[1]=0.0                                    
        end
    end
    sol= solve(odeprob,1.5,qss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(sol)
    
   #=  sol= solve(odeprob,1.0,qss3(),saveat(0.01),0.0,1e-6,1e-3)
      save_Sol(sol) =#
end
test2() =#
