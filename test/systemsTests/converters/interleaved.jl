using qss
using BenchmarkTools
#= using Plots;
gr(); =#

function test()
    odeprob = @NLodeProblem begin
        name=(interleaved,)
        #= parameter1=3000.0# cache can be dynamic....parameters take this feature
        parameter2=0.00001 =#
        C=1e4;L=1e4;R=0.1;uu=24.0;T=1e-4;DC=0.5;ROn = 1e-5; ROff = 1e5;N=4.0
        u = [0.0,0.0,0.0,0.0,0.0]
        discrete = [1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e-8,0.0]
        du[1] =C*(u[2]+u[3]-R*u[1]+u[4]+u[5])
        #= for k in 2:5 
            du[k]=(((uu/discrete[k+3]) - u[k]) * (discrete[k+3]*discrete[k-1]/(discrete[k+3]+discrete[k-1])) - u[1])*L;
        end  =#
        
            du[2]=(((uu/discrete[5]) - u[2]) * (discrete[5]*discrete[1]/(discrete[5]+discrete[1])) - u[1])*L
            du[3]=(((uu/discrete[6]) - u[3]) * (discrete[6]*discrete[2]/(discrete[6]+discrete[2])) - u[1])*L
            du[4]=(((uu/discrete[7]) - u[4]) * (discrete[7]*discrete[3]/(discrete[7]+discrete[3])) - u[1])*L
            du[5]=(((uu/discrete[8]) - u[5]) * (discrete[8]*discrete[4]/(discrete[8]+discrete[4])) - u[1])*L
        

            if t-discrete[9]>0.0 
                discrete[10]=discrete[9] 
                discrete[9]=discrete[9]+T                           
            end

            if t-discrete[10]-0.01*T>0.0
                discrete[5] = ROn
                discrete[1] = ROff
            end
            #zcf3
            if t-discrete[10]-T/4.0-0.01*T>0.0
                discrete[6] = ROn
                discrete[2] = ROff
            end
            if t-discrete[10]-T/2.0-0.01*T>0.0
                discrete[7] = ROn
                discrete[3] = ROff
            end
            if t-discrete[10]-T*3.0/4.0-0.01*T>0.0
                discrete[8] = ROn
                discrete[4] = ROff
            end
            #zcf6
            if t - discrete[10]-DC*T/N-0.01*T>0 
                discrete[5]  = ROff
                discrete[1]  = ROn
            end 

            if t -  discrete[10]-T/4.0-DC*T/N-0.01*T>0 
                discrete[6]  = ROff
                discrete[2]  = ROn
            end 
            if t -  discrete[10]-T/2.0-DC*T/N-0.01*T>0 
                discrete[7]  = ROff
                discrete[3] = ROn
            end 
            if t -  discrete[10]-T*3.0/N-DC*T/N-0.01*T>0 
                discrete[8] = ROff
                discrete[4]  = ROn
            end 

            #zcf10
            if -u[2]>0 
                discrete[1] = ROff
            end 
            if -u[3]>0 
                discrete[2] = ROff
            end 

            if -u[4]>0 
                discrete[3] = ROff
            end 
            if -u[5]>0 
                discrete[4] = ROff
            end 
    end
    tspan=(0.0,0.0005)
   sol= solve(odeprob,nmliqss2(),abstol=1e-4,reltol=1e-3,tspan)
  # @show sol
   save_Sol(sol)
  # save_Sol(sol,xlims=(0.0,15.0) ,ylims=(-2.04e-1,2.06e-1))
end
#@time 
test()

