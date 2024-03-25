

using QuantizedSystemSolver
using BenchmarkTools
#= using Plots;
gr(); =#
#include("D:/models/bball.jl")
function test()
    odeprob = @NLodeProblem begin
        name=(railgun1,)
        m=20.0 #mass 20g
        I0=1000.0 #max current 1000A
        w=60.0#freq of electric source
        τ=0.05#time constant
        l=1.0# 1meter bullet distance to exit
        u = [1.0,0.0,0.0,0.0]
        y=u[1];vy=u[2];x=u[3];vx=u[4];
        discrete = [1.0,0.0,0.0,0.1,1.0,0.0]
        
        ka1=100.0 #magnetic force coeff in [0,t1] #another study to find these
        ka2=100.0 #magnetic force coeff in [t1,t2]
        ka3=100.0 #magnetic force coeff in [t3,...]
        currSwitch1=discrete[1] #magnetic force coeff in [0,t1] #initially this is switched on
        currSwitch2=discrete[2] #magnetic force coeff in [t1,t2]
        currSwitch3=discrete[3] #magnetic force coeff in [t2,...]
        f=discrete[4] #friction coefficient initially=0.1
        moveBullet=discrete[5] #initially=1
       F=currSwitch1*ka1*I0*I0*sin(w*t)+currSwitch2*ka2*I0*I0+currSwitch3*ka3*I0*I0*exp(-2.0*t/τ)
       exit=discrete[6]
      # F2=currSwitch1*ka1*I0*I0*w*cos(w*t)+currSwitch3*ka3*I0*I0*(-2.0/τ)*exp(-2.0*t/τ)
       # F=currSwitch1*ka1*I0*I0*sin(w*t)+currSwitch2*ka2*I0*I0+currSwitch3*ka3*I0*I0
       # F=exp(2*t)
        du[1] =moveBullet*vy              #u1==y
        du[2] =(-9.8*(exit)-0.001)*moveBullet            #u2==vy # how to compute dvy to account for 
        du[3]=moveBullet*vx               #u3==x
        du[4]=moveBullet*(F-f*vx)/m       #u4=vx
       
        if t-1.0e-5>0   #10ns switch current to mode 2 # constant
           
            currSwitch1=0.0 
            currSwitch2=1.0   
          
           #=  @show currSwitch1,currSwitch2,currSwitch3
            @show moveBullet*(F-f*vx)/m    =#                         
        end
        if t-5.0e-5>0   #50ns switch current to mode 3 #drop
         
            currSwitch2=0.0 
            currSwitch3=1.0  
           
           #=  @show currSwitch1,currSwitch2,currSwitch3
            @show moveBullet*(F-f*vx)/m     =#                          
        end
        if x-l>0 #bullet exited gun
            currSwitch3=0.0 #no more magnetic force
            f=10.0 # friction to account only for air resistance drag... later f=coef*area/2*v*v
            exit=1.0
          #  @show currSwitch1,currSwitch2
          #  @show moveBullet*(F-f*vx)/m  
          # @show u,t
        end
        if -y>0 #bullet touched the ground
            moveBullet=0.0
        end
        if -vx>0 #vx=0
            moveBullet=0.0
        end

      
    end
    tspan=(0.0,2.5)
   # @show odeprob
   
  sol= solve(odeprob,nmliqss2()#= ,abstol=absTol,saveat=0.01,reltol=relTol =#,tspan)
 # save_Sol(sol,3,note="saveindex:l1_x current in Mode1:rise",xlims=(0.0,1.1e-5),ylims=(0.0,0.00025)) #look before [0,t1]
# save_Sol(sol,3,note="saveindex:l1_x afterCurrentSwithedToMode2:constant",xlims=(1.1e-5,5.5e-5),ylims=(0.00019,0.01)) #look at [t1,t2]
 # save_Sol(sol,3,note="saveindex:l1_x afterCurrentSwithedToMode3:drop",xlims=(5.0e-5,0.00099),ylims=(0.004,1.0)) #look at t2 
#  save_Sol(sol,3,note="saveindex:l1_x afterBulletExited:noForce",xlims=(0.0006,0.001),ylims=(0.95,5.0)) #look after bullet exited
  save_Sol(sol,3,note="saveindex:l1_x followTrajectory",xlims=(0.0,2.5)#= ,ylims=(0.0,5000.0) =#) #followTrajectory
  save_Sol(sol,1,note="saveindex:l1_y followTrajectory",xlims=(0.0,2.5)#= ,ylims=(0.0,5000.0) =#) #followTrajectory
# save_SolDer(sol,4,note="saveall:l1_dvx current in Mode1:rise",xlims=(0.0,0.9e-5)#= ,ylims=(0.0,0.025) =#) #look after t1
#  save_SolDer(sol,4,note="saveall:l1_dvx afterCurrentSwithedToMode2:constant",xlims=(1.1e-5,5.9e-5)#= ,ylims=(0.025,350.0) =#) #look after t1
 # save_SolDer(sol,4,note="saveall:l1_dvx afterCurrentSwithedToMode3:drop",xlims=(5.5e-5,0.0006)#= ,ylims=(0.0,0.03) =#) #look after t2 
 # save_Sol(sol,4,note="saveindex:l1_vx afterBulletExited:noForce",xlims=(0.0006,0.001)#= ,ylims=(0.0,2.0) =#) #look after bullet exited
  save_Sol(sol,4,note="saveindex:l1_vx followTrajectory",xlims=(0.0,2.5)#= ,ylims=(0.0,5000.0) =#) #followTrajectory
  save_Sol(sol,2,note="saveindex:l1_vy followTrajectory",xlims=(0.0,2.5)#= ,ylims=(0.0,5000.0) =#) #followTrajectory


#=   save_SolDer(sol,3,note="saveindex:l1_dx afterCurrentSwithedToMode2:constant",xlims=(0.0,1.5e-5),ylims=(0.0,12.0)) #look after t1
 save_SolDer(sol,3,note="saveindex:l1_dx afterCurrentSwithedToMode3:drop",xlims=(0.0,5.5e-5),ylims=(0.0,2000.0)) #look after t2 
 save_SolDer(sol,3,note="saveindex:l1_dx afterBulletExited:noForce",xlims=(0.0,8.0e-4),ylims=(0.0,2.0e4)) #look after bullet exited


   save_SolDer(sol,4,note="saveindex:l1_dvx afterCurrentSwithedToMode2:constant",xlims=(0.0,1.5e-5),ylims=(0.0,1.0e5)) #look after t1
  save_SolDer(sol,4,note="saveindex:l1_dvx afterCurrentSwithedToMode3:drop",xlims=(0.0,5.5e-5),ylims=(1.0e7,1.0e8)) #look after t2 
  save_SolDer(sol,4,note="saveindex:l1_dvx afterBulletExited:noForce",xlims=(0.0,8.0e-4),ylims=(0.0,1.0e8)) #look after bullet exited =#






end
#@time 
test()

