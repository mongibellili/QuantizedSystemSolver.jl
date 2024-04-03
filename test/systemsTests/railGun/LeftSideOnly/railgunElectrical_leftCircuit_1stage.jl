using QuantizedSystemSolver
function test()
 
    odeprob = @NLodeProblem begin
          name=(RGElectrical,)
         ROn = 1e-5;ROff = 1e5;
         # Lpr = 4.2*1e-9#0.48*1e-6#0.45 * 1e-6
L1 = 0.6*1e-4 #28.0*1e-6#0.6*1e-6
L2 = 4.0e-6;L3 = 1.1*1e-6

 R1= 4.0e-3; R2 = 0.28*1e-3;R3 = 3.6*1e-3

C = 3.08*1e-3#3.08*1e-3

#= #C = 0.865 *1e-6
γ = 50.0*1e6; w = 15.0*1e-3#25.0*1e-3#15.0*1e-3 
μ = 4.0*3.14*1e-7
Rpr0 = 15.5*1e-6
#FNmec = 680.0; α = 0.154

t0=1.0 =#

#Rpr=Rpr0*(sqrt(t0/(t+1e-4))+(t/t0)^16)/(1.0+(t/t0)^16)

          discrete = [1e5,1e-5,0.0];u = [0.0,10000.75,0.0]
          rd1=discrete[1];rs1=discrete[2];diodeon=discrete[3]#LR = discrete[3];nextT=discrete[4];manualIntg=discrete[5]
          is1=u[1] ;uc1=u[2]; il1=u[3] 
          #α=LR+Lpr;
          RR=4.0
          #RR=manualIntg*2.0*sqrt(μ/(3.14*γ))/w
         # dRR=1.0*2.0*sqrt(μ/(3.14*γ))/w
          #dRpr=Rpr
          #β=(-rd1-R2-R3-RR-Rpr)*dil1+rd1*dis1-(dRR+Rpr)*(il1);
         
         # Δ=(L2+L3)*(L2+L3+2.0*α)
        #=   du[1] =dis1
          du[2]=(-(R1+rs1+rd1)*dis1+rd1*dil1+is1/C)/L1
          du[3]=dil1
          du[4]=(θ*rd1*is1-γ*rd2*is2-(θ*α-γ*(rd2+β))*il2-(θ*(rd1+β)-α*γ)*il1)/Δ
          du[5]=-is1/C =#
          du[1] =(-(R1+rs1+rd1)*is1+rd1*il1+uc1)/L1
          du[2]=-is1/C
          du[3]=((-rd1-R2-R3-RR)*il1+rd1*is1)/(L2+L3)#((L2+L3+α)*β)/Δ
          

          #= if t-0.00025>0.0 
            rs1=ROn
           
          end =#
          if -uc1>0.0 
            rs1=ROff;rd1=ROn
          end 
        #=   if -is1>0.0 
            rs1=ROff;rd1=ROn
          end   =#   
                               
          if diodeon*((il1-is1))+(1.0-diodeon)*((il1-is1)*rd1)>0
            rd1=ROn;diodeon=1.0
          else
            rd1=ROff;diodeon=0.0
          end 
         #=  if t-0.2*1e-3>0.0
            #LR = (0.416*1e-6-0.397*1e-6)/(0.2*1e-3-1.5*1e-3)*(t-0.2*1e-3)+0.416*1e-6
            LR = (0.427*1e-6-0.397*1e-6)/(0.2*1e-3-1.5*1e-3)*(t[0]-0.2*1e-3)+0.427*1e-6
          end
          if t-1.5*1e-3>0.0
            LR = 0.397*1e-6
          end
          if t-nextT>0
            manualIntg=manualIntg+0.5
            nextT=nextT+1e-2
          end =#
          



          
    end
    #@show odeprob
    tspan = (0.0, 0.6)
    sol= solve(odeprob,nmliqss2(),tspan,abstol=1e-4,reltol=1e-3)    
    save_Sol(sol,1)
    save_Sol(sol,2)
    save_Sol(sol,3)
    
   
end
test()
