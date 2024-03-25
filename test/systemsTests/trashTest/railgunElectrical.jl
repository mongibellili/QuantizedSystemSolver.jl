using QuantizedSystemSolver
function test()
 
    odeprob = @NLodeProblem begin
          name=(RGElectrical,)
         ROn = 1e-5;ROff = 1e5;
          Lpr = 4.2*1e-9#0.48*1e-6#0.45 * 1e-6
L1 = 0.6*1e-6 #28.0*1e-6#0.6*1e-6
L2 = 4.0e-6;L3 = 1.1*1e-6

 R1= 4.0e-3; R2 = 0.28*1e-3;R3 = 3.6*1e-3

C = 3.08*1e-3#3.08*1e-3

#C = 0.865 *1e-6
γ = 50.0*1e6; w = 15.0*1e-3#25.0*1e-3#15.0*1e-3 
μ = 4.0*3.14*1e-7
Rpr0 = 15.5*1e-6
#FNmec = 680.0; α = 0.154
t0 = 0.001
t0=1.0

Rpr=Rpr0*(sqrt(t0/(t+1e-4))+(t/t0)^16)/(1.0+(t/t0)^16)

          discrete = [1e5,1e5,1e5,1e5,0.427e-6,1e-2,0.0,1.0];u = [10.1,0.0,0.0,0.0,10000.75,10.1,0.0,0.0,0.0,10000.75]
          rd1=discrete[1];rs1=discrete[2];rd2=discrete[3];rs2=discrete[4];LR = discrete[5];nextT=discrete[6];manualIntg=discrete[7]
          is1=u[1] ;dis1=u[2]; il1=u[3] ;dil1=u[4];uc1=u[5];is2=u[6] ;dis2=u[7]; il2=u[8] ;dil2=u[9];uc2=u[10];
          α=LR+Lpr;
          
          RR=manualIntg*2.0*sqrt(μ/(3.14*γ))/w
          dRR=1.0*2.0*sqrt(μ/(3.14*γ))/w
          #dRpr=Rpr
          β=(-rd1-R2-R3-RR-Rpr)*dil1+rd1*dis1-(dRR+Rpr)*(il1+il2)-(RR+Rpr)*dil2;
          θ=(-rd2-R2-R3-RR-Rpr)*dil2+rd2*dis2-(dRR+Rpr)*(il1+il2)-(RR+Rpr)*dil1
          Δ=(L2+L3)*(L2+L3+2.0*α)
        #=   du[1] =dis1
          du[2]=(-(R1+rs1+rd1)*dis1+rd1*dil1+is1/C)/L1
          du[3]=dil1
          du[4]=(θ*rd1*is1-γ*rd2*is2-(θ*α-γ*(rd2+β))*il2-(θ*(rd1+β)-α*γ)*il1)/Δ
          du[5]=-is1/C =#
          du[1] =dis1*discrete[8]
          du[2]=((-(R1+rs1+rd1)*dis2+rd1*dil1+is1/C)/L1)*discrete[8]
          du[3]=dil1*discrete[8]     #  =dil1 --->u[3]=il1
          #du[4]=(θ*rd1*is1-γ*rd2*is2-(θ*α-γ*(rd2+β))*il2-(θ*(rd1+β)-α*γ)*il1)/Δ
          du[4]=(((L2+L3+α)*β-α*θ)/Δ)*discrete[8] #  =ddil1 --->u[4]=dil1
          du[5]=(-is1/C)*discrete[8]
          du[6]=dis2*discrete[8]
          du[7]=((-(R1+rs2+rd2)*dis2+rd2*dil2+is2/C)/L1)*discrete[8]
          du[8]=dil2 *discrete[8]                                     #(θ*rd2*is2-γ*rd1*is1-(θ*α-γ*(rd1+β))*il1-(θ*(rd2+β)-α*γ)*il2)/Δ
          du[9]=(((L2+L3+α)*θ-α*β)/Δ)*discrete[8]  #  =ddil2 --->u[9]=dil1
          du[10]=(-is2/C)*discrete[8]

          if t-0.0000025>0.0 
            rs1=ROn
            rs2=ROn
          end
          if -uc1>0.0 
            rs1=ROff;rd1=ROn
          end     
          if -uc2>0.0 
            rs2=ROff;rd2=ROn
          end   
          if -il1*il2>0   
            discrete[8]=0.0
          end    
                      
       
          if t-0.2*1e-3>0.0
            #LR = (0.416*1e-6-0.397*1e-6)/(0.2*1e-3-1.5*1e-3)*(t-0.2*1e-3)+0.416*1e-6
            LR = (0.427*1e-6-0.397*1e-6)/(0.2*1e-3-1.5*1e-3)*(t[0]-0.2*1e-3)+0.427*1e-6
          end
          if t-1.5*1e-3>0.0
            LR = 0.397*1e-6
          end
          if t-nextT>0
            manualIntg=manualIntg+0.5
            nextT=nextT+1e-2
          end
          



          
    end
    #@show odeprob
    tspan = (0.0, 0.006)
    sol= solve(odeprob,nmliqss2(),tspan,abstol=1e-4,reltol=1e-3)    
    save_Sol(sol,3)
    save_Sol(sol,5)
    save_Sol(sol,8)
    save_Sol(sol,10)
end
test()
