#= function get_Sol_LTI(alg,tspan,  abstol, reltol)
    function f!(du, u, p, t)
        du[1] = - 10.0*u[1] + u[2]+1.0
        du[2] = - 0.1*u[2]+0.0001*u[1]+2.0
         if t>p[1]
            p[1]+=1.0
            u[2] = 10.0
            u[1] = 10.0
            
        end
    end
    u0 = [10.0, 10.0]
    p = [1.0]
    prob = ODEProblem(f!, u0, tspan, p)
   
    sol = solve(prob, alg, abstol=abstol, reltol=reltol)
   # @show sol.stats
    #=     p= plot(sol, marker=:circle, title="LTI_system", markersize=1)
    savefig(p, "$(alg)_LTI_system_tol_$(abstol)_$(tspan[2]).png") =#
    
    return sol
end
=#
function get_Sol_Bball(alg,tspan,  abstol, reltol)
    function bbal(du, u, p, t)
        du[1] =u[2]
        du[2] =-9.8- 0.01*u[2]
        if u[1]<0.0
            if  u[2]<-0.1  
                u[2]=-u[2]*0.98
            else
                t=Inf  # stop the simulation
            end
        end  
    end
   
    
    
    u0 = [10.0,0.0]
    p=[0]
    prob = ODEProblem(bbal, u0, tspan,p)
    sol = solve(prob, alg,reltol=reltol,abstol=abstol)
    @show sol.stats
        p= plot(sol, marker=:circle, markersize=1,)
    savefig(p, "$(alg)_Bball_system3_tol_$(abstol)_$(tspan[2]).png") 
    
    return sol
end 
 

function get_Sol_TsodyksMarkram(alg,tspan,  abstol, reltol)
    
        alpha_n(v) = (0.02 * (v - 25.0)) / (1.0 - exp((-1.0 * (v - 25.0)) / 9.0))
        beta_n(v) = (-0.002 * (v - 25.0)) / (1.0 - exp((v - 25.0) / 9.0))

        # Sodium ion-channel rate functions
        alpha_m(v) = (0.182 * (v + 35.0)) / (1.0 - exp((-1.0 * (v + 35.0)) / 9.0))
        beta_m(v) = (-0.124 * (v + 35.0)) / (1.0 - exp((v + 35.0) / 9.0))

        alpha_h(v) = 0.25 * exp((-1.0 * (v + 90.0)) / 12.0)
        beta_h(v) = (0.25 * exp((v + 62.0) / 6.0)) / exp((v + 90.0) / 12.0)
        function HH!(du, u, p, t)
            gK, gNa, gL, EK, ENa, EL, C, I, tau, tau_u, tau_R, u0, gmax, Esyn,timer = p
            v, n, m, h, uu, R, gsyn = u

            alpha_n = (0.02 * (v - 25.0)) / (1.0 - exp((-1.0 * (v - 25.0)) / 9.0))
            beta_n = (-0.002 * (v - 25.0)) / (1.0 - exp((v - 25.0) / 9.0))
            
            # Sodium ion-channel rate functions
            alpha_m = (0.182 * (v + 35.0)) / (1.0 - exp((-1.0 * (v + 35.0)) / 9.0))
            beta_m = (-0.124 * (v + 35.0)) / (1.0 - exp((v + 35.0) / 9.0))
            
            alpha_h = 0.25 * exp((-1.0 * (v + 90.0)) / 12.0)
            beta_h = (0.25 * exp((v + 62.0) / 6.0)) / exp((v + 90.0) / 12.0)

            du[1] = ((gK * (n^4.0) * (EK - v)) + (gNa * (m^3.0) * h * (ENa - v)) + (gL * (EL - v)) +
                    I + gsyn * (Esyn - v)) / C

            du[2] = (alpha_n * (1.0 - n)) - (beta_n * n)
            du[3] = (alpha_m * (1.0 - m)) - (beta_m * m)
            du[4] = (alpha_h * (1.0 - h)) - (beta_h * h)

            # Synaptic variables
            du[5] = -(uu / tau_u)
            du[6] = (1 - R) / tau_R
            du[7] = -(gsyn / tau)

            if t>=timer
                timer+=7.2
                u[5] += p[12] * (1 - u[5])
                u[7] += p[13] * u[5] * u[6]
                u[6] -= u[5] * u[6]
            end

        end



        n_inf(v) = alpha_n(v) / (alpha_n(v) + beta_n(v))
        m_inf(v) = alpha_m(v) / (alpha_m(v) + beta_m(v))
        h_inf(v) = alpha_h(v) / (alpha_h(v) + beta_h(v))
        p = [35.0, 40.0, 0.3, -77.0, 55.0, -65.0, 1, 0, 30, 1000, 50, 0.5, 0.005, 0,5.0]
        u0 = [-60, n_inf(-60), m_inf(-60), h_inf(-60), 0.0, 1.0, 0.0]
        prob = ODEProblem(HH!, u0, tspan, p)
        sol = solve(prob,alg,reltol=reltol,abstol=abstol );
    @show sol.stats
      p1= plot(sol,idxs=[1], marker=:circle, markersize=1)
      p2= plot(sol,idxs=[2], marker=:circle, markersize=1)
      p3= plot(sol,idxs=[3], marker=:circle, markersize=1)
        p4= plot(sol,idxs=[4], marker=:circle, markersize=1)
        p56= plot(sol,idxs=[5,6], marker=:circle, markersize=1)
        p7= plot(sol,idxs=[7], marker=:circle, markersize=1)
      p=plot(p1,p2,p3,p4,p56,p7,layout=(2,3))
    savefig(p, "$(alg)_Tsodyks-Markram_system_tol_$(abstol)_$(tspan[2]).png")
   # @btime solve($prob,$alg ,reltol=$reltol,abstol=$abstol);
    #@time solve(prob,alg ,reltol=reltol,abstol=abstol);
    return sol
end  

function get_Sol_Interleaved(alg,tspan,  abstol, reltol)
    function interleaved(du, u, p, t)
        C=1e4;L=1e4;R=0.1;T=1e-4;DC=0.5;ROn = 1e-5; ROff = 1e5;N=4.0
        du[1] =C*(u[2]+u[3]-R*u[1]+u[4]+u[5])
        for k in 2:5 
            du[k]=(((24.0/p[k+3]) - u[k]) * (p[k+3]*p[k-1]/(p[k+3]+p[k-1])) - u[1])*L;
        end 
            if t-p[9]>0.0 
                p[10]=p[9] 
                p[9]=p[9]+T                           
            end
            if t-p[10]-0.01*T>0.0
                p[5] = ROn
                p[1] = ROff
            end
            if t-p[10]-T/4.0-0.01*T>0.0
                p[6] = ROn
                p[2] = ROff
            end
            if t-p[10]-T/2.0-0.01*T>0.0
                p[7] = ROn
                p[3] = ROff
            end
            if t-p[10]-T*3.0/4.0-0.01*T>0.0
                p[8] = ROn
                p[4] = ROff
            end
            if t - p[10]-DC*T/N-0.01*T>0 
                p[5]  = ROff
                p[1]  = ROn
            end 
            if t -  p[10]-T/4.0-DC*T/N-0.01*T>0 
                p[6]  = ROff
                p[2]  = ROn
            end 
            if t -  p[10]-T/2.0-DC*T/N-0.01*T>0 
                p[7]  = ROff
                p[3] = ROn
            end 
            if t -  p[10]-T*3.0/N-DC*T/N-0.01*T>0 
                p[8] = ROff
                p[4]  = ROn
            end 
            if -u[2]>0 
                p[1] = ROff
            end 
            if -u[3]>0 
                p[2] = ROff
            end 
            if -u[4]>0 
                p[3] = ROff
            end 
            if -u[5]>0 
                p[4] = ROff
            end 
    end
    u0 = [0.0,0.0,0.0,0.0,0.0]
    p = [1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e-4,0.0]
    prob = ODEProblem(interleaved, u0, tspan, p)    
    sol = solve(prob,alg ,reltol=reltol,abstol=abstol);
     @show sol.stats
     # p= plot(sol, marker=:circle, title="Interleaved_system", markersize=1)
    #savefig(p, "$(alg)_Interleaved_system_tol_$(abstol)_$(tspan[2]).png")
   # @btime solve($prob,$alg,reltol=$reltol,abstol=$abstol);

    return sol
end   


function get_Sol_IZH(alg,tspan,  abstol, reltol)
   
    function izh!(du, u, p, t)
        a, b, c, d, I = p
        du[1] = 0.04 * u[1]^2 + 5.0 * u[1] + 140.0 - u[2] + I
        du[2] = a * (b * u[1] - u[2])
        if u[1] >= 30.0
            u[1] = c
            u[2] += d
        end
        if t >= 50.0
            p[5] += 10.0
        end
    end
    p = [0.02, 0.2, -50.0, 2.0, 0.0]
    u0 = [-60.0, p[2] * -60.0]

    prob = ODEProblem(izh!, u0, tspan, p)
    sol = solve(prob,alg,reltol=reltol,abstol=abstol);
     @show sol.stats
    #=   p= plot(sol, marker=:circle, title="IzhikevichchModel_system", markersize=1)
    savefig(p, "$(alg)_IzhikevichchModel_system_tol_$(abstol)_$(tspan[2]).png") =#
   # @btime solve($prob,$alg,reltol=$reltol,abstol=$abstol);
    return sol
end


function get_Sol_Railgun(alg,tspan,  abstol, reltol)
    
    function rg(du,u,p,t)
        L1 = 1e4 
        R1= 1e-3; 
        LR2= 2e3; 
        C = 100.0;
        _m=8.0
        FNmec = 2.0;
        period=0.0004
          
          operate1=p[1];charge1=p[2]; 
          operate2=p[3];charge2=p[4]; 
          operate3=p[5];charge3=p[6]; 
          operate4=p[7];charge4=p[8]; 
          operate5=p[9];charge5=p[10]; 
          operate6=p[11];charge6=p[12]; 
          operate7=p[13];charge7=p[14]; 
          operate8=p[15];charge8=p[16]; 

          uc1=u[1]; il1=u[2]    ;uc2=u[3]; il2=u[4] ;    uc3=u[5]; il3=u[6]   ;uc4=u[7];  il4=u[8] ;
          uc5=u[9]; il5=u[10]  ;uc6=u[11]; il6=u[12] ;  uc7=u[13]; il7=u[14]  ;uc8=u[15]; il8=u[16] ;   
          x=u[17]; v=u[18] 
          Il=il1+il2+il3+il4-1.0+il5+il6+il7+il8+1.0
          F=1e-2*Il*Il-FNmec
          du[1]=(-il1*C)*charge1
          du[2]= operate1*(charge1*L1*(uc1-R1*il1)+(1.0-charge1)*(-LR2*il1))
      
          du[3]=(-il2*C)*charge2
          du[4]=operate2*(charge2*L1*(uc2-R1*il2)+(1.0-charge2)*(-LR2*il2))

          du[5]=(-il3*C)*charge3
          du[6]= operate3*(charge3*L1*(uc3-R1*il3)+(1.0-charge3)*(-LR2*il3))
      
          du[7]=(-il4*C)*charge4
          du[8]=operate4*(charge4*L1*(uc4-R1*il4)+(1.0-charge4)*(-LR2*il4))

          du[9]=(-il5*C)*charge5
          du[10]=operate5*(charge5*L1*(uc5-R1*il5)+(1.0-charge5)*(-LR2*il5))

          du[11]=(-il6*C)*charge6
          du[12]=operate6*(charge6*L1*(uc6-R1*il6)+(1.0-charge6)*(-LR2*il6))

          du[13]=(-il7*C)*charge7
          du[14]=operate7*(charge7*L1*(uc7-R1*il7)+(1.0-charge7)*(-LR2*il7))

          du[15]=(-il8*C)*charge8
          du[16]=operate8*(charge8*L1*(uc8-R1*il8)+(1.0-charge8)*(-LR2*il8))
         
          du[17]=v
          du[18]=F*_m
          #events
          if t-period>0.0 
            operate2=1.0
            charge2=1.0
          end 
          if t-period*2>0.0 
            operate3=1.0
            charge3=1.0
          end 
          if t-period*3>0.0 
            operate4=1.0
            charge4=1.0
          end 
          if t-period*4>0.0 
            operate5=1.0
            charge5=1.0
          end 
          if t-period*5>0.0 
            operate6=1.0
            charge6=1.0
          end
          if t-period*6>0.0 
            operate7=1.0
            charge7=1.0
          end 
          if t-period*7>0.0 
            operate8=1.0
            charge8=1.0
          end
        
          if uc1<1.0
            charge1=0.0 
            uc1=0.0
          end 
          if uc2<1.0
            charge2=0.0
            uc2=0.0
          end 

          if uc3<1.0
            charge3=0.0
            uc3=0.0
          end
          if uc4<1.0
            charge4=0.0
            uc4=0.0
          end
          if uc5<1.0
            charge5=0.0
            uc5=0.0
          end
          if uc6<1.0
            charge6=0.0
            uc6=0.0
          end
          if uc7<1.0
            charge7=0.0
            uc7=0.0
          end
          if uc8<1.0
            charge8=0.0
            uc8=0.0
          end

          if il1<1.1
            il1=0.0
            operate1=0.0
          end
           if il2<1.1
            il2=0.0
            operate2=0.0
          end
          if il3<1.1
            il3=0.0
            operate3=0.0
          end
          if il4<1.1
            il4=0.0
            operate4=0.0
          end
          if il5<1.1
            il5=0.0
            operate5=0.0
          end
           if il6<1.1
            il6=0.0
            operate6=0.0
          end
          if il7<1.1
            il7=0.0
            operate7=0.0
          end
          if il8<1.1
            il8=0.0
            operate8=0.0
          end
    end

    charge1=1.0;charge2=0.0;charge3=0.0;charge4=0.0;charge5=0.0;charge6=0.0;charge7=0.0;charge8=0.0;
    operate1=1.0;operate2=0.0;operate3=0.0;operate4=0.0;operate5=0.0;operate6=0.0;operate7=0.0;operate8=0.0;
    p0 = [operate1, charge1, operate2, charge2, operate3, charge3, operate4, charge4, operate5, charge5, operate6, charge6, operate7, charge7, operate8, charge8];
    u0 = [12.0,0.0,12.0,0.0,12.0,0.0,12.0,0.0,12.0,0.0,12.0,0.0,12.0,0.0,12.0,0.0,0.0,0.2]
    
    prob = ODEProblem(rg,u0,tspan,p0)
    sol = solve(prob, alg,reltol=reltol,abstol=abstol)
    @show sol.stats
   #=  p7=plot(sol,idxs=[17],marker=:circle,markersize=1,title="vars 17 (x position)",xlabel="time (s)",ylabel="x (m)")
      p8=plot(sol,idxs=[18],marker=:circle,markersize=1,title="vars 18 (v velocity)",xlabel="time (s)",ylabel="v (m/s)")
      p=plot(p7,p8,layout=(2,1))
      savefig(p,"$(alg)simple_period4_rg_8_alg_vars1718_tf$(tspan[2])_tol23.png") =#
       p_i=plot(sol,idxs=[2,4,6,8,10,12,14,16],marker=:circle,markersize=1,xlabel="time (s)",ylabel="currents (A)")
        savefig(p_i,"$(alg)simple_period4_rg_8_alg_currents_tf$(tspan[2])_tol23.png")

    #  @btime solve($prob, $alg,reltol=$reltol,abstol=$abstol )
    
    return sol
end 


function get_Sol_Interleaved(alg,tspan,  abstol, reltol)
    function interleaved(du, u, p, t)
        C=1e4;L=1e4;R=0.1;T=1e-4;DC=0.5;ROn = 1e-5; ROff = 1e5;N=4.0
        du[1] =C*(u[2]+u[3]-R*u[1]+u[4]+u[5])
        #= for k in 2:5 
            du[k]=(((24.0/p[k+3]) - u[k]) * (p[k+3]*p[k-1]/(p[k+3]+p[k-1])) - u[1])*L;
        end  =#
            du[2]=(((24.0/p[5]) - u[2]) * (p[5]*p[1]/(p[5]+p[1])) - u[1])*L;
            du[3]=(((24.0/p[6]) - u[3]) * (p[6]*p[2]/(p[6]+p[2])) - u[1])*L;
            du[4]=(((24.0/p[7]) - u[4]) * (p[7]*p[3]/(p[7]+p[3])) - u[1])*L;
            du[5]=(((24.0/p[8]) - u[5]) * (p[8]*p[4]/(p[8]+p[4])) - u[1])*L;
            if t-p[9]>0.0 
                p[10]=p[9] 
                p[9]=p[9]+T                           
            end
            if t-p[10]-0.01e-4>0.0
                p[5] = ROn
                p[1] = ROff
            end
            if t-p[10]-0.26e-4>0.0
                p[6] = ROn
                p[2] = ROff
            end
            if t-p[10]-0.51e-4>0.0
                p[7] = ROn
                p[3] = ROff
            end
            if t-p[10]-0.76e-4>0.0
                p[8] = ROn
                p[4] = ROff
            end
            if t - p[10]-0.135e-4>0 
                p[5]  = ROff
                p[1]  = ROn
            end 
            if t - p[10]-0.385e-4>0 
                p[6]  = ROff
                p[2]  = ROn
            end 
            if t - p[10]-0.635e-4>0 
                p[7]  = ROff
                p[3] = ROn
            end 
            if t - p[10]-0.885e-4>0 
                p[8] = ROff
                p[4]  = ROn
            end 
            if -u[2]>0 
                p[1] = ROff
            end 
            if -u[3]>0 
                p[2] = ROff
            end 
            if -u[4]>0 
                p[3] = ROff
            end 
            if -u[5]>0 
                p[4] = ROff
            end 
    end
    u0 = [0.0,0.0,0.0,0.0,0.0]
    p = [1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e-4,0.0]
    prob = ODEProblem(interleaved, u0, tspan, p)    
     sol = solve(prob,alg ,reltol=reltol,abstol=abstol);
     @show sol.stats
     p= plot(sol, marker=:circle, title="Interleaved_system", markersize=1)
    savefig(p, "$(alg)_Interleaved_system_tol_$(abstol)_$(tspan[2]).png")
    @btime solve($prob,$alg,reltol=$reltol,abstol=$abstol);
   
    return sol
end  
 