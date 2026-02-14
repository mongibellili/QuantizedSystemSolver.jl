#= function get_Sol_LTI(alg,tspan,  abstol, reltol)
    function f!(du, u, p, t)
        du[1] = - 10.0*u[1] + u[2]+1.0
        du[2] = - 0.1*u[2]+0.0001*u[1]+2.0
    end
    u0 = [10.0, 10.0]
    p = nothing
    prob = ODEProblem(f!, u0, tspan, p)
    cb = PresetTimeCallback(1:1:tspan[2], function (integrator)
        integrator.u[2] = 10.0
        integrator.u[1] = 10.0
    end)
    sol = solve(prob, alg, callback=cb, abstol=abstol, reltol=reltol, saveat=0.01)
  #  @show sol.stats
  #=   p= plot(sol, marker=:circle, title="LTI_system", markersize=1)
    savefig(p, "$(typeof(alg).name.name)_LTI_system_tol_$(abstol)_$(tspan[2]).png") =#
    
    return sol
end =#

function get_Sol_Bball(alg,tspan,  abstol, reltol)
    function bbal(du, u, p, t)
        du[1] =u[2]
        du[2] =-9.8- 0.01*u[2]    
    end
    function condition(u, t, integrator)
        5.0-u[1]
    end
    function bounce(integrator)
     # integrator.p[1]+=1
        u=integrator.u
        if u[2]<-0.1
            u[2] = - u[2] * 0.98
        else
            integrator.stop!()
        end
    end
    cb=ContinuousCallback(condition, bounce; save_positions=(false, false))
    
    
    u0 = [10.0,0.0]
    p=[0]
    prob = ODEProblem(bbal, u0, tspan,p,callback=cb)
    ts= collect(0.0:0.01:tspan[2])
    sol = solve(prob, alg,reltol=reltol,abstol=abstol,save_everystep=false,dense=false, saveat=ts)
    @show length(sol.t)
   # @show sol.stats
  #=   p= plot(sol, marker=:circle, title="Bball_system", markersize=1)
    savefig(p, "$(typeof(alg).name.name)_Bball_system_tol_$(abstol)_$(tspan[2]).png") =#
   # @btime solve($prob,$alg ,reltol=$reltol,abstol=$abstol,dense=false);
    
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
        gK, gNa, gL, EK, ENa, EL, C, I, tau, tau_u, tau_R, u0, gmax, Esyn = p
        v, n, m, h, u, R, gsyn = u
        du[1] = ((gK * (n^4.0) * (EK - v)) + (gNa * (m^3.0) * h * (ENa - v)) + (gL * (EL - v)) +
                I + gsyn * (Esyn - v)) / C
        du[2] = (alpha_n(v) * (1.0 - n)) - (beta_n(v) * n)
        du[3] = (alpha_m(v) * (1.0 - m)) - (beta_m(v) * m)
        du[4] = (alpha_h(v) * (1.0 - h)) - (beta_h(v) * h)
        # Synaptic variables
        du[5] = -(u / tau_u)
        du[6] = (1 - R) / tau_R
        du[7] = -(gsyn / tau)
    end

    function epsp!(integrator)
        integrator.u[5] += integrator.p[12] * (1 - integrator.u[5])
        integrator.u[7] += integrator.p[13] * integrator.u[5] * integrator.u[6]
        integrator.u[6] -= integrator.u[5] * integrator.u[6]
    end

    epsp_ts = PresetTimeCallback(5:7.2:300, epsp!; save_positions = (false, false))

    n_inf(v) = alpha_n(v) / (alpha_n(v) + beta_n(v))
    m_inf(v) = alpha_m(v) / (alpha_m(v) + beta_m(v))
    h_inf(v) = alpha_h(v) / (alpha_h(v) + beta_h(v))
    p = [35.0, 40.0, 0.3, -77.0, 55.0, -65.0, 1, 0, 30, 1000, 50, 0.5, 0.005, 0]
    u0 = [-60, n_inf(-60), m_inf(-60), h_inf(-60), 0.0, 1.0, 0.0]


    prob = ODEProblem(HH!, u0, tspan, p)
   # sol = solve(prob,alg,reltol=reltol,abstol=abstol ,dense=false; callback=epsp_ts, saveat=0.1);
     ts = collect(0.0:0.1:tspan[2])
    sol = solve(prob,alg ,reltol=reltol,abstol=abstol ,save_everystep=false, saveat=ts, dense=false,callback=epsp_ts);
    @show sol.stats
    @show length(sol.u)
    #display(sol.alg) 
    #display(sol.interp) 
       # display(sol.retcode )              
    #=   p1=plot(sol, idxs = 7,marker=:circle,markersize=1)
        savefig(p1, "Tsodyks-Markram_dense$(typeof(alg).name.name)_ft$(tspan[2])_abst$(abstol).png") =#
         @show sol.stats
      p1= plot(sol,idxs=[1], marker=:circle, markersize=1)
      p2= plot(sol,idxs=[2], marker=:circle, markersize=1)
      p3= plot(sol,idxs=[3], marker=:circle, markersize=1)
        p4= plot(sol,idxs=[4], marker=:circle, markersize=1)
        p56= plot(sol,idxs=[5,6], marker=:circle, markersize=1)
        p7= plot(sol,idxs=[7], marker=:circle, markersize=1)
      p=plot(p1,p2,p3,p4,p56,p7,layout=(2,3))
    savefig(p, "$(typeof(alg).name.name)_Tsodyks-Markram_system_tol_$(abstol)_$(tspan[2]).png")
   # @btime solve($prob,$alg ,reltol=$reltol,abstol=$abstol,dense=false ; callback=$epsp_ts);
 
    return sol  
end  

function get_Sol_Interleaved(alg,tspan,  abstol, reltol)
    k=[0]
    function interleaved(du, u, p, t)
        C=1e4;L=1e4;R=0.1;uu=24.0;T=1e-4;DC=0.5;ROn = 1e-5; ROff = 1e5;N=4.0
        du[1] =C*(u[2]+u[3]-R*u[1]+u[4]+u[5])
        for k in 2:5 
            du[k]=(((uu/p[k+3]) - u[k]) * (p[k+3]*p[k-1]/(p[k+3]+p[k-1])) - u[1])*L;
        end          
    end
    function condition(out, u, t, integrator) 
        T=1e-4;DC=0.5;N=4.0
        out[1] = t-integrator.p[9]
        out[2] =  t-integrator.p[10]-0.01*T
        out[3] = t - integrator.p[10]-T/4.0-0.01*T
        out[4] = t - integrator.p[10]-T/2.0-0.01*T
        out[5] = t - integrator.p[10]-T*3.0/4.0-0.01*T
        out[6] = t - integrator.p[10]-DC*T/N-0.01*T
        out[7] = t - integrator.p[10]-T/4.0-DC*T/N-0.01*T
        out[8] = t - integrator.p[10]-T/2.0-DC*T/N-0.01*T
        out[9] = t - integrator.p[10]-T*3.0/N-DC*T/N-0.01*T
        out[10] = -u[2]
        out[11] = -u[3]
        out[12] = -u[4]
        out[13] = -u[5]
    end

    function affect!(integrator, idx)
        T=1e-4;ROn = 1e-5; ROff = 1e5
        if idx == 1
            integrator.p[10]=integrator.p[9] 
            integrator.p[9]=integrator.p[9]+T   
            k[1] += 1
        elseif idx == 2
            integrator.p[5] = ROn
            integrator.p[1] = ROff
            k[1] += 1
        elseif idx == 3
            integrator.p[6] = ROn
            integrator.p[2] = ROff
            k[1] += 1
        elseif idx == 4
            integrator.p[7] = ROn
            integrator.p[3] = ROff
            k[1] += 1
        elseif idx == 5
            integrator.p[8] = ROn
            integrator.p[4] = ROff
            k[1] += 1
        elseif idx == 6
            integrator.p[5]  = ROff
            integrator.p[1]  = ROn
            k[1] += 1
        elseif idx == 7
            integrator.p[6]  = ROff
            integrator.p[2]  = ROn
            k[1] += 1
        elseif idx == 8
            integrator.p[7]  = ROff
            integrator.p[3] = ROn
            k[1] += 1
        elseif idx == 9
            integrator.p[8] = ROff
            integrator.p[4]  = ROn
            k[1] += 1
        elseif idx == 10
           # println("derivative du2 before reset: ",(((24.0/integrator.p[5]) - integrator.u[2]) * (integrator.p[5]*integrator.p[1]/(integrator.p[5]+integrator.p[1])) - integrator.u[1])*1e4; )
            integrator.p[1] = ROff
            #k[1] += 1
           # println("derivative du2 after reset: ",(((24.0/integrator.p[5]) - integrator.u[2]) * (integrator.p[5]*integrator.p[1]/(integrator.p[5]+integrator.p[1])) - integrator.u[1])*1e4; )
        elseif idx == 11
           # println("derivative du3 before reset: ",(((24.0/integrator.p[6]) - integrator.u[3]) * (integrator.p[6]*integrator.p[2]/(integrator.p[6]+integrator.p[2])) - integrator.u[1])*1e4; )
            integrator.p[2] = ROff
            #k[1] += 1
        elseif idx == 12
           # println("derivative du4 before reset: ",(((24.0/integrator.p[7]) - integrator.u[4]) * (integrator.p[7]*integrator.p[3]/(integrator.p[7]+integrator.p[3])) - integrator.u[1])*1e4; )
            integrator.p[3] = ROff
           # k[1] += 1
        elseif idx == 13
           # println("derivative du5 before reset: ",(((24.0/integrator.p[8]) - integrator.u[5]) * (integrator.p[8]*integrator.p[4]/(integrator.p[8]+integrator.p[4])) - integrator.u[1])*1e4; )
            integrator.p[4] = ROff
           # k[1] += 1
        end
    end


    cbs = VectorContinuousCallback(condition, affect!, 13,save_positions = (false, false) )
      

  
    u0 = [0.0,0.0,0.0,0.0,0.0]
    p = [1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e-4,0.0]
   
    prob = ODEProblem(interleaved, u0, tspan,p, callback=cbs)
   # sol = solve(prob,alg ,reltol=reltol,abstol=abstol,dense=false#= , saveat=1e-6 =#);
    tspan=(0.0,0.0005)
   ts = collect(0.0:1e-6:tspan[2])
    sol = solve(prob,alg ,reltol=reltol,abstol=abstol ,save_everystep=false, saveat=ts, dense=false);
     @show sol.stats
     @show k[1]
     @show length(sol.t)
   #=    p= plot(sol, marker=:circle, title="Interleaved_system", markersize=1)
    savefig(p, "$(typeof(alg).name.name)_Interleaved_densesystem_tol_$(abstol)_$(tspan[2]).png")
    =#

   # @btime solve($prob,$alg ,reltol=$reltol,abstol=$abstol,dense=false#= ; saveat=1e-6 =#);
   
    return sol
end  

function get_Sol_IZH(alg,tspan,  abstol, reltol)
    function izh!(du, u, p, t)
        a, b, c, d, I = p
        du[1] = 0.04 * u[1]^2 + 5.0 * u[1] + 140.0 - u[2] + I
        du[2] = a * (b * u[1] - u[2])

    end
    function thr(u, t, integrator)
        u[1] -30.0
    end
    function reset!(integrator)
        k[1] += 1
       # println("derivative du1 before reset: ", 0.04 * integrator.u[1]^2 + 5.0 * integrator.u[1] + 140.0 - integrator.u[2] + integrator.p[5])
        integrator.u[1] = integrator.p[3]
        integrator.u[2] += integrator.p[4]
        #println("derivative du1 after reset: ", 0.04 * integrator.u[1]^2 + 5.0 * integrator.u[1] + 140.0 - integrator.u[2] + integrator.p[5])
    end
    cb1 =ContinuousCallback(thr, reset!; save_positions = (false, false))
    current_step = PresetTimeCallback(50, integrator -> integrator.p[5] = 10.0; save_positions = (false, false))
    cb = CallbackSet(cb1, current_step)

    p = [0.02, 0.2, -50, 2.0, 0.0]
    u0 = [-60.0, p[2] * -60.0]
    prob = ODEProblem(izh!, u0, tspan, p , callback = cb)
     ts= collect(0.0:0.1:tspan[2])
    sol = solve(prob, alg,reltol=reltol,abstol=abstol,save_everystep=false,dense=false, saveat=ts)
    @show length(sol.t)

    #sol = solve(prob,alg,reltol=reltol,abstol=abstol,saveat=0.1);
    @show sol.stats
   #=  p1=plot(sol, idxs = [1,2],marker=:circle,markersize=1)
    savefig(p1, "IzhikevichchModel$(typeof(alg).name.name)_24nvautosave_.png") =#

    
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
    end
    
    function condition(out, u, t, integrator) # Event when condition(out,u,t,integrator) == 0
            period=0.0004
            out[1] = t-period*1
            out[2] = t-period*2
            out[3] = t-period*3
            out[4] = t-period*4
            out[5] = t- period*5
            out[6] = t- period*6
            out[7] = t- period*7

            out[8] = 1.0-u[1]
            out[9] = 1.0-u[3]
            out[10] = 1.0-u[5]
            out[11] = 1.0-u[7]
            out[12] = 1.0-u[9]
            out[13] = 1.0-u[11]
            out[14] = 1.0-u[13]
            out[15] = 1.0-u[15]

            out[16] = 1.1-u[2]
            out[17] = 1.1-u[4]
            out[18] = 1.1-u[6]
            out[19] = 1.1-u[8]
            out[20] = 1.1-u[10]
            out[21] = 1.1-u[12]
            out[22] = 1.1-u[14]
            out[23] = 1.1-u[16]
    end
    
    function affect!(integrator, idx)
            if idx == 1
                integrator.p[3]=1.0
                integrator.p[4]=1.0
            elseif idx == 2
                integrator.p[5]=1.0
                integrator.p[6]=1.0
            elseif idx == 3
                integrator.p[7]=1.0
                integrator.p[8]=1.0
            elseif idx == 4
                integrator.p[9]=1.0
                integrator.p[10]=1.0
            elseif idx == 5
                integrator.p[11]=1.0
                integrator.p[12]=1.0
            elseif idx == 6
                integrator.p[13]=1.0
                integrator.p[14]=1.0
            elseif idx == 7
                integrator.p[15]=1.0
                integrator.p[16]=1.0
            elseif idx == 8
                integrator.p[2]=0.0
                integrator.u[1]=0.0
            elseif idx == 9
                integrator.p[4]=0.0
                integrator.u[3]=0.0
            elseif idx == 10
                integrator.p[6]=0.0
                integrator.u[5]=0.0
            elseif idx == 11
                integrator.p[8]=0.0
                integrator.u[7]=0.0    
            elseif idx == 12
                integrator.p[10]=0.0
                integrator.u[9]=0.0
            elseif idx == 13
                integrator.p[12]=0.0
                integrator.u[11]=0.0
            elseif idx == 14
                integrator.p[14]=0.0
                integrator.u[13]=0.0
            elseif idx == 15
                integrator.p[16]=0.0
                integrator.u[15]=0.0 

            elseif idx == 16
                integrator.p[1]=0.0
                integrator.u[2]=0.0
            elseif idx == 17
                integrator.p[3]=0.0
                integrator.u[4]=0.0
            elseif idx == 18
                integrator.p[5]=0.0
                integrator.u[6]=0.0
            elseif idx == 19
                integrator.p[7]=0.0
                integrator.u[8]=0.0


            elseif idx == 20
                integrator.p[9]=0.0
                integrator.u[10]=0.0
            elseif idx == 21
                integrator.p[11]=0.0
                integrator.u[12]=0.0    
            elseif idx == 22
                integrator.p[13]=0.0
                integrator.u[14]=0.0
            elseif idx == 23
                integrator.p[15]=0.0
                integrator.u[16]=0.0    

            end
    end
    
        cbs = VectorContinuousCallback(condition, affect!,nothing, 23; save_positions = (false, false) )


        charge1=1.0;charge2=0.0;charge3=0.0;charge4=0.0;charge5=0.0;charge6=0.0;charge7=0.0;charge8=0.0;
        operate1=1.0;operate2=0.0;operate3=0.0;operate4=0.0;operate5=0.0;operate6=0.0;operate7=0.0;operate8=0.0;
        p0 = [operate1, charge1, operate2, charge2, operate3, charge3, operate4, charge4, operate5, charge5, operate6, charge6, operate7, charge7, operate8, charge8];
        u0 = [12.0,0.0,12.0,0.0,12.0,0.0,12.0,0.0,12.0,0.0,12.0,0.0,12.0,0.0,12.0,0.0,0.0,1.0]
        
        prob = ODEProblem(rg,u0,tspan,p0,callback=cbs)
       # sol = solve(prob, alg,reltol=reltol,abstol=abstol,dtmin=1e-12,force_dtmin=true  , saveat=1e-5)
         ts= collect(0.0:1e-5:tspan[2])
        sol = solve(prob, alg,reltol=reltol,abstol=abstol,save_everystep=false,dense=false, saveat=ts)
        @show sol.stats
       #=  p7=plot(sol,idxs=[17],marker=:circle,markersize=1,title="vars 17 (x position)",xlabel="time (s)",ylabel="x (m)")
        p8=plot(sol,idxs=[18],marker=:circle,markersize=1,title="vars 18 (v velocity)",xlabel="time (s)",ylabel="v (m/s)")
        p=plot(p7,p8,layout=(2,1))
        savefig(p,"$(typeof(alg).name.name)simple__15_12_25_period4_rg_8_alg_vars1718_tf$(tspan[2])_tol23.png") =#

      #  @btime solve($prob, $alg,reltol=$reltol,abstol=$abstol ,dtmin=1e-12,force_dtmin=true )
        
        return sol


end 