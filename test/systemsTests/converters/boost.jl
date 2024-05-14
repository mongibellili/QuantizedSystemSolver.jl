using QuantizedSystemSolver
using BenchmarkTools


function test()
    odeprob = @NLodeProblem begin
        name=(boost,)
        C = 1e-4; L = 1e-4; R = 10.0;U = 24.0; T = 1e-4; DC = 0.5; ROn = 1e-5;ROff = 1e5;
        discrete = [1e5,1e-5,1e-4,0.0,0.0]
        u = [0.0,0.0]
        rd=discrete[1];rs=discrete[2];nextT=discrete[3];lastT=discrete[4];diodeon=discrete[5]
        il=u[1] ;uc=u[2]

        id=(il*rs-uc)/(rd+rs) # diode's current
        ir=uc/R               # ohm's law at the load
        is=il-id              # kirchoff's 1st law left node
        ic=id-ir              # kirchoff's 1st law right node
        ul=U-rs*(is)          # kirchoff's 2nd law in left mesh
        du[1] =ul/L           # voltage-current formula in Inductance
        du[2]=ic/C            # voltage-current formula in Capacitor

        if t-nextT>0.0 # time to close the switch
            lastT=nextT
            nextT=nextT+T
            rs=ROn
        end

        if t-lastT-DC*T>0.0 #time to open the switch
            rs=ROff
        end                          
        
      if diodeon*(id)+(1.0-diodeon)*(id*rd)-0.6>0
        rd=ROn
        diodeon=1.0
      end 

      if -(diodeon*(id)+(1.0-diodeon)*id)>0
        rd=ROff
        diodeon=0.0
      end
    #=   if diodeon*(id)+(1.0-diodeon)*(id*rd)>0
        rd=ROn
        diodeon=1.0
      else
        rd=ROff
        diodeon=0.0
      end  =#
           
    end
    tspan=(0.0,0.0025)
   sol= solve(odeprob,nmliqss2(),abstol=1e-4,reltol=1e-3,tspan)
            
  save_Sol(sol)
 # save_Sol(sol,xlims=(0.0,0.0006) ,ylims=(-2.04e-1,40.0))
end
#@btime 
test()
