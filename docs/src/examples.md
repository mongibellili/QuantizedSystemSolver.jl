
# Examples:

### Systems of 2 Linear Time Invariant Differential equations
```julia
odeprob = NLodeProblem(quote
     name=(sysb53,)
    u = [-1.0, -2.0]
    du[1] = -20.0*u[1]-80.0*u[2]+1600.0
    du[2] =1.24*u[1]-0.01*u[2]+0.2
end)  
tspan=(0.0,1.0)
```

```julia
sol=solve(odeprob,qss1(),tspan)
save_Sol(sol)
```
![plot_sysb53_qss1](./assets/img/plot_sysb53_qss1.png)
```julia
sol=solve(odeprob,qss2(),tspan)
save_Sol(sol)
```
![plot_sysb53_qss2](./assets/img/plot_sysb53_qss2.png)
```julia
sol=solve(odeprob,liqss1(),tspan)
save_Sol(sol)
```
![plot_sysb53_liqss1](./assets/img/plot_sysb53_liqss1.png)
```julia
sol=solve(odeprob,liqss2(),tspan)
save_Sol(sol)
```
![plot_sysb53_liqss2](./assets/img/plot_sysb53_liqss2.png)
```julia
sol=solve(odeprob,nmliqss1(),tspan)
save_Sol(sol)
```
![plot_sysb53_nmliqss1](./assets/img/plot_sysb53_nmliqss1.png)
```julia
sol=solve(odeprob,nmliqss2(),tspan)
save_Sol(sol)
```
![plot_sysb53_nmliqss2](./assets/img/plot_sysb53_nmliqss2.png)


### The Tyson Model
```julia
function test(solvr,absTol,relTol)
odeprob = NLodeProblem(quote
    name=(tyson,)
    u = [0.0,0.75,0.25,0.0,0.0,0.0]
    du[1] = u[4]-1e6*u[1]+1e3*u[2]
    du[2] =-200.0*u[2]*u[5]+1e6*u[1]-1e3*u[2]
    du[3] = 200.0*u[2]*u[5]-u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)
    du[4] =u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)-u[4]
    du[5] = 0.015-200.0*u[2]*u[5]
    du[6] =u[4]-0.6*u[6]
end ) 
println("start tyson solving")
tspan=(0.0,25.0)
sol=solve(odeprob,solvr,abstol=absTol,reltol=relTol,tspan)
println("start saving plot")
save_Sol(sol)
end

absTol=1e-5
relTol=1e-2
solvrs=[qss1(),liqss1(),nmliqss1(),nmliqss2()]
for solvr in solvrs
    test(solvr,absTol,relTol)
end
```

![plot_tyson_qss1](./assets/img/plot_tyson_qss1.png)


![plot_tyson_liqss1](./assets/img/plot_tyson_liqss1.png)


![plot_tyson_nmliqss1](./assets/img/plot_tyson_nmliqss1.png)


![plot_tyson_nmliqss2](./assets/img/plot_tyson_nmliqss2.png)


###    Oregonator; Vanderpl; Loktavoltera

![plot_oregonator_mliqss1](./assets/img/plot_oregonator_mliqss1.png)

![plot_vanderpol_mliqss2](./assets/img/plot_vanderpol_mliqss2.png)

![plot_loktavoltera_mliqss3](./assets/img/plot_loktavoltera_mliqss3.png)





