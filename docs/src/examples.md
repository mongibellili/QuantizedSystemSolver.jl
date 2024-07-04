
```julia
odeprob = NLodeProblem(quote
    name=(sysb1,)
    u = [-1.0, -2.0]
    du[1] = -2.0
    du[2] =1.24*u[1]-0.01*u[2]+0.2
end)  
tspan=(0.0,1.0)
```

```julia
sol=solve(odeprob,qss1(),tspan)
```

```julia
sol=solve(odeprob,qss2(),tspan)
```

```julia
sol=solve(odeprob,liqss1(),tspan)
```

```julia
sol=solve(odeprob,liqss2(),tspan)
```

```julia
sol=solve(odeprob,nmliqss1(),tspan)
```

```julia
sol=solve(odeprob,nmliqss2(),tspan)
```