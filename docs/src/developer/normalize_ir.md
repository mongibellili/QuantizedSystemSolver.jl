# Normalize_IR


## allow composite if_statments and normalize

allow event conditions to be composite (contain &&, || ). 

### Example

```Julia
if u[1]+u[2] > u[3] && u < 10 || u > 5 && u < 15 || u > 20 
    a = 0
end
```

### Solution


case1)
```Julia
if A>0 && B>0
   X
end
```
<-->
```Julia
if A>0 
    if  B>0
      X
   end
end   
if B>0 
    if  A>0
      X
   end
end  
```
case2)
```Julia
if A>0 && B>0
   X
else
  Y
end
```
<-->
```Julia
if A>0 
    if  B>0
      X
   end
else
   Y
end   
if B>0 
    if  A>0
      X
   end
else
   Y
end  
```
case3)
```Julia
if A>0 || B>0
   X
end
```
<-->
```Julia
if A>0 
      X
end   
if B>0 
      X
end 
```
case4)
```Julia
if A>0 || B>0
   X
else
  Y
end
```
<-->
```Julia
if A>0 
   X
else
      if  !(B>0)
         Y
      end
end   
if B>0 
      X
else
     if !(A>0)
       Y
    end
end  
```
case5) multi
```Julia
if A>0 || B>0 && C>0 #.....
      #(with or without else : the logic is the same)
end
```
<-->
```Julia
if A>0 
      #whole user if-statment and its body
end
if B>0 
      #whole user if-statment and its body
end
if C>0 
     #whole user if-statment and its body
end
```



```@docs
QuantizedSystemSolver.changeVarNames_params(ex::Expr,stateVarName::Symbol,discrParamName::Symbol,muteVar::Symbol,param::Dict{Symbol,Union{Float64,Int64,Expr,Symbol}},helperFunSymSet::Set{Symbol})
```


```@docs
QuantizedSystemSolver.changeVarNames_params(element::Symbol,stateVarName::Symbol,discrParamName::Symbol,muteVar::Symbol,param::Dict{Symbol,Union{Float64,Int64,Expr,Symbol}},helperFunSymSet::Set{Symbol})
```

```@docs
QuantizedSystemSolver.recurse(e::Expr,flattened::Vector{Expr})
```

```@docs
QuantizedSystemSolver.decompose_condition(cond::Expr)
```


```@docs
QuantizedSystemSolver.to_zcf(expr::Expr)
```


```@docs
QuantizedSystemSolver.process_if_condition(cond, stateVarName, discrParamName, param, helperFunSymSet)
```

```@docs
QuantizedSystemSolver.process_if_block(block_expr::Expr, stateVarName, discrParamName, param, helperFunSymSet)
```

```@docs
QuantizedSystemSolver.process_if_expr(statement,stateVarName,discrParamName,param,helperFunSymSet)
```

```@docs
QuantizedSystemSolver.normalize_ir(ir, stateVarName::Symbol, discrParamName::Symbol)
```

