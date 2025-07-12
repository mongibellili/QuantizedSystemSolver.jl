# IR


## Goal: separate the problem parsing from the problem construction.

problem parsing is delegated to a new Module named **SimpleModelIR**:
```Julia
module SimpleModelIR
using MacroTools: postwalk 
include("utils.jl")
include("types.jl")
include("build_ir.jl")
include("normalize_ir.jl")
export
    AbstractODEStatement,
    AssignStatement,
    IfStatement,
    ForStatement,
    WhileStatement,
    ExprStatement,
    ODEFunctionIR,
    build_ir,
    normalize_ir,
    problem_to_normalized_ir,
end
```
The _build_ir_ function loops through the user code expression and builds an ODEFunctionIR(statement). i.e the expression will be changed to a vector of statements. Each statement is an object of type ....<: AbstractODEStatement
The _normalize_ir_ function loops over the vectors and call personalized functions such as _changeVarNames_params_ or _normalize_if_statement!_ to change composite if-statements to simple if-statements.

## Example:
```Julia
  function simpleModel(dy,y,p,t)# api requires four args
      U = 24.0; 
      rd,rs=p;
      il=y[1] 
      id=(il*rs-U)/(rd+rs)
      dy[1] =id
      dy[2]=il
      if t>0 || y[2]>0 && y[1]>y[2]
        y[1]=0.0
      else
        y[2]=0.0
      end       
   end
```
build_ir_ -->

```Julia
ir = SimpleModelIR.ODEFunctionIR(SimpleModelIR.AbstractODEStatement[SimpleModelIR.AssignStatement(:U, 24.0), 
SimpleModelIR.AssignStatement(:((rd, rs)), :p), SimpleModelIR.AssignStatement(:il, :(y[1])), SimpleModelIR.AssignStatement(:id, :
((il * rs - U) / (rd + rs))), SimpleModelIR.AssignStatement(:(dy[1]), :id), SimpleModelIR.AssignStatement(:(dy[2]), :il), 
SimpleModelIR.IfStatement(:(t > 0 || y[2] > 0 && y[1] > y[2]), :(if t > 0 || y[2] > 0 && y[1] > y[2]
      y[1] = 0.0
  else
      y[2] = 0.0
  end))])
```
normalize_ir_ -->

```Julia
ir = SimpleModelIR.ODEFunctionIR(SimpleModelIR.AbstractODEStatement[SimpleModelIR.AssignStatement(:U, 24.0), 
SimpleModelIR.AssignStatement(:((rd, rs)), :p), SimpleModelIR.AssignStatement(:il, :(q[1])), SimpleModelIR.AssignStatement(:id, :
((q[1] * p[2] - 24.0) / (p[1] + p[2]))), SimpleModelIR.AssignStatement(:(dy[1]), :((q[1] * p[2] - 24.0) / (p[1] + p[2]))), 
SimpleModelIR.AssignStatement(:(dy[2]), :(q[1])), SimpleModelIR.IfStatement(:(t - 0.0), :(if t - 0.0
      if t > 0.0 || q[2] > 0.0 && q[1] > q[2]
          q[1] = 0.0
      else
          q[2] = 0.0
      end
  else
      if t > 0.0 || q[2] > 0.0 && q[1] > q[2]
          q[1] = 0.0
      else
          q[2] = 0.0
      end
  end)), SimpleModelIR.IfStatement(:(q[2] - 0.0), :(if q[2] - 0.0
      if t > 0.0 || q[2] > 0.0 && q[1] > q[2]
          q[1] = 0.0
      else
          q[2] = 0.0
      end
  else
      if t > 0.0 || q[2] > 0.0 && q[1] > q[2]
          q[1] = 0.0
      else
          q[2] = 0.0
      end
  end)), SimpleModelIR.IfStatement(:(q[1] - q[2]), :(if q[1] - q[2]
      if t > 0.0 || q[2] > 0.0 && q[1] > q[2]
          q[1] = 0.0
      else
          q[2] = 0.0
      end
  else
      if t > 0.0 || q[2] > 0.0 && q[1] > q[2]
          q[1] = 0.0
      else
          q[2] = 0.0
      end
  end))])
```






```@docs
QuantizedSystemSolver.problem_to_normalized_ir(expr::Expr, stateVarName::Symbol, discrParamName::Symbol)
```

```@docs
QuantizedSystemSolver.build_ir(expr::Expr)
```



```@docs
QuantizedSystemSolver.AbstractODEStatement
```


```@docs
QuantizedSystemSolver.AssignStatement
```


```@docs
QuantizedSystemSolver.IfStatement
```



```@docs
QuantizedSystemSolver.ForStatement
```



```@docs
QuantizedSystemSolver.ExprStatement
```





