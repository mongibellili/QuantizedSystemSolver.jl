# Algorithm Extension 

Currently only QSS1,2,3 ; LiQSS1,2,3 ; and mLiQSS1,2,3 exist. Any new algorithm can be added via the name N which is of type Val, with O is the order which also of type Val. For example qss1 is created by:
```julia
qss1()=QSSAlgorithm(Val(:qss),Val(1))
```

```@docs
 QuantizedSystemSolver.ALGORITHM{N,O}
```