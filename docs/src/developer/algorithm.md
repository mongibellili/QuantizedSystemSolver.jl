# QSS Algorithms


The Algorithm type is parametric on the solver name (N) and order (O). Currently there are three specific algorithm names: (QSS Algorithm, LiQSS Algorithm, and mLiQSS Algorithm), and two orders (1 and 2).
## Algorithm extension
Any new algorithm can be added via the name N which is of type Val, and with the order O which also of type Val according to this type:

```@docs
 QuantizedSystemSolver.Algorithm{N,O}
```

```@docs
 QSSAlgorithm{N,O}
```

### What is needed with a new algorithm:
  - A new algorithm with the same old name and a different order (>2) requires an addition only of the quantizer methods.
  - A new algorithm with a new name and an old order (1,2) requires an addition of an `integrate` method only.
  - A new algorithm with a new name and a new order requires an addition of an `integrate` method and addition of the quantizer methods.
  - Creating a new superclass similar to the `QSSAlgorithm` requires an addition of a new `solve` function.

### Example
```julia
my_new_algorithm()=QSSAlgorithm(Val(:my_new_algorithm),Val(3))
```
This is a new algorithm name of order 3 belongs to the family of the `QSSAlgorithm`. Only an addition of the `integrate` and the quantizer methods is needed.