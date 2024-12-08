# Taylor Variables

## Purpose
QSS methods rely heavily on root finding and derivative computation. The approximation through Taylor variables transforms any complicated equations to polynomials, which makes root finding cheaper. In addition, the use of Taylor variables provide easily the different derivatives for complex expressions. The implementation includes various constructors and methods to manipulate and evaluate the Taylor series.

In the quantized system solver, and similar to the Taylor struct from the package TaylorSeries.jl, the Taylor0 struct serves as a convenient representation of Taylor series approximations while avoiding memory allocation. 

New functions are created to match each old function but with a different name and added caches in the parameters. These new functions are designed to optimize performance by avoiding memory allocation during their execution where arithmetic operations and mathematical functions leverage the existing cached data rather than creating new instances of taylor variables. In addition, the old Taylor is kept, with minimum functionalities, as a fallback in case there is an expression that does not use the available cache vector. 

**Arithmetic operations**: The provided code introduces personalized arithmetic operations such as addsub, subsub, and mulT that are designed to operate directly on existing data structures, utilizing a caching mechanism to store intermediate results.

**Mathematical functions**: Similar to arithmetic operations, new mathematical functions are designed to leverage in-place operations. Functions like exp, log, sin, and cos, sqrt, power iterate over the input Taylor series and apply the corresponding operation without allocating additional arrays.

**Transformation of expressions**: The `transformF` function is designed to translate userdefined mathematical expressions into optimized forms that leverage the custom arithmetic and function implementations. By traversing the expression tree with the prewalk function, it identifies operations such as addition, subtraction, multiplication, division, and specific mathematical functions (like exponential and logarithmic functions). For each identified operation, the code modifies the expression to call specialized versions (e.g., subT, addT, mulT). Additionally, it tracks the number of caches needed for these operations.

## Example:
```julia
function (+)(a::Taylor0, b::Taylor0)
    v = similar(a.coeffs)
    @__dot__ v = (+)(a.coeffs, b.coeffs)
    return Taylor0(v, a.order)
end

function addT(a::Taylor0, b::Taylor0, cache::Taylor0)
  @__dot__ cache.coeffs = (+)(a.coeffs, b.coeffs)
  return cache
end
```
The first function (+) from the TaylorSeries.jl creates a new Taylor variable every time, while the customized arithmetic operation addT uses an existing cache.

## Taylor references
Taylor0 is defined for many functions. However, other functions can also be added.


```@docs
Taylor0
```
### non-allocation operations
```@docs
createT(a::Taylor0, cache::Taylor0)
```

```@docs
createT(a::T, cache::Taylor0) where {T<:Number}
```
```@docs
addT(a::Taylor0, b::Taylor0, cache::Taylor0)
```
```@docs
subT(a::Taylor0, b::Taylor0, cache::Taylor0)
```
```@docs
mulT(a::Taylor0, b::Taylor0, cache1::Taylor0)
```
```@docs
divT(a::Taylor0, b::Taylor0, cache1::Taylor0)
```
```@docs
addsub(a::Taylor0, b::Taylor0, c::Taylor0, cache::Taylor0)
```
```@docs
addsub(a::Taylor0, b::Taylor0, c::T, cache::Taylor0) where {T<:Number}
```
```@docs
addsub(a::T, b::Taylor0, c::Taylor0, cache::Taylor0) where {T<:Number}
```
```@docs
subsub(a::Taylor0, b::Taylor0,c::Taylor0,cache::Taylor0) 
```

```@docs
mulsub(a::P,b::Q,c::R,cache1::Taylor0) where {P,Q,R <:Union{Taylor0,Number}}
```

```@docs
muladdT(a::P,b::Q,c::R,cache1::Taylor0) where {P,Q,R <:Union{Taylor0,Number}}
```


```@docs
negateT(a::Taylor0, cache::Taylor0)
```

```@docs
powerT(a::T, r::S, cache1::Taylor0) where {S<:Real,T<:Number}
```
```@docs
QuantizedSystemSolver.square(a::Taylor0, cache1::Taylor0)
```
 
### non-allocation functions
 
```@docs
sqrt(a::Taylor0, cache1::Taylor0)
```
```@docs
sqrt(a::T, cache1::Taylor0) where {T<:Number}
```
```@docs
exp(a::Taylor0, c::Taylor0)
```
```@docs
exp(a::T, c::Taylor0) where {T<:Number}
```
```@docs
log(a::Taylor0, c::Taylor0)
```
```@docs
sin(a::T, s::Taylor0, c::Taylor0) where {T<:Number}
```
```@docs
cos(a::T, s::Taylor0, c::Taylor0) where {T<:Number}
```
```@docs
tan(a::Taylor0, c::Taylor0, c2::Taylor0)
```
```@docs
asin(a::Taylor0, c::Taylor0, r::Taylor0, cache3::Taylor0)
```
```@docs
acos(a::Taylor0, c::Taylor0, r::Taylor0, cache3::Taylor0)
```
```@docs
atan(a::Taylor0, c::Taylor0, r::Taylor0)
```
```@docs
abs(a::Taylor0, cache1::Taylor0)
```


### Internals

All internals functions of the TaylorSeries.jl are used with the exception of the 
following two functions, which are rewritten to avoid one allocation.
```@docs
QuantizedSystemSolver.asin!
```
```@docs
QuantizedSystemSolver.acos!
```


!!! note "Allocating TaylorSeries functions"

    As a fallback all TaylorSeries functions are used in case an expression
    is not transformed to one of the existing personalized functions above.

## Index

```@index
Pages = ["taylor.md"]
Order = [:type, :function]
```