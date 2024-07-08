# Taylor0 
These are just some examples. Taylor0 is defined for many other functions. However, other functions can also be added.

```@docs
Taylor0
```

```@docs
createT(a::T,cache::Taylor0) where {T<:Number}
```
```@docs
addT(a::Taylor0, b::Taylor0,cache::Taylor0) 
```
```@docs
subT(a::Taylor0, b::Taylor0,cache::Taylor0)
```
```@docs
mulT(a::Taylor0, b::Taylor0,cache1::Taylor0)
```
```@docs
divT(a::Taylor0, b::Taylor0,cache1::Taylor0) 
```
```@docs
addsub(a::Taylor0, b::Taylor0,c::Taylor0,cache::Taylor0) 
```
```@docs
addsub(a::T, b::Taylor0,c::Taylor0,cache::Taylor0) where {T<:Number} 
```
```@docs
negateT(a::Taylor0,cache::Taylor0)
```
```@docs
subsub(a::Taylor0, b::Taylor0,c::Taylor0,cache::Taylor0) 
```
```@docs
muladdT(a::P,b::Q,c::R,cache1::Taylor0) where {P,Q,R <:Union{Taylor0,Number}}
```
```@docs
mulsub(a::P,b::Q,c::R,cache1::Taylor0) where {P,Q,R <:Union{Taylor0,Number}}
```
