using QuantizedSystemSolver
using Test
t1=Taylor0([1.0,1.0,0.0],2)
t2=Taylor0([0.5,2.0,3.0],2)
t3=Taylor0([2.0,1.0,0.0],2)
cache1=Taylor0([1.0,1.0,1.0],2)
cache2=Taylor0([1.0,1.0,1.0],2)
testTaylor(t1)
@test  exp(t2)[0]≈1.6487212707001282
@test  log(t2)[0]≈-0.6931471805599453
@test  sin(t2)[0]≈0.479425538604203
@test  cos(t2)[0]≈0.8775825618903728
@test  tan(t2)[0]≈0.5463024898437905
@test  asin(t2)[0]≈0.5235987755982989
@test  acos(t2)[0]≈1.0471975511965979
@test  atan(t2)[0]≈0.4636476090008061

@test  sinh(t2)[0]≈0.5210953054937474
@test  cosh(t2)[0]≈1.1276259652063807
@test  tanh(t2)[0]≈0.46211715726000974
@test  asinh(t2)[0]≈0.48121182505960347
@test  acosh(t3)[0]≈1.3169578969248166
@test  atanh(t2)[0]≈0.5493061443340549
@test  abs(t2)[0]≈0.5
@test  abs(-t2)[0]≈0.5

cache3=Taylor0([0.0,0.0,0.0],2)
@test   exp(t2,cache1)[0]≈1.6487212707001282
@test   log(t2,cache1)[0]≈-0.6931471805599453
@test   sin(t2,cache1,cache2)[0]≈0.479425538604203
@test   cos(t2,cache1,cache2)[0]≈0.8775825618903728
@test   tan(t2,cache1,cache2)[0]≈0.5463024898437905
@test   asin(t2,cache1,cache2,cache3)[0]≈0.5235987755982989
@test   acos(t2,cache1,cache2,cache3)[0]≈1.0471975511965979 
@test   atan(t2,cache1,cache2)[0]≈0.4636476090008061
@test   abs(t2,cache1)[0]≈0.5

@test   (t2^0)[0]≈1.0
@test   (t2^1)[0]≈0.5
@test   (t2^2)[0]≈0.25
@test   (t2^3)[0]≈0.125
@test   (t2^0.0)[0]≈1.0
@test   (t2^1.0)[0]≈0.5
@test   (t2^2.0)[0]≈0.25
@test   (t2^3.0)[0]≈0.125

@test   (t2^0.5)[0]≈0.7071067811865476
@test   sqrt(t2)[0]≈0.7071067811865476
@test   powerT(t2,0,cache1)[0]≈1.0
@test   powerT(t2,1,cache1)[0]≈0.5
@test   powerT(t2,2,cache1)[0]≈0.25
@test   powerT(t2,3.0,cache1)[0]≈0.125
@test   sqrt(t2,cache1)[0]≈0.7071067811865476
t2=1.0
@test    exp(t2,cache1)[0]≈2.718281828459045
@test    log(t2,cache1)[0]≈0.0
@test    sin(t2,cache1,cache2)[0]≈0.8414709848078965
@test    cos(t2,cache1,cache2)[0]≈0.5403023058681398
@test    tan(t2,cache1,cache2)[0]≈1.5574077246549023
@test    abs(t2,cache1)[0]≈1.0
t4=Taylor0([-2.0,1.0,0.0],2)
@test    abs(t4,cache1)[0]≈2.0


@test    (t2^2)≈1.0
@test    (t2^3.0)≈1.0
@test    sqrt(t2)≈1.0
@test    powerT(t2,2,cache1)[0]≈1.0
@test    powerT(t2,3.0,cache1)[0]≈1.0
@test    sqrt(t2,cache1)[0]≈1.0