

using Test
using QuantizedSystemSolver
using LinearAlgebra

t1=Taylor0([1.0,1.0,0.0],2)
t2=Taylor0([0.5,2.0,3.0],2)
t3=Taylor0([2.0,1.0,0.0],2)
t4=Taylor0([0.0,0.1,0.5],2)
vectTaylor1=[t1,t2,t3]
@test norm(vectTaylor1) ≈ 2.29128784747792
matrixTaylor=[t1 t2;t3 t4]
@test transpose(matrixTaylor) == [1.0 2.0; 0.5 0.0]
@test adjoint(matrixTaylor)  == [1.0 2.0; 0.5 0.0]
vectTaylor2=[t2,t1,t3-t4]
@test dot(vectTaylor1,vectTaylor2) ≈ 5.0
@test dot(vectTaylor2,vectTaylor1) ≈ 5.0
vectFloats1=[1.0,2.0,3.0]
@test dot(vectFloats1,vectTaylor1) ≈ 8.0
@test dot(vectTaylor1,vectFloats1) ≈ 8.0
vectFloats2=[3.0,4.0]
vectTaylor3=[t1,t2]
@test *(matrixTaylor,vectFloats2) == [5.0, 6.0]
@test *(matrixTaylor,vectTaylor3) == [1.25, 2.0]
@test cross(vectTaylor1,vectTaylor2) == [-1.0, -1.0, 0.75]
@test cross(vectTaylor2,vectTaylor1) == [1.0, 1.0, -0.75]
@test cross(vectFloats1,vectTaylor1) == [2.5, 1.0, -1.5]
@test cross(vectTaylor1,vectFloats1) == [-2.5, -1.0, 1.5]



