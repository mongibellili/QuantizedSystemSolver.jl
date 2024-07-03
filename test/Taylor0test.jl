using QuantizedSystemSolver
#using Test
t1=Taylor0([1.0,1.0,0.0],2)
t2=Taylor0([2.0,2.0,3.0],2)
t3=Taylor0([1.0,1.0,0.0],2)
cache1=Taylor0([1.0,1.0,1.0],2)
cache2=Taylor0([1.0,1.0,1.0],2)
@test zero(t2)==Taylor0([0.0,0.0,0.0],2)
@test one(t2)==Taylor0([1.0,0.0,0.0],2)
@test t1==t3
zero(t1,cache1)
@test cache1[0]==0.0
one(t1,cache1)
@test cache1[0]==1.0
t4=t1+t3
@test t4==Taylor0([2.0,2.0,0.0],2)
@test 5.0+t4==Taylor0([7.0,2.0,0.0],2)
@test t4+5.0==Taylor0([7.0,2.0,0.0],2)
@test t4-t3==t1
@test -t4==Taylor0([-2.0,-2.0,0.0],2)
@test t4-5.0==Taylor0([-3.0,2.0,0.0],2)
@test 5.0-t4==Taylor0([3.0,-2.0,-0.0],2)
@test 5.0*t4==Taylor0([10.0,10.0,0.0],2)
@test t4/2.0==Taylor0([1.0,1.0,0.0],2)
@test t1*t2 ==  Taylor0([2.0, 4.0, 5.0], 2)
@test t4/t1 ==Taylor0([2.0, 0.0, 0.0], 2)
createT(3.6,cache1)
@test cache1[0]==3.6
createT(t2,cache1)
@test cache1[0]==2.0
addT(4.3,5.2,cache1)
@test cache1[0]==9.5
addT(t1,5.6,cache1)
@test cache1[0]==6.6
addT(t1,t2,cache1)
@test cache1[0]==3.0
addT(t1,t2,4.8,cache1)
@test cache1[0]==7.8
addT(t1,t2,t3,cache1)
@test cache1[0]==4.0
addT(t1,t2,t3,1.2,cache1)
@test cache1[0]==5.2
#mulT,mulTT,addsub,negateT,subsub,subadd,subT,muladdT,mulsub,divT