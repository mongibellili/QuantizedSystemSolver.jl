

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
@test +t1==Taylor0([1.0,1.0,0.0],2)
@test t4==Taylor0([2.0,2.0,0.0],2)
@test 5.0+t4==Taylor0([7.0,2.0,0.0],2)
@test t4+5.0==Taylor0([7.0,2.0,0.0],2)
@test t4-t3==t1
@test -t4==Taylor0([-2.0,-2.0,0.0],2)
@test t4-5.0==Taylor0([-3.0,2.0,0.0],2)
@test 5.0-t4==Taylor0([3.0,-2.0,-0.0],2)
@test 5.0*t4==Taylor0([10.0,10.0,0.0],2)
@test t4*5.0==Taylor0([10.0,10.0,0.0],2)
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
addT(5.6,t1,cache1)
@test cache1[0]==6.6
addT(t1,t2,cache1)
@test cache1[0]==3.0
addT(t1,t2,4.8,cache1)
@test cache1[0]==7.8
addT(1.0,1.0,4.8,cache1)
@test cache1[0]==6.8
addT(t1,4.8,t2,cache1)
@test cache1[0]==7.8
addT(4.8,t2,t1,cache1)
@test cache1[0]==7.8
addT(4.8,t2,1.0,cache1)
@test cache1[0]==7.8
addT(4.8,1.0,t2,cache1)
@test cache1[0]==7.8
addT(t1,t2,t3,cache1)
@test cache1[0]==4.0
addT(t1,t2,t3,1.2,cache1)
@test cache1[0]==5.2
addT(t1,t2,t3,1.2,1.0,cache1)
@test cache1[0]==6.2
addT(t1,t2,t3,1.2,1.0,1.0,cache1)
@test cache1[0]==7.2
addT(t1,t2,t3,1.2,1.0,1.0,t1,cache1)
@test cache1[0]==8.2
addT(t1,t2,t3,1.2,1.0,1.0,t1,1.0,cache1)
@test cache1[0]==9.2
addT(t1,t2,t3,1.2,1.0,1.0,t1,1.0,0.3,cache1)
@test cache1[0]==9.5
addT(t1,0.5,t2,t3,1.2,1.0,1.0,t1,1.0,0.3,cache1)
@test cache1[0]==10.0
subT(4.3,5.3,cache1)
@test cache1[0]==-1.0

subT(t2,t1,cache1)
@test cache1==Taylor0([1.0, 1.0, 3.0], 2)
cache1=Taylor0([0.0,0.0,0.0],2)
subT(2.0,t1,cache1)
@test cache1==Taylor0([1.0, -1.0, 0.0], 2)
subT(t1,2.0,cache1)
@test cache1==Taylor0([-1.0, 1.0, 0.0], 2)
negateT(t1,cache1)
@test cache1==Taylor0([-1.0, -1.0, 0.0], 2)
cache1=Taylor0([0.0,0.0,0.0],2)
negateT(2.0,cache1)
@test cache1==Taylor0([-2.0, 0.0, 0.0], 2)
addsub(5.0,1.5,2.0,cache1)
@test cache1[0]==4.5
addsub(t1,t2,t3,cache1)
@test cache1[0]==2.0
addsub(t1,5.0,t2,cache1)
@test cache1[0]==4.0
addsub(0.5,5.0,t2,cache1)
@test cache1[0]==3.5
addsub(t2,1.0,0.5,cache1)
@test cache1[0]==2.5
addsub(1.0,t2,0.5,cache1)
@test cache1[0]==2.5

subadd(5.0,1.5,2.0,cache1)
@test cache1[0]==5.5

subsub(5.0,1.5,2.0,cache1)
@test cache1[0]==1.5
subsub(t1,t2,t3,cache1)
@test cache1[0]==-2.0
subsub(t1,5.0,t2,cache1)
@test cache1[0]==-6.0
subsub(5.0,t1,t2,cache1)
@test cache1[0]==2.0
subsub(0.5,5.0,t2,cache1)
@test cache1[0]==-6.5
subsub(t1,1.0,1.0,cache1)
@test cache1[0]==-1.0
#mulT
mulT(4.3,5.2,cache1)
@test cache1[0]==22.36
mulT(t1,5.6,cache1)
@test cache1[0]==5.6
mulT(t1,t2,cache1)
@test cache1[0]==2.0

#mulTT
mulTT(5.0,1.5,2.0,cache1,cache2)
@test cache1[0]==15.0
mulTT(t1,t2,t3,cache1,cache2)
@test cache1[0]==2.0
mulTT(t1,5.0,t2,cache1,cache2)
@test cache1[0]==10.0
mulTT(0.5,5.0,t2,cache1,cache2)
@test cache1[0]==5.0
mulTT(0.5,5.0,t2,1.0,cache1,cache2)
@test cache1[0]==5.0
mulTT(0.5,5.0,t2,1.0,t1,cache1,cache2)
@test cache1[0]==5.0
mulTT(0.5,5.0,t2,1.0,t1,0.5,cache1,cache2)
@test cache1[0]==2.5
mulTT(0.5,5.0,t2,1.0,t1,0.5,1.0,cache1,cache2)
@test cache1[0]==2.5

#muladdT
muladdT(5.0,1.5,2.0,cache1)
@test cache1[0]==9.5
muladdT(t1,t2,t3,cache1)
@test cache1[0]==3.0
muladdT(t1,5.0,t2,cache1)
@test cache1[0]==7.0
muladdT(0.5,5.0,t2,cache1)
@test cache1[0]==4.5
#mulsub
mulsub(5.0,1.5,2.0,cache1)
@test  cache1[0]==5.5
mulsub(t1,t2,t3,cache1)
@test  cache1[0]==1.0
mulsub(t1,5.0,t2,cache1)
@test  cache1[0]==3.0
mulsub(0.5,5.0,t2,cache1)
@test  cache1[0]==0.5

#divT
divT(4.3,2.0,cache1)
@test  cache1[0]==2.15
divT(t1,2.0,cache1)
@test  cache1[0]==0.5
divT(5.6,t1,cache1)
@test  cache1[0]==5.6
divT(t1,t2,cache1)
@test  cache1[0]==0.5