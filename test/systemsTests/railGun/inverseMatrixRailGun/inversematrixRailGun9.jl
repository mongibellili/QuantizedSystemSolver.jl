using Symbolics, LinearAlgebra
import Base.:+
+(a::T, b::Vector{T}) where {T <:Union{Float64,Num}}=a+b[1]
@show 2
@variables x,il1,il2,il3,il4,il5,il6,il7,il8,il9,rd1,rd2,rd3,rd4,rd5,rd6,rd7,rd8,rd9,rr,rpr,v,is1,is2,is3,is4,is5,is6,is7,is8,is9,Il
b=4.2e-3+0.453*x;a=5.1+4.2e-3+0.453*x  ;λ=a-b;u=[b;b;b;b;b;b;b;b;b];vec=[1 1 1 1 1 1 1 1 1]
mm=I/λ-(u*vec)/(λ*(λ+vec*u))
#display(mm) 
C=[-(rd1+3.88e-3)*il1-(rr+rpr+0.453e-6*v)*Il+rd1*is1;-(rd2+3.88e-3)*il2-(rr+rpr+0.453e-6*v)*Il+rd2*is2;-(rd3+3.88e-3)*il3-(rr+rpr+0.453e-6*v)*Il+rd3*is3;-(rd4+3.88e-3)*il4-(rr+rpr+0.453e-6*v)*Il+rd4*is4;-(rd5+3.88e-3)*il5-(rr+rpr+0.453e-6*v)*Il+rd5*is5;-(rd6+3.88e-3)*il6-(rr+rpr+0.453e-6*v)*Il+rd6*is6;-(rd7+3.88e-3)*il7-(rr+rpr+0.453e-6*v)*Il+rd7*is7;-(rd8+3.88e-3)*il8-(rr+rpr+0.453e-6*v)*Il+rd8*is8;-(rd9+3.88e-3)*il9-(rr+rpr+0.453e-6*v)*Il+rd9*is9]
#display(C)
D=mm*C
display(D) 
