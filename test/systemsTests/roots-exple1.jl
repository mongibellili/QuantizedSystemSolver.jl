
using Polynomials
using Roots
p1=Polynomial([-2,1,1])
p2=Polynomial([-2,1,2])
p3=Polynomial([-2,1,2])
#------------ tutor1 find_zeros and find_zero-------------from the roots packagae
#real sols are -2 and 1"""
#s=find_zeros(p1, 0, 2) #[1.0] # 
#s=find_zeros(p1, -3, 0)#[-2.0]
#s=find_zeros(p1, -3, 2)#[-2.0, 1.0]
#s=find_zero(p1, 0)  #1.0
#s2=find_zeros(p2, -3, 1) #[-1.2807764064044151, 0.7807764064044151]
#s2=find_zero(p2, 0)  #0.7807764064044151
#= s2=find_zeros(p2, 0, 1)#[0.7807764064044151]
println(s2) =#
#------------- tutor2 fzeros and fzero-------------"""
s=fzeros(p1, 0.5, 2) #[1.0]
println(s)
#s=fzeros(p1, -3, 0)#[-2.0]
#s=fzeros(p1, -3, 2)#[-2.0, 1.0]
#s=fzero(p1, 0)  #1.0
#s2=fzeros(p2, -3, 1) #[-1.2807764064044151, 0.7807764064044151]
#s2=fzero(p2, 0)  #0.7807764064044151
#= s2=fzeros(p2, 0, 1)#[0.7807764064044151]
println(s2) =#
#------------- tutor3 roots-------------just from the Polynomials package
#= a = roots(p1)
println(a)
b = roots(p2)
println(b)
r = filter(v->abs(imag(v)) < 1.0e-15 && real(v)>=0.0, [a..., b...])
println(r) =#
println( ChebyshevT([1, 0, 3, 4]))