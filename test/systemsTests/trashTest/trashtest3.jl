

#= counterRootj=1
for h in (resj1,resj2,resj3,resj4,resj5,resj6,resj7,resj8,resj9,resj10,resj11,resj12)
  if h>0
    cacheRootsj[counterRootj]=h
    counterRootj+=1
  end
end
sort!(cacheRootsj);


(coefi2[1]*h^6+coefi2[2]*h^5+coefi2[3]*h^4+coefi2[4]*h^3+coefi2[5]*h^2+coefi2[6]*h+coefi2[7]) =#
#= a=BigFloat(12.05e1)
@show a
b=0.0005+a
@show b
c=Float64(b)
@show c =#

#= a=1/0
b=0/0
@show a,b =#

#= coeffs=NTuple{7,Float64}((29160.956861496, 67.56376290717117, 0.014452560693316113,56.2,1250.12,565.21,0.002))
coeff=coeffs[2:end]
@show coeff
@show typeof(coeff) =#

#@show BigFloat <: Real
#= using InteractiveUtils
P=Union{Float64,Int64}
function test(x::P)
  @show x
  @show typeof(x)
end
y=12
display(@code_typed test(y))
function test2(x::Float64)
  @show x
  @show typeof(x)
end
z=12.2
display(@code_typed test(z)) =#

#= coefH2=2.0*(-aii*uij-aij*uji-uij2)

coefH3=2.0*((ajj*uij2-aij*uji2))
coefH4=0.0
coefH5=-ajj*aiijj_*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)-aij*aiijj*(aii*uji2-uji*aiijj_-aji*uij2)+aij*aiijj_*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_)
coefH6=ajj*aiijj_*(ajj*uij2-uij*aiijj_-aij*uji2)-aij*aiijj_*(aii*uji2-uji*aiijj_-aji*uij2) =#
#= aii=ajj=aij=aji=uij=uji=uij2=uji2=xi=xjaux=2.0
aiijj=aiijj_=0.0
#= coefH2=2.0*(-aii*uij-aij*uji-uij2+xi*(-3.0*aiijj_+2.0*aiijj*aiijj+aiijj*ajj)-aij*aiijj*xjaux)
#coefH3=2.0*((-uij*aiijj_+ajj*uij2-aij*uji2)+aiijj*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)-xi*aiijj*(2.0*aii*ajj-aji*aij+ajj*ajj)+ajj*aiijj_*xi+aij*aiijj*aiijj*xjaux-aij*aiijj_*xjaux)
coefH3=2.0*((-uij*aiijj_+ajj*uij2-aij*uji2)+aiijj*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)-xi*aiijj*(ajj*aiijj-aiijj_)+ajj*aiijj_*xi+aij*aiijj*aiijj*xjaux-aij*aiijj_*xjaux)
coefH4=2.0*(aiijj*(uij*aiijj_-ajj*uij2+aij*uji2)-0.5*(ajj*aiijj-aiijj_)*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)-ajj*aiijj_*aiijj*xi+0.5*aij*aiijj*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_)+aij*aiijj_*aiijj*xjaux)
coefH5=(ajj*aiijj-aiijj_)*(ajj*uij2-uij*aiijj_-aij*uji2)-ajj*aiijj_*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)-aij*aiijj*(aii*uji2-uji*aiijj_-aji*uij2)+aij*aiijj_*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_)
coefH6=ajj*aiijj_*(ajj*uij2-uij*aiijj_-aij*uji2)-aij*aiijj_*(aii*uji2-uji*aiijj_-aji*uij2) =#

coefH2j=2.0*(-ajj*uji-aji*uij-uji2+xjaux*(-2.0*aiijj_+2.0*aiijj*aiijj+aii*aiijj-aiijj_)-aji*aiijj*xi)
coefH3j=2.0*((-uji*aiijj_+aii*uji2-aji*uij2)+aiijj*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_)-xjaux*aiijj*(aii*aiijj-aiijj_)+aii*aiijj_*xjaux+aji*aiijj*aiijj*xi-aji*aiijj_*xi)
coefH4j=2.0*(aiijj*(uji*aiijj_-aii*uji2+aji*uij2)-0.5*(aii*aiijj-aiijj_)*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_)-aii*aiijj_*aiijj*xjaux+0.5*aji*aiijj*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)+aji*aiijj_*aiijj*xi)
coefH5j=(aii*aiijj-aiijj_)*(aii*uji2-uji*aiijj_-aji*uij2)-aii*aiijj_*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_)-aji*aiijj*(ajj*uij2-uij*aiijj_-aij*uji2)+aji*aiijj_*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)
coefH6j=aii*aiijj_*(aii*uji2-uji*aiijj_-aji*uij2)-aji*aiijj_*(ajj*uij2-uij*aiijj_-aij*uji2)


@show coefH2j,coefH3j,coefH4j,coefH5j,coefH6j =#

#= using PolynomialRoots
setprecision(BigFloat,1024) =#
#= using Polynomials
using Roots =#
#40.80630971014806897771247, 178.48731823943597794466468, 0.0001425406365094045,
#coefi=[(-0.00057016254603761695740615778), (-1140.3250641304001419575054), (-8.5524376124873223589469262e+08), BigFloat(-2.8508122708869649897226691e+14), (1.7114669426922051413515583e+13), (1.0474193043858423245087585e+18), (-2.5668078295588949092296481e+16)]
#coefi=[BigFloat(-0.00057016254603761695740615778), BigFloat(-1140.3250641304001419575054), BigFloat(-8.5524376124873223589469262e+08), BigFloat(-2.8508122708869649897226691e+14), BigFloat(1.7114669426922051413515583e+13), BigFloat(1.0474193043858423245087585e+18), BigFloat(-2.5668078295588949092296481e+16)]
#coefi=BigFloat[(-0.00057016254603761695740615778), (-1140.3250641304001419575054), (-8.5524376124873223589469262e+08), (-2.8508122708869649897226691e+14), (1.7114669426922051413515583e+13), (1.0474193043858423245087585e+18), (-2.5668078295588949092296481e+16)]
#coefi=BigFloat[-0.00057016254603761695740615778, -1140.3250641304001419575054, -8.5524376124873223589469262e+08, -2.8508122708869649897226691e+14, 1.7114669426922051413515583e+13, 1.0474193043858423245087585e+18, -2.5668078295588949092296481e+16]
#@show coefi
#= 
y=BigFloat(-1140.3250641304)
@show y =#
#p=Polynomial(coefi)
 #@show p
#=  r=roots(coefi)
 #r=find_zeros(p, BigFloat(0.0), BigFloat(100.0))
 #@show r
#= r1=[Float64(real(a)) for a in r]
r2=[Float64(imag(a)) for a in r] =#

#= r1=map(a->Float64(real(a)),r)
r2=map(a->Float64(imag(a)),r) =#
#= @show r1
@show r2 =#
#@show r
 r3=[]
for h in r
  if abs(imag(h)) < 1.0e-15 && real(h)>0.0
    push!(r3,real(h))
  end
end
@show r3
for h in r3
  if h>2.0
  x=h =#
#=   h=x=BigFloat(40.80630971014806897771247)
  zci=coefi[7]*h^6+coefi[6]*h^5+coefi[5]*h^4+coefi[4]*h^3+coefi[3]*h^2+coefi[2]*h+coefi[1]
  @show h,zci,coefi
  @show typeof(h),typeof(zci),typeof(coefi) =#
#=   zc1=coefi[3]*x^2+coefi[4]*x^3+coefi[5]*x^4+coefi[6]*x^5+coefi[7]*x^6+coefi[1]+coefi[2]*x
  #= x=40.80630971014807
  x=40.806309710148064 =#
  
  zc3=coefi[1]+coefi[2]*x+coefi[3]*x^2+coefi[4]*x^3+coefi[5]*x^4+coefi[6]*x^5+coefi[7]*x^6
  zc4=(coefi[1]+coefi[2]*x+coefi[3]*x^2)+coefi[4]*x^3+coefi[5]*x^4+(coefi[6]*x^5+coefi[7]*x^6)
  #@show coefi[1],coefi[2],coefi[3],coefi[4],coefi[5],coefi[6],coefi[7]

  zc2=-0.000570162546037617 - 1140.3250641304*x - 8.552437612487322e8*x^2 - 2.850812270886965e14*x^3 + 1.711466942692205e13*x^4 + 1.0474193043858423e18*x^5 - 2.5668078295588948e16*x^6
 =#
  
#=   end
end  =#
#40.80630971014806897771247, 178.48731823943597794466468, 0.0001425406365094045, BigFloat[-0.00057016254603761695740615778, -1140.3250641304001419575054, -8.5524376124873223589469262e+08, -2.8508122708869649897226691e+14, 1.7114669426922051413515583e+13, 1.0474193043858423245087585e+18, -2.5668078295588949092296481e+16]



ex=:(0*(t/6))
dump(ex)