#= 
setprecision(BigFloat,100)

aii, aij, aji, ajj, h, xi, xjaux, uij, uij2, uji, uji2 = BigFloat(-1.0e6), BigFloat(1000.0), BigFloat(1.0e6), BigFloat(-1000.0157261488598), BigFloat(8.858170228361118), BigFloat(0.0009738922631926892), BigFloat(0.9538344383858633), BigFloat(0.024830940874835505), BigFloat(-0.02117642588838642), BigFloat(-5.135234459885396e-10), BigFloat(0.00014385727589605324)

#aii, aij, aji, ajj, h, xi, xjaux, uij, uij2, uji, uji2 = -1.0e6, 1000.0, 1.0e6, -1000.0157261488598, 8.858170228361118, 0.0009738922631926892, 0.9538344383858633, 0.024830940874835505, -0.02117642588838642, -5.135234459885396e-10, 0.00014385727589605324
h_2=h*h;h_3=h_2*h;h_4=h_3*h;h_5=h_4*h;h_6=h_5*h
aiijj=aii+ajj
aiijj_=aij*aji-aii*ajj
#Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
Δ1=1.0-h*(aiijj)-h_2*(aiijj_)
if abs(Δ1)==0.0
  Δ1=1e-30
  @show Δ1
end

#=  qpar2=(h_2*aij*aiijj+h_3*aij*aiijj_)*(xjaux-h*xjaux*aiijj-0.5*h_2*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_)+0.5*h_3*(aii*uji2-uji*aiijj_-aji*uij2))
 qpar21=(h_2*aij*aiijj)*(xjaux-h*xjaux*aiijj-0.5*h_2*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_)+0.5*h_3*(aii*uji2-uji*aiijj_-aji*uij2))
 qpar22=(h_3*aij*aiijj_)*(xjaux-h*xjaux*aiijj-0.5*h_2*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_)+0.5*h_3*(aii*uji2-uji*aiijj_-aji*uij2))
 qpar_2=h_2*aij*aiijj*xjaux-h_3*aij*aiijj*aiijj*xjaux-0.5*h_4*aij*aiijj*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_)+0.5*h_5*aij*aiijj*(aii*uji2-uji*aiijj_-aji*uij2)+h_3*aij*aiijj_*xjaux-h_4*aij*aiijj_*xjaux*aiijj-0.5*h_5*aii*aiijj_*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_)+0.5*h_6*aij*aiijj_*(aii*uji2-uji*aiijj_-aji*uij2)
 qpar_21=h_2*aij*aiijj*xjaux-h_3*aij*aiijj*aiijj*xjaux-0.5*h_4*aij*aiijj*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_)+0.5*h_5*aij*aiijj*(aii*uji2-uji*aiijj_-aji*uij2)
  =#
 
 qpar22=(h_3*aij*aiijj_)*(xjaux-h*xjaux*aiijj-0.5*h_2*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_)+0.5*h_3*(aii*uji2-uji*aiijj_-aji*uij2))
 qpar221=(h_3*aij*aiijj_)*(xjaux-h*xjaux*aiijj-0.5*h_2*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_))
 #qpar222=(h_3*aij*aiijj_)*(xjaux-h*xjaux*aiijj-0.5*h_2*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_)+0.5*h_3*(aii*uji2-uji*aiijj_-aji*uij2))
 qpar_22=h_3*aij*aiijj_*xjaux-h_4*aij*aiijj_*xjaux*aiijj-0.5*h_5*aii*aiijj_*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_)+0.5*h_6*aij*aiijj_*(aii*uji2-uji*aiijj_-aji*uij2)
 qpar_221=h_3*aij*aiijj_*xjaux-h_4*aij*aiijj_*xjaux*aiijj-0.5*h_5*aij*aiijj_*(ajj*uji+aji*uij+uji2+2.0*xjaux*aiijj_)
 @show qpar221,qpar_221 =#