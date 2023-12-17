#= using SymEngine
x=4.0
d=[1.0,1.0,2.0,3.0]
#m=quote 3+5x+d*x end

macro test(m)
Base.remove_linenums!(m)
str=string(m)
@show str
basi = convert(Basic, m)
coef = diff(basi, :x)
@show coef
end

@test 3+5x+d[g]*x =#
#= using Plots
function test()
q(h)=-7.961146634883693-32887.10978074606*h
x(h)=-7.9611478450843896-32887.10978145532*h+0.001082080416381359*h*h
#= xd=x(0.033857263959124986)
qd=q(0.033857263959124986)
@show qd-xd =#
xd=x(0.0003)
qd=q(0.000000)
@show qd-xd


#display(plot!(x))
#= display(plot!(x,xlims=(0.0,0.034) ,ylims=(-7.9612,-7.96114)))
#display(plot!(x,xlims=(0.0,1e-9) ,ylims=(-7.9612,-7.96114)))
display(plot!(q))

println("press enter to exit")
    readline()  =#
end
test() =#
#= using qss
import qss: sin
t = Taylor0(zeros(3), 2)
    t[1]=1.0
   # t[0]=1e-5
t1=60.0*t
@show sin(t) =#
#= bb=1010.4995370527463
cc=184882.09294849282
delta=6032.086629145888
@show (bb+cc)/delta
@show bb/delta+cc/delta =#
#setprecision(BigFloat,80)
#= delta=BigFloat(1.135700289678955e7)
h=9.852699092995836
uij, uij2 = 0.057013364027397984, -0.04983356413099216
d=h * uij + (h * h * uij2) / 2
@show d
d2=delta*(h * uij + (h * h * uij2) / 2)

@show d2

res=d2/delta
d2_=Float64(res)
@show d2_
 =#
#=  aii,aij,aji,ajj,h,xi,xjaux,uij,uij2,uji,uji2= -1.0e6, 1000.0, 1.0e6, -1000.0020103844182, 24.999993861351562, 0.0007398552897831782, 0.7492601358736745, 1.0113325856764277e-5, 0.01846103323623538, 6.396817298082169e-9, -1.2535862317308784
 aii,aij,aji,ajj,h,xi,xjaux,uij,uij2,uji,uji2=BigFloat(aii),BigFloat(aij),BigFloat(aji),BigFloat(ajj),BigFloat(h),BigFloat(xi),BigFloat(xjaux),BigFloat(uij),BigFloat(uij2),BigFloat(uji),BigFloat(uji2)
 Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji


 αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
 # αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
  αij=aij/Δ1
  αji=aji/Δ1
  αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1


  βji_=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
  βji=h*h*aji*(aii+ajj)/(2*Δ1)+h*h*h*(aji*aji*aij-aji*ajj*aii)/(2*Δ1)


@show βji_,βji, abs(βji_-βji) =#



#=  λii=(h*h*aii/2-h)*(1-h*ajj)+h*h*h*aji*aij/2
 λij=(h*h*aii/2-h)*h*aij+h*h*aij*(1-h*aii)/2
 λji=h*h*aji/2*(1-h*ajj)+(h*h*ajj/2-h)*h*aji
 λjj=h*h*h*aij*aji/2+(h*h*ajj/2-h)*(1-h*aii)

 parti_=((λii*(uij+h*uij2)+λij*(uji+h*uji2))/Δ1)+(xi+h*uij+h*h*uij2/2)#part1[1]+xpart2[1]#
  parti=((xi-h*xi*(aii+ajj))-h*h*((aii)*uij+aij*uji+uij2-2*xi*(aii*ajj-aij*aji))/2+h*h*h*(uij*(aii*ajj-aij*aji)+uij2*ajj-aij*uji2)/2)/Δ1
 @show parti_,parti =#


#=  x=1e20*(1.0+1e-16)-1e20
 y=1e20+1e4-1e20
 z=1e4
 @show x,y,z =#

#=  xx=BigFloat(1e20*(1.0+1e-16)-1e20)
 yy=BigFloat(1e20+1e4-1e20)
 zz=BigFloat(1e4)
 @show xx,yy,zz =#

#=  xx=BigFloat(1e20)*(BigFloat(1.0)+BigFloat(1e-16))-BigFloat(1e20)
 yy=BigFloat(1e20)+BigFloat(1e4)-BigFloat(1e20)
 zz=BigFloat(1e4)
 @show xx,yy,zz =#
#= # =#

a=0.0

c=0.0

d=a/c
@show d
@show isnan(d)