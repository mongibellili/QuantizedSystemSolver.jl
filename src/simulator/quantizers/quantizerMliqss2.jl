"""
    isCycle_simulUpdate(aii::Float64,ajj::Float64,aij::Float64,aji::Float64,trackSimul,::Val{2},::Val{M},index::Int,j::Int,dirI::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},dxaux::Vector{MVector{2,Float64}},qaux::Vector{MVector{2,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64) where {M}

Performs a simultaneous update of two quantized variables `qi` and `qj` and their derivatives if cycle conditions are met.
This is similar to the [`isCycle_simulUpdate`](@ref) order 1 function.
"""
function isCycle_simulUpdate(aii::Float64,ajj::Float64,aij::Float64,aji::Float64,trackSimul,::Val{2},::Val{M},index::Int,j::Int,dirI::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},dxaux::Vector{MVector{2,Float64}},qaux::Vector{MVector{2,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64) where {M}

  if isnan(aii)
    @warn("a is NaN: The Jacobian is not defined at this instant $(simt). This may be due to  an undefined operation. Consider computing the Jacobian coefficient manually by Explicitly writing jac_mode = :approximate in ODEProblem.")
    aii = 0.0
  end
  if isnan(ajj)
    @warn("a is NaN: The Jacobian is not defined at this instant $(simt). This may be due to  an undefined operation. Consider computing the Jacobian coefficient manually by Explicitly writing jac_mode = :approximate in ODEProblem.")
    ajj = 0.0
  end
  if isnan(aij)
    @warn("a is NaN: The Jacobian is not defined at this instant $(simt). This may be due to  an undefined operation. Consider computing the Jacobian coefficient manually by Explicitly writing jac_mode = :approximate in ODEProblem.")
    aij = 0.0
  end
  if isnan(aji)
    @warn("a is NaN: The Jacobian is not defined at this instant $(simt). This may be due to  an undefined operation. Consider computing the Jacobian coefficient manually by Explicitly writing jac_mode = :approximate in ODEProblem.")
    aji = 0.0
  end
  uii=dxaux[index][1]-aii*qaux[index][1]
  ui2=dxaux[index][2]-aii*qaux[index][2]
  xi=x[index][0];xj=x[j][0];qi=q[index][0];qj=q[j][0];qi1=q[index][1];qj1=q[j][1];xi1=x[index][1];xi2=2*x[index][2];xj1=x[j][1];xj2=2*x[j][2]
  quanj=quantum[j];quani=quantum[index];
  e1 = simt - tx[j];e2 = simt - tq[j];
  prevXj=x[j][0]
  x[j][0]= x[j](e1);xj=x[j][0];tx[j]=simt
  recentjDir=(xj-prevXj)
  qj=qj+e2*qj1  ;qaux[j][1]=qj;tq[j] = simt    ;q[j][0]=qj  
  xj1=x[j][1]+e1*xj2;
  newDiff=(xj-prevXj)
  ujj=xj1-ajj*qj
  uji=ujj-aji*qaux[index][1]# 
  uj2=xj2-ajj*qj1###################################################-----------------------
  uji2=uj2-aji*qaux[index][2]#
  dxj=aji*qi+ajj*qaux[j][1]+uji
  ddxj=aji*qi1+ajj*qj1+uji2
  if abs(ddxj)==0.0 ddxj=1e-30 end
  iscycle=false
  qjplus=xj-sign(ddxj)*quanj
  hj=sqrt(2*quanj/abs(ddxj))#
  α1=1-hj*ajj
  if abs(α1)==0.0 α1=1e-30 end
  dqjplus=(aji*(qi+hj*qi1)+ajj*qjplus+uji+hj*uji2)/α1
  uij=uii-aij*qaux[j][1]
  uij2=ui2-aij*qj1#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
  dxithrow=aii*qi+aij*qj+uij
  ddxithrow=aii*qi1+aij*qj1+uij2
  dxi=aii*qi+aij*qjplus+uij
  ddxi=aii*qi1+aij*dqjplus+uij2
  if abs(ddxi)==0.0 ddxi=1e-30 end
  hi=sqrt(2*quani/abs(ddxi))
  βidir=dxi+hi*ddxi/2
  βjdir=dxj+hj*ddxj/2
  βidth=dxithrow+hi*ddxithrow/2
  αidir=xi1+hi*xi2/2
  βj=xj1+hj*xj2/2
  ########cycle detection condition

  iscycle=detect2(Val(M),xi1,dxi,dxithrow,xi2,ddxi,ddxithrow,βidir,βidth,xj1,dxj,xj2,ddxj,dqjplus,recentjDir,dirI) 

  if iscycle
        aiijj=aii+ajj
        aiijj_=aij*aji-aii*ajj
        coefΔh1=-8.0*aiijj
        coefΔh2=6.0*aiijj*aiijj-4.0*aiijj_
        coefΔh3=6.0*aiijj*aiijj_-2.0*aiijj*aiijj*aiijj
        coefΔh4=aiijj_*aiijj_-4.0*aiijj_*aiijj*aiijj
        coefΔh5=-3.0*aiijj*aiijj_*aiijj_
        coefΔh6=-aiijj_*aiijj_*aiijj_
        coefH1=-8.0*xi*aiijj
        coefH2=2.0*(-aii*uij-aij*uji-uij2+xi*(-2.0*aiijj_+2.0*aiijj*aiijj+2.0*aii*ajj-aji*aij+ajj*ajj)-aij*aiijj*xj)
        coefH3=2.0*((-uij*aiijj_+ajj*uij2-aij*uji2)+aiijj*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)-xi*aiijj*(2.0*aii*ajj-aji*aij+ajj*ajj)+ajj*aiijj_*xi+aij*aiijj*aiijj*xj-aij*aiijj_*xj)
        coefH4=2.0*(aiijj*(uij*aiijj_-ajj*uij2+aij*uji2)-0.5*(2.0*aii*ajj-aji*aij+ajj*ajj)*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)-ajj*aiijj_*aiijj*xi+0.5*aij*aiijj*(ajj*uji+aji*uij+uji2+2.0*xj*aiijj_)+aij*aiijj_*aiijj*xj)
        coefH5=(2.0*aii*ajj-aji*aij+ajj*ajj)*(ajj*uij2-uij*aiijj_-aij*uji2)-ajj*aiijj_*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)-aij*aiijj*(aii*uji2-uji*aiijj_-aji*uij2)+aij*aiijj_*(ajj*uji+aji*uij+uji2+2.0*xj*aiijj_)
        coefH6=ajj*aiijj_*(ajj*uij2-uij*aiijj_-aij*uji2)-aij*aiijj_*(aii*uji2-uji*aiijj_-aji*uij2)
        coefH1j=-8.0*xj*aiijj
        coefH2j=2.0*(-ajj*uji-aji*uij-uji2+xj*(-2.0*aiijj_+2.0*aiijj*aiijj+2.0*ajj*aii-aij*aji+aii*aii)-aji*aiijj*xi)
        coefH3j=2.0*((-uji*aiijj_+aii*uji2-aji*uij2)+aiijj*(ajj*uji+aji*uij+uji2+2.0*xj*aiijj_)-xj*aiijj*(2.0*ajj*aii-aij*aji+aii*aii)+aii*aiijj_*xj+aji*aiijj*aiijj*xi-aji*aiijj_*xi)
        coefH4j=2.0*(aiijj*(uji*aiijj_-aii*uji2+aji*uij2)-0.5*(2.0*ajj*aii-aij*aji+aii*aii)*(ajj*uji+aji*uij+uji2+2.0*xj*aiijj_)-aii*aiijj_*aiijj*xj+0.5*aji*aiijj*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)+aji*aiijj_*aiijj*xi)
        coefH5j=(2.0*ajj*aii-aij*aji+aii*aii)*(aii*uji2-uji*aiijj_-aji*uij2)-aii*aiijj_*(ajj*uji+aji*uij+uji2+2.0*xj*aiijj_)-aji*aiijj*(ajj*uij2-uij*aiijj_-aij*uji2)+aji*aiijj_*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)
        coefH6j=aii*aiijj_*(aii*uji2-uji*aiijj_-aji*uij2)-aji*aiijj_*(ajj*uij2-uij*aiijj_-aij*uji2)
        h = ft-simt
        Δ1=1.0-h*(aiijj)-h*h*(aiijj_)
        
        if abs(Δ1)!=0.0 
          h_2=h*h;h_3=h_2*h;h_4=h_3*h;h_5=h_4*h;h_6=h_5*h
          Δ22=4.0+h*coefΔh1+h_2*(coefΔh2)+h_3*(coefΔh3)+h_4*(coefΔh4)+h_5*coefΔh5+h_6*coefΔh6
          if abs(Δ22)!=0.0
            qi=(4.0*xi+h*coefH1+h_2*coefH2+h_3*coefH3+h_4*coefH4+h_5*coefH5+h_6*coefH6)/Δ22
            qj=(4.0*xj+h*coefH1j+h_2*coefH2j+h_3*coefH3j+h_4*coefH4j+h_5*coefH5j+h_6*coefH6j)/Δ22
          end
        end

        if (abs(qi - xi) > 2.0*quani || abs(qj - xj) > 2.0*quanj) 
          h1 = sqrt(abs(2*quani/xi2));
          h2 = sqrt(abs(2*quanj/xj2));   #later add derderX =1e-12 when x2==0?

          h=min(h1,h2)
          if h!=Inf
              Δ1=1.0-h*(aiijj)-h*h*(aiijj_)
              if abs(Δ1)!=0.0 
                h_2=h*h;h_3=h_2*h;h_4=h_3*h;h_5=h_4*h;h_6=h_5*h
                Δ22=4.0+h*coefΔh1+h_2*(coefΔh2)+h_3*(coefΔh3)+h_4*(coefΔh4)+h_5*coefΔh5+h_6*coefΔh6
                if abs(Δ22)!=0.0
                  qi=(4.0*xi+h*coefH1+h_2*coefH2+h_3*coefH3+h_4*coefH4+h_5*coefH5+h_6*coefH6)/Δ22
                  qj=(4.0*xj+h*coefH1j+h_2*coefH2j+h_3*coefH3j+h_4*coefH4j+h_5*coefH5j+h_6*coefH6j)/Δ22
                end
              end
          end

        end

        maxIter=10
        while (abs(qi - xi) > 2.0*quani || abs(qj - xj) > 2.0*quanj) && (maxIter>0)
          maxIter-=1
          if maxIter>1  #ie simul step failed
            return false 
          end
          h1 = h * sqrt(quani / abs(qi - xi));
          h2 = h * sqrt(quanj / abs(qj - xj));
          h=min(h1,h2)
          if h!=Inf
   
              Δ1=1.0-h*(aiijj)-h*h*(aiijj_)
              if abs(Δ1)!=0.0 
                h_2=h*h;h_3=h_2*h;h_4=h_3*h;h_5=h_4*h;h_6=h_5*h
                Δ22=4.0+h*coefΔh1+h_2*(coefΔh2)+h_3*(coefΔh3)+h_4*(coefΔh4)+h_5*coefΔh5+h_6*coefΔh6
                if abs(Δ22)!=0.0
                  qi=(4.0*xi+h*coefH1+h_2*coefH2+h_3*coefH3+h_4*coefH4+h_5*coefH5+h_6*coefH6)/Δ22
                  qj=(4.0*xj+h*coefH1j+h_2*coefH2j+h_3*coefH3j+h_4*coefH4j+h_5*coefH5j+h_6*coefH6j)/Δ22
                end
              end
          end
        end
 
        if  0<h<1e-20  #ie simul step failed # h==0 is allowed because is like an equilibrium
          println("simulstep ord2 failed small h=",h)
          return false
        end
        if Δ1==0.0 # should never happen when maxIter>0...before deployement, group the 2 prev if statments
          println("iter ord2 simul Δ1 ==0")
          return false
        end
        
        q[index][0]=qi# store back helper vars
        q[j][0]=qj     
        q1parti=aii*qi+aij*qj+uij+h*uij2
        q1partj=aji*qi+ajj*qj+uji+h*uji2
        q[index][1]=((1-h*ajj)/Δ1)*q1parti+(h*aij/Δ1)*q1partj# store back helper vars
        q[j][1]=(h*aji/Δ1)*q1parti+((1-h*aii)/Δ1)*q1partj
        trackSimul[1]+=1 # do not have to recomputeNext if qi never changed
      end #end if iscycle
  return iscycle
end

