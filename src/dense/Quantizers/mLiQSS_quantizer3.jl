

function nmisCycle_and_simulUpdate(::Val{3},index::Int,j::Int,prevStepVal::Float64,direction::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},a::Vector{Vector{Float64}},u::Vector{Vector{MVector{O,Float64}}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},olddxSpec::Vector{MVector{O,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64,qminus::Vector{Float64})where{O}
   aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index]
  xi=x[index][0];xj=x[j][0];qi=q[index][0];qj=q[j][0];qi1=q[index][1];qi2=2*q[index][2];qj1=q[j][1];qj2=2*q[j][2]
  uii=u[index][index][1];uij=u[index][j][1];ujj=u[j][j][1]#;uji=u[j][index][1]#;uji2=u[j][index][2]
  quanj=quantum[j];quani=quantum[index];xi1=x[index][1];xi2=2*x[index][2];xi3=6*x[index][3];xj1=x[j][1];xj2=2*x[j][2];xj3=6*x[j][3]
  e1 = simt - tx[j];e2 = simt - tq[j]; tx[j]=simt;
  x[j][0]= x[j](e1);xjaux=x[j][0]

  xj1=xj1+e1*xj2+e1*e1*xj3/2;olddxSpec[j][1]=xj1;olddx[j][1]=xj1; xj2=xj2+e1*xj3
  qj=qj+e2*qj1+e2*e2*qj2/2  ;qj1=qj1+e2*qj2; qaux[j][1]=qj ;qaux[j][2]=qj1

  newDiff=(x[j][0]-prevStepVal)
  dir=direction[j]
  if newDiff*dir <0.0
    direction[j]=-dir

  elseif newDiff==0 && dir!=0.0
    direction[j]=0.0

  elseif newDiff!=0 && dir==0.0
    direction[j]=newDiff
  else
  #do not update direction
  end          

  ujj=ujj+e1*u[j][j][2]+e1*e1*u[j][j][3]/2  
  #ujj=xj1-ajj*qj
  u[j][j][1]=ujj
  #u[j][j][2]=u[j][j][2]+e1*u[j][j][3]
  u[j][j][2]=xj2-ajj*qj1
  u[j][j][3]=xj3-ajj*qj2###################################
  u[j][index][1]=ujj-aji*qaux[index][1]# 
  uji=u[j][index][1]
  u[j][index][2]=u[j][j][2]-aji*qaux[index][2]#less cycles but with a bump at 1.5...ft20: smooth with some bumps
  uji2=u[j][index][2]
   u[j][index][3]=u[j][j][3]-aji*qaux[index][3]#
   uji3=u[j][index][3]

  dxj=aji*qi+ajj*qaux[j][1]+uji
  ddxj=aji*qi1+ajj*qj1+uji2
  dddxj=aji*qi2+ajj*qj2+uji3
 
  iscycle=false

  u[index][j][1]=u[index][index][1]-aij*qaux[j][1]
  uij=u[index][j][1]
  u[index][j][2]=u[index][index][2]-aij*qj1
  uij2=u[index][j][2]
  u[index][j][3]=u[index][index][3]-aij*qj2
  uij3=u[index][j][3]

  # if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2) || abs(dddxj-xj3)>(abs(dddxj+xj3)/2))  
    
  h1=cbrt(abs(6*quanj/dddxj))
  qjplus=xjaux+h1*h1*h1*dddxj/6
  #@show h,qjplus,xjaux
  λ=ajj*qjplus+aji*qi+uji+h1*(aji*qi1+uji2)+h1*h1*(aji*qi2+uji3)/2
  α=(h1-h1*h1/2)*(aji*(qi1+h1*qi2)+uji2+h1*uji3)/(1-h1*ajj)
  β=1-h1*ajj+(h1-h1*h1/2)*ajj/(1-h1*ajj)
  dqjplus=(λ-α)/β
  ddqjplus=(aji*(qi1+h1*qi2)+ajj*dqjplus+uji2+h1*uji3)/(1-h1*ajj)
  ### 
  #u[index][j][1]=u[index][index][1]-a[index][j]*q[j][0]  # shifts down at 18
  
  #=  u[index][j][1]=u[index][index][1]-aij*qaux[j][1]
  uij=u[index][j][1]
  u[index][j][2]=u[index][index][2]-aij*qj1#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
  uij2=u[index][j][2] =#
  dxi=aii*qi+aij*qjplus+uij
  ddxi=aii*qi1+aij*dqjplus+uij2
  dddxi=aii*qi2+aij*ddqjplus+uij3
  
   #=  #dqjplus-qj1 is enough: normally (qjplus+h*dqjlpus-qj+h*qj1) . dqjplus-qj1 is better than dqjplus alone cuz the case dqjplus=0
  if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  || (dqjplus+h1*ddqjplus)*direction[j]<0.0 #(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
    h2=cbrt(abs(6*quani/dddxi))
    β=dxi+ddxi*h2/2+dddxi*h2*h2/6
    if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)) || β*direction[index]<0.0 =#

  if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2) || abs(dddxj-xj3)>(abs(dddxj+xj3)/2))  || (dqjplus+h1*ddqjplus)*direction[j]<0.0#(dqjplus+h1*ddqjplus-qj1-e1*qj2)*direction[j]<0.0
    h2=cbrt(abs(6*quani/dddxi))
    h2=cbrt(abs(6*quani/dddxi))
    β=dxi+ddxi*h2/2+dddxi*h2*h2/6
    if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)  || abs(dddxi-xi3)>(abs(dddxi+xi3)/2))|| β*direction[index]<0.0

   # if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)  || abs(dddxi-xi3)>(abs(dddxi+xi3)/2))

        iscycle=true
        h = ft-simt;
        qi,qj,Mi,Mj,Linvii,Linvij,Linvji,Linvjj,Nii,Nij,Nji,Njj=simulQ3(aii,aij,aji,ajj,h,xi,xjaux,uij,uij2,uij3,uji,uji2,uji3)
        if (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) # removing this did nothing...check @btime later
          h1 = cbrt(abs(6*quani/dddxi));h2 = cbrt(abs(6*quanj/dddxj));   #later add derderX =1e-12 when x2==0
          h=min(h1,h2)
          qi,qj,Mi,Mj,Linvii,Linvij,Linvji,Linvjj,Nii,Nij,Nji,Njj=simulQ3(aii,aij,aji,ajj,h,xi,xjaux,uij,uij2,uij3,uji,uji2,uji3)
        end  
        maxIter=600000
        while (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) && (maxIter>0) && (h!=0.0)
          maxIter-=1
         # h1 = h * (0.98*quani / abs(qi - xi));
          h1 = h *cbrt(quani / abs(qi - xi))
        #  h2 = h * (0.98*quanj / abs(qj - xjaux));
          h2 = h *cbrt(quanj / abs(qj - xjaux))
          h=min(h1,h2)
          qi,qj,Mi,Mj,Linvii,Linvij,Linvji,Linvjj,Nii,Nij,Nji,Njj=simulQ3(aii,aij,aji,ajj,h,xi,xjaux,uij,uij2,uij3,uji,uji2,uji3)  
        end                                                          
        if maxIter < 1
        println("maxtiter simult_val{3}")
        @show maxIter 
      end
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
      
        MAQi=Mi+aii*qi+aij*qj+uij+h*uij2+h*h*uij3/2;MAQj=Mj+ajj*qj+aji*qi+uji+h*uji2+h*h*uji3/2
        # MAQ=(M+A*Q+U+h*U2+(h*h/2)*U3)# Q1=Linv*MAQ
        q1i=Linvii*MAQi+Linvij*MAQj;q1j=Linvji*MAQi+Linvjj*MAQj
        q[index][1]=q1i#Q1[1]# store back helper vars
        q[j][1]=q1j#Q1[2]

        AQUi=aii*q1i+aij*q1j+uij2+h*uij3;AQUj=ajj*q1j+aji*q1i+uji2+h*uji3
        q2i=Nii*AQUi+Nij*AQUj;q2j=Nji*AQUi+Njj*AQUj
        q[index][2]=q2i/2#Q2[1]/2# store back helper vars: /2 for taylor standard storage
        q[j][2]=q2j/2#Q2[2]/2

        tq[j]=simt
      end #end second dependecy check
  end # end outer dependency check
  #= if debug  
    println("end of iscycle function")
   end =#
  return iscycle
 
end 



#######################################################################################################################################################



@inline function simulQ3(#= simul3Helper::MVector{10,Float64}, =#aii::Float64,aij::Float64,aji::Float64,ajj::Float64,h::Float64,xi::Float64,xjaux::Float64,uij::Float64,uij2::Float64,uij3::Float64,uji::Float64,uji2::Float64,uji3::Float64)
  #=
  the following matrix-based code is changed below (no matrices)
  N=inv(I-h*A)
  R=h*I-(h*h/2)*A+((h*h/2)*N*A-(h*h*h/6)*A*N*A)
  Linv=inv(I-h*A+h*N*A-(h*h/2)*A*N*A)
  M=((h*h/2)*A-h*I)*N*(U2+h*U3)
  P=X+(h*U+(h*h/2)*U2+(h*h*h/6)*U3)
  S=-R*Linv*(M+U+h*U2+(h*h/2)*U3)-(((h*h/2)*I-(h*h*h/6)*A)*N*(U2+h*U3))+P
  Q=inv(I-h*A+R*Linv*A)*S =# 
   h_2=h*h/2;h_6=h*h*h/6
  #  N=inv(I-h*A)
    Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
    Nii=(1-h*ajj)/Δ1;Nij=(h*aij)/Δ1;Nji=(aji*h)/Δ1;Njj=(1-h*aii)/Δ1#N=inv(I-h*A)
    NAii=Nii*aii+Nij*aji;NAij=Nii*aij+Nij*ajj;NAji=Nji*aii+Njj*aji;NAjj=Nji*aij+Njj*ajj
    ANAii=aii*NAii+aij*NAji;ANAij=aii*NAij+aij*NAjj;ANAji=aji*NAii+ajj*NAji;ANAjj=aji*NAij+ajj*NAjj
    Rii=h-(h_2)*aii+(h_2)*NAii-(h_6)*ANAii;Rij=-(h_2)*aij+(h_2)*NAij-(h_6)*ANAij;Rji=-(h_2)*aji+(h_2)*NAji-(h_6)*ANAji;Rjj=h-(h_2)*ajj+(h_2)*NAjj-(h_6)*ANAjj#R=h*I-(h*h/2)*A+((h*h/2)*N*A-(h*h*h/6)*A*N*A)
   # R=[Rii Rij;Rji Rjj]
    Lii=1-h*aii+h*NAii-(h_2)*ANAii;Lij=-h*aij+h*NAij-(h_2)*ANAij;Lji=-h*aji+h*NAji-(h_2)*ANAji;Ljj=1-h*ajj+h*NAjj-(h_2)*ANAjj
    ΔL=Lii*Ljj-Lij*Lji
    Linvii=Ljj/ΔL;Linvij=-Lij/ΔL;Linvji=-Lji/ΔL;Linvjj=Lii/ΔL## Linv=inv(I-h*A+h*N*A-(h*h/2)*A*N*A)
   # Linv=[Linvii Linvij;Linvji Linvjj]
   # M=((h*h/2)*A-h*I)*N*(U2+h*U3)
    U23i=uij2+h*uij3
    U23j=uji2+h*uji3
    NU23i=Nii*U23i+Nij*U23j
    NU23j=Nji*U23i+Njj*U23j
    Mi=(h_2*aii-h)*NU23i+h_2*aij*NU23j
    Mj=h_2*aji*NU23i+(h_2*ajj-h)*NU23j
   # M=[Mi;Mj]#M=((h*h/2)*A-h*I)*N*(U2+h*U3)
    Pi=xi+(h*uij+(h_2)*uij2+(h_6)*uij3);Pj=xjaux+(h*uji+(h_2)*uji2+(h_6)*uji3)# P=X+(h*U+(h*h/2)*U2+(h*h*h/6)*U3)
   # P=[Pi;Pj]
    S2i=(h_2-h_6*aii)*NU23i-h_6*aij*NU23j;S2j=(h_2-h_6*ajj)*NU23j-h_6*aji*NU23i# S2=(((h_2)*I-(h_6)*A)*N*(U2+h*U3))
   # S2=[S2i;S2j]
    S1i=(Mi+uij+h*uij2+(h_2)*uij3);S1j=(Mj+uji+h*uji2+(h_2)*uji3)# S1=(M+U+h*U2+(h_2)*U3)
    LinvS1i=Linvii*S1i+Linvij*S1j;LinvS1j=Linvji*S1i+Linvjj*S1j
   # LinvS1=[LinvS1i;LinvS1j]
    RLinvS1i=Rii*LinvS1i+Rij*LinvS1j;RLinvS1j=Rji*LinvS1i+Rjj*LinvS1j
   # RLinvS1=[RLinvS1i;RLinvS1j]
   # S=P-RLinvS1-S2
    Si=Pi-RLinvS1i-S2i
    Sj=Pj-RLinvS1j-S2j
    #S=[Si;Sj]
    LinvAii=Linvii*aii+Linvij*aji;LinvAij=Linvii*aij+Linvij*ajj;LinvAji=Linvji*aii+Linvjj*aji;LinvAjj=Linvji*aij+Linvjj*ajj
    RLinvAii=Rii*LinvAii+Rij*LinvAji;RLinvAij=Rii*LinvAij+Rij*LinvAjj;RLinvAji=Rji*LinvAii+Rjj*LinvAji;RLinvAjj=Rji*LinvAij+Rjj*LinvAjj
    #α=I_hARLinvA#
    αii=1-h*aii+RLinvAii; αij=-h*aij+RLinvAij;αji=-h*aji+RLinvAji; αjj=1-h*ajj+RLinvAjj
    Δα=αii*αjj-αij*αji
    InvAlphaii=αjj/Δα;InvAlphaij=-αij/Δα;InvAlphaji=-αji/Δα;InvAlphajj=αii/Δα;
   # InvAlpha=[InvAlphaii InvAlphaij;InvAlphaji InvAlphajj]
    qi=InvAlphaii*Si+InvAlphaij*Sj
    qj=InvAlphaji*Si+InvAlphajj*Sj
   return (qi,qj,Mi,Mj,Linvii,Linvij,Linvji,Linvjj,Nii,Nij,Nji,Njj)
end