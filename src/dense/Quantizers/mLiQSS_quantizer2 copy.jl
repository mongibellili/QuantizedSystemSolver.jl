

function nmisCycle_and_simulUpdate(cacheRootsi::Vector{Float64},cacheRootsj::Vector{Float64},acceptedi::Vector{Vector{Float64}},acceptedj::Vector{Vector{Float64}},aij::Float64,aji::Float64,respp::Ptr{Float64}, pp::Ptr{NTuple{2,Float64}},trackSimul,::Val{2},index::Int,j::Int,dirI::Float64,firstguessH::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{2,Float64}},qaux::Vector{MVector{2,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
  

  cacheA[1]=0.0;exacteA(q,d,cacheA,index,index)
   aii=cacheA[1]
   cacheA[1]=0.0; exacteA(q,d,cacheA,j,j)
   ajj=cacheA[1]


   uii=dxaux[index][1]-aii*qaux[index][1]
   ui2=dxaux[index][2]-aii*qaux[index][2]

  #=  uii=x[index][1]-aii*q[index][0]
   ui2=x[index][2]-aii*q[index][1] =#

  xi=x[index][0];xj=x[j][0];qi=q[index][0];qj=q[j][0];qi1=q[index][1];qj1=q[j][1];xi1=x[index][1];xi2=2*x[index][2];xj1=x[j][1];xj2=2*x[j][2]
  #uii=u[index][index][1];ujj=u[j][j][1]#;uij=u[index][j][1];uji=u[j][index][1]#;uji2=u[j][index][2]
  quanj=quantum[j];quani=quantum[index];
  e1 = simt - tx[j];e2 = simt - tq[j];#= e3=simt - tu[j];tu[j]=simt;  =#
  prevXj=x[j][0]
  x[j][0]= x[j](e1);xj=x[j][0];tx[j]=simt

  qj=qj+e2*qj1  ;qaux[j][1]=qj;tq[j] = simt    ;q[j][0]=qj  

  xj1=x[j][1]+e1*xj2;
  newDiff=(xj-prevXj)

  ujj=xj1-ajj*qj
  #u[j][j][1]=ujj
  uji=ujj-aji*qaux[index][1]# 
  #uji=u[j][index][1]
  uj2=xj2-ajj*qj1###################################################-----------------------
  uji2=uj2-aji*qaux[index][2]#
 #u[j][index][2]=u[j][j][2]-ajj*qaux[index][1] # from article p20 line25 more cycles ...shaky with no bumps
 # uji2=u[j][index][2] 
  dxj=aji*qi+ajj*qaux[j][1]+uji
  ddxj=aji*qi1+ajj*qj1+uji2
  if abs(ddxj)==0.0
    ddxj=1e-30
    if DEBUG2  @show ddxj end
  end
  iscycle=false


    qjplus=xj-sign(ddxj)*quanj
    hj=sqrt(2*quanj/abs(ddxj))#
    α1=1-hj*ajj
    if abs(α1)==0.0
      α1=1e-30
      @show α1
    end
    dqjplus=(aji*(qi+hj*qi1)+ajj*qjplus+uji+hj*uji2)/α1
 

    uij=uii-aij*qaux[j][1]
    # uij=u[index][j][1]
     uij2=ui2-aij*qj1#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
       #β=dxi+sqrt(abs(ddxi)*quani/2)
      #h2=sqrt(2*quani/abs(ddxi))
    
      #uij2=u[index][j][2]
      dxithrow=aii*qi+aij*qj+uij
      ddxithrow=aii*qi1+aij*qj1+uij2
      dxi=aii*qi+aij*qjplus+uij
      ddxi=aii*qi1+aij*dqjplus+uij2
      if abs(ddxi)==0.0
        ddxi=1e-30
        if DEBUG2 @show ddxi end
      end
      hi=sqrt(2*quani/abs(ddxi))
      βidir=dxi+hi*ddxi/2
      βjdir=dxj+hj*ddxj/2
      βidth=dxithrow+hi*ddxithrow/2
      αidir=xi1+hi*xi2/2
      βj=xj1+hj*xj2/2
########condition:Union 
  #=   if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  || dqjplus*newDiff<0.0 #(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
    
      if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)) || βidir*dirI<0.0
        iscycle=true
      end
    end =#
 ########condition:Union i
  #=   if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  || dqjplus*newDiff<0.0 #(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
    
      if (abs(dxi-dxithrow)>(abs(dxi+dxithrow)/2) || abs(ddxi-ddxithrow)>(abs(ddxi+ddxithrow)/2)) || βidir*dirI<0.0
        iscycle=true
      end
    end =#
 ########condition:combineDer Union i


 if abs(βjdir-βj)>(abs(βjdir+βj)/2)  || dqjplus*newDiff<0.0 
    
          if abs(βidir-βidth)>(abs(βidir+βidth)/2)  || βidir*dirI<0.0
            iscycle=true
          end
        coef_signig=100
          if (abs(dxi)>(abs(xi1)*coef_signig) #= && (abs(ddxi)>(abs(xi2)*coef_signig)||abs(xi2)>(abs(ddxi)*coef_signig)) =#) || (abs(xi1)>(abs(dxi)*coef_signig) #= &&  (abs(ddxi)>(abs(xi2)*coef_signig)||abs(xi2)>(abs(ddxi)*coef_signig)) =#)
            iscycle=true
          end
        #=  if abs(βidir-αidir)>(abs(βidir+αidir)/2)  
            iscycle=true
          end =#
   end


if (ddxithrow>0 && dxithrow<0 && qi1>0) || (ddxithrow<0 && dxithrow>0 && qi1<0)
  xmin=xi-dxithrow*dxithrow/ddxithrow
  if abs(xmin-xi)/abs(qi-xi)>0.8 # xmin is close to q and a change in q  might flip its sign of dq
    iscycle=false #
  end
end


#= if (ddxi>0 && dxi<0 && qi1>0) || (ddxi<0 && dxi>0 && qi1<0)
  xmin=xi-dxi*dxi/ddxi
  if abs(xmin-xi)/abs(qi-xi)>0.8 # xmin is close to q and a change in q  might flip its sign of dq
  
    iscycle=false #
  end

end =#

########condition:kinda signif alone i
  #=   if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2)) # || dqjplus*newDiff<0.0 
    
      if (abs(dxi-dxithrow)>(abs(dxi+dxithrow)/2) || abs(ddxi-ddxithrow)>(abs(ddxi+ddxithrow)/2)) #|| βidir*dirI<0.0
        iscycle=true
      end
    end =#

     ########condition:cond1 
      #=    if  dqjplus*newDiff<0.0 #= || (dqjplus==0.0 && newDiff!=0.0)  =##(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
    
                if βidir*dirI<0.0
                  iscycle=true
                end
              end

              if  dqjplus*newDiff<0.0 #= || (dqjplus==0.0 && newDiff!=0.0)  =##(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
              
                if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)) 
                  iscycle=true
                end
              end

              if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))   #= || (dqjplus==0.0 && newDiff!=0.0)  =##(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
              
                if  βidir*dirI<0.0
                  iscycle=true
                end
    end =#
 


     ########condition:cond1 i
     #=     if  dqjplus*newDiff<0.0 #= || (dqjplus==0.0 && newDiff!=0.0)  =##(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
    
              if βidir*dirI<0.0
                iscycle=true
              end
            end

            if  dqjplus*newDiff<0.0 #= || (dqjplus==0.0 && newDiff!=0.0)  =##(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
            
              if (abs(dxi-dxithrow)>(abs(dxi+dxithrow)/2) || abs(ddxi-ddxithrow)>(abs(ddxi+ddxithrow)/2)) 
                iscycle=true
              end
            end

            if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))   #= || (dqjplus==0.0 && newDiff!=0.0)  =##(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
            
              if  βidir*dirI<0.0
                iscycle=true
              end
    end =#

       #=  if index==1 && j==4 && 0.1<simt<1.1
      iscycle=true
    end =#



  if iscycle
   #=  quani=1.7*quani
    quanj=1.7*quanj =#
        for i =1:7# 3 ord1 ,7 ord2
          acceptedi[i][1]=0.0; acceptedi[i][2]=0.0
          acceptedj[i][1]=0.0; acceptedj[i][2]=0.0
        end 
        for i=1:12
          cacheRootsi[i]=0.0
          cacheRootsj[i]=0.0
        end

          trackSimul[1]+=1 
          aiijj=(aii+ajj)
          aiijj_=(aij*aji-aii*ajj)



          coefΔh0=4.0
          coefΔh1=-8.0*aiijj
          coefΔh2=6.0*aiijj*aiijj-4.0*aiijj_
          coefΔh3=6.0*aiijj*aiijj_-2.0*aiijj*aiijj*aiijj
          coefΔh4=aiijj_*aiijj_-4.0*aiijj_*aiijj*aiijj
          coefΔh5=-3.0*aiijj*aiijj_*aiijj_
          coefΔh6=-aiijj_*aiijj_*aiijj_
          
        
        
          coefH0=4.0*xi
          coefH1=-8.0*xi*aiijj
          #coefH2=2.0*(-aii*uij-aij*uji-uij2+xi*(-2.0*aiijj_+2.0*aiijj*aiijj+2.0*aii*ajj-aji*aij+ajj*ajj)-aij*aiijj*xj)
          coefH2=2.0*(-aii*uij-aij*uji-uij2+xi*(-3.0*aiijj_+2.0*aiijj*aiijj+aiijj*ajj)-aij*aiijj*xj)
          #coefH3=2.0*((-uij*aiijj_+ajj*uij2-aij*uji2)+aiijj*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)-xi*aiijj*(2.0*aii*ajj-aji*aij+ajj*ajj)+ajj*aiijj_*xi+aij*aiijj*aiijj*xj-aij*aiijj_*xj)
          coefH3=2.0*((-uij*aiijj_+ajj*uij2-aij*uji2)+aiijj*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)-xi*aiijj*(ajj*aiijj-aiijj_)+ajj*aiijj_*xi+aij*aiijj*aiijj*xj-aij*aiijj_*xj)
          coefH4=2.0*(aiijj*(uij*aiijj_-ajj*uij2+aij*uji2)-0.5*(2.0*aii*ajj-aji*aij+ajj*ajj)*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)-ajj*aiijj_*aiijj*xi+0.5*aij*aiijj*(ajj*uji+aji*uij+uji2+2.0*xj*aiijj_)+aij*aiijj_*aiijj*xj)
          coefH5=(2.0*aii*ajj-aji*aij+ajj*ajj)*(ajj*uij2-uij*aiijj_-aij*uji2)-ajj*aiijj_*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)-aij*aiijj*(aii*uji2-uji*aiijj_-aji*uij2)+aij*aiijj_*(ajj*uji+aji*uij+uji2+2.0*xj*aiijj_)
          coefH6=ajj*aiijj_*(ajj*uij2-uij*aiijj_-aij*uji2)-aij*aiijj_*(aii*uji2-uji*aiijj_-aji*uij2)
          
          coefH0j=4.0*xj
          coefH1j=-8.0*xj*aiijj
          coefH2j=2.0*(-ajj*uji-aji*uij-uji2+xj*(-2.0*aiijj_+2.0*aiijj*aiijj+2.0*ajj*aii-aij*aji+aii*aii)-aji*aiijj*xi)
          coefH3j=2.0*((-uji*aiijj_+aii*uji2-aji*uij2)+aiijj*(ajj*uji+aji*uij+uji2+2.0*xj*aiijj_)-xj*aiijj*(2.0*ajj*aii-aij*aji+aii*aii)+aii*aiijj_*xj+aji*aiijj*aiijj*xi-aji*aiijj_*xi)
          coefH4j=2.0*(aiijj*(uji*aiijj_-aii*uji2+aji*uij2)-0.5*(aii*aiijj-aiijj_)*(ajj*uji+aji*uij+uji2+2.0*xj*aiijj_)-aii*aiijj_*aiijj*xj+0.5*aji*aiijj*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)+aji*aiijj_*aiijj*xi)
          coefH5j=(aii*aiijj-aiijj_)*(aii*uji2-uji*aiijj_-aji*uij2)-aii*aiijj_*(ajj*uji+aji*uij+uji2+2.0*xj*aiijj_)-aji*aiijj*(ajj*uij2-uij*aiijj_-aij*uji2)+aji*aiijj_*(aii*uij+aij*uji+uij2+2.0*xi*aiijj_)
          coefH6j=aii*aiijj_*(aii*uji2-uji*aiijj_-aji*uij2)-aji*aiijj_*(ajj*uij2-uij*aiijj_-aij*uji2)


        #=   coefi=NTuple{7,Float64}((coefH6-coefΔh6*(xi+quani),coefH5-coefΔh5*(xi+quani),coefH4-coefΔh4*(xi+quani),coefH3-coefΔh3*(xi+quani),coefH2-coefΔh2*(xi+quani),coefH1-coefΔh1*(xi+quani),coefH0-coefΔh0*(xi+quani)))
          coefi2=NTuple{7,Float64}((coefH6-coefΔh6*(xi-quani),coefH5-coefΔh5*(xi-quani),coefH4-coefΔh4*(xi-quani),coefH3-coefΔh3*(xi-quani),coefH2-coefΔh2*(xi-quani),coefH1-coefΔh1*(xi-quani),coefH0-coefΔh0*(xi-quani)))

          coefj=NTuple{7,Float64}((coefH6j-coefΔh6*(xj+quanj),coefH5j-coefΔh5*(xj+quanj),coefH4j-coefΔh4*(xj+quanj),coefH3j-coefΔh3*(xj+quanj),coefH2j-coefΔh2*(xj+quanj),coefH1j-coefΔh1*(xj+quanj),coefH0j-coefΔh0*(xj+quanj)))
          coefj2=NTuple{7,Float64}((coefH6j-coefΔh6*(xj-quanj),coefH5j-coefΔh5*(xj-quanj),coefH4j-coefΔh4*(xj-quanj),coefH3j-coefΔh3*(xj-quanj),coefH2j-coefΔh2*(xj-quanj),coefH1j-coefΔh1*(xj-quanj),coefH0j-coefΔh0*(xj-quanj)))
          resi1,resi2,resi3,resi4,resi5,resi6,resi7,resi8,resi9,resi10,resi11,resi12=0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
          resj1,resj2,resj3,resj4,resj5,resj6,resj7,resj8,resj9,resj10,resj11,resj12=0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
          Base.GC.enable(false)

          unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2);unsafe_store!(respp, -1.0, 3);unsafe_store!(respp, -1.0, 4);unsafe_store!(respp, -1.0, 5);unsafe_store!(respp, -1.0, 6)
          allrealrootintervalnewtonregulafalsi(coefi,respp,pp)
          resi1,resi2,resi3,resi4,resi5,resi6=unsafe_load(respp,1),unsafe_load(respp,2),unsafe_load(respp,3),unsafe_load(respp,4),unsafe_load(respp,5),unsafe_load(respp,6) 
      
          unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2);unsafe_store!(respp, -1.0, 3);unsafe_store!(respp, -1.0, 4);unsafe_store!(respp, -1.0, 5);unsafe_store!(respp, -1.0, 6)
          allrealrootintervalnewtonregulafalsi(coefi2,respp,pp)
          resi7,resi8,resi9,resi10,resi11,resi12=unsafe_load(respp,1),unsafe_load(respp,2),unsafe_load(respp,3),unsafe_load(respp,4),unsafe_load(respp,5),unsafe_load(respp,6) 
      
          unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2);unsafe_store!(respp, -1.0, 3);unsafe_store!(respp, -1.0, 4);unsafe_store!(respp, -1.0, 5);unsafe_store!(respp, -1.0, 6)
          allrealrootintervalnewtonregulafalsi(coefj,respp,pp)
          resj1,resj2,resj3,resj4,resj5,resj6=unsafe_load(respp,1),unsafe_load(respp,2),unsafe_load(respp,3),unsafe_load(respp,4),unsafe_load(respp,5),unsafe_load(respp,6) 
          
          unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2);unsafe_store!(respp, -1.0, 3);unsafe_store!(respp, -1.0, 4);unsafe_store!(respp, -1.0, 5);unsafe_store!(respp, -1.0, 6)
          allrealrootintervalnewtonregulafalsi(coefj2,respp,pp)
          resj7,resj8,resj9,resj10,resj11,resj12=unsafe_load(respp,1),unsafe_load(respp,2),unsafe_load(respp,3),unsafe_load(respp,4),unsafe_load(respp,5),unsafe_load(respp,6) 
          
          Base.GC.enable(true) =#
            #to be deleted: test
     #=        for h in (resi1,resi2,resi3,resi4,resi5,resi6,resi7,resi8,resi9,resi10,resi11,resi12)
              if h>0
                h_2=h*h;h_3=h_2*h;h_4=h_3*h;h_5=h_4*h;h_6=h_5*h
                Δ22=4.0+h*coefΔh1+h_2*(coefΔh2)+h_3*(coefΔh3)+h_4*(coefΔh4)+h_5*coefΔh5+h_6*coefΔh6
                if abs(Δ22)!=0.0
                  qi=(4.0*xi+h*coefH1+h_2*coefH2+h_3*coefH3+h_4*coefH4+h_5*coefH5+h_6*coefH6)/Δ22
                  if abs(qi-xi)>2quani
                    error("your own h not working for you, i")
                  end
                end
              end
            end
            for h in (resj1,resj2,resj3,resj4,resj5,resj6,resj7,resj8,resj9,resj10,resj11,resj12)
              if h>0
                h_2=h*h;h_3=h_2*h;h_4=h_3*h;h_5=h_4*h;h_6=h_5*h
                Δ22=4.0+h*coefΔh1+h_2*(coefΔh2)+h_3*(coefΔh3)+h_4*(coefΔh4)+h_5*coefΔh5+h_6*coefΔh6
                if abs(Δ22)!=0.0
                  qj=(4.0*xj+h*coefH1j+h_2*coefH2j+h_3*coefH3j+h_4*coefH4j+h_5*coefH5j+h_6*coefH6j)/Δ22
                  if abs(qj-xj)>2quanj
                    error("your own h not working for you, j")
                  end
              end
              end
            end

            #filter roots
            counterRooti=1
            for h in (resi1,resi2,resi3,resi4,resi5,resi6,resi7,resi8,resi9,resi10,resi11,resi12)
              if h>0 && h!=Inf
                cacheRootsi[counterRooti]=h
                counterRooti+=1
              end
            end
            sort!(cacheRootsi);

            counterRootj=1
            for h in (resj1,resj2,resj3,resj4,resj5,resj6,resj7,resj8,resj9,resj10,resj11,resj12)
              if h>0 && h!=Inf
                cacheRootsj[counterRootj]=h
                counterRootj+=1
              end
            end =#

            coefi=[coefH0-coefΔh0*(xi+quani),coefH1-coefΔh1*(xi+quani),coefH2-coefΔh2*(xi+quani),coefH3-coefΔh3*(xi+quani),coefH4-coefΔh4*(xi+quani),coefH5-coefΔh5*(xi+quani),coefH6-coefΔh6*(xi+quani)]
            coefj=[coefH0j-coefΔh0*(xj+quanj),coefH1j-coefΔh1*(xj+quanj),coefH2j-coefΔh2*(xj+quanj),coefH3j-coefΔh3*(xj+quanj),coefH4j-coefΔh4*(xj+quanj),coefH5j-coefΔh5*(xj+quanj),coefH6j-coefΔh6*(xj+quanj)]
            coefi2=[coefH0-coefΔh0*(xi-quani),coefH1-coefΔh1*(xi-quani),coefH2-coefΔh2*(xi-quani),coefH3-coefΔh3*(xi-quani),coefH4-coefΔh4*(xi-quani),coefH5-coefΔh5*(xi-quani),coefH6-coefΔh6*(xi-quani)]
            coefj2=[coefH0j-coefΔh0*(xj-quanj),coefH1j-coefΔh1*(xj-quanj),coefH2j-coefΔh2*(xj-quanj),coefH3j-coefΔh3*(xj-quanj),coefH4j-coefΔh4*(xj-quanj),coefH5j-coefΔh5*(xj-quanj),coefH6j-coefΔh6*(xj-quanj)]
          
         #=   coefi=[BigFloat(coefH0-coefΔh0*(xi+quani)),BigFloat(coefH1-coefΔh1*(xi+quani)),BigFloat(coefH2-coefΔh2*(xi+quani)),BigFloat(coefH3-coefΔh3*(xi+quani)),BigFloat(coefH4-coefΔh4*(xi+quani)),BigFloat(coefH5-coefΔh5*(xi+quani)),BigFloat(coefH6-coefΔh6*(xi+quani))]
           coefj=[BigFloat(coefH0j-coefΔh0*(xj+quanj)),BigFloat(coefH1j-coefΔh1*(xj+quanj)),BigFloat(coefH2j-coefΔh2*(xj+quanj)),BigFloat(coefH3j-coefΔh3*(xj+quanj)),BigFloat(coefH4j-coefΔh4*(xj+quanj)),BigFloat(coefH5j-coefΔh5*(xj+quanj)),BigFloat(coefH6j-coefΔh6*(xj+quanj))]
          coefi2=[BigFloat(coefH0-coefΔh0*(xi-quani)),BigFloat(coefH1-coefΔh1*(xi-quani)),BigFloat(coefH2-coefΔh2*(xi-quani)),BigFloat(coefH3-coefΔh3*(xi-quani)),BigFloat(coefH4-coefΔh4*(xi-quani)),BigFloat(coefH5-coefΔh5*(xi-quani)),BigFloat(coefH6-coefΔh6*(xi-quani))]
          coefj2=[BigFloat(coefH0j-coefΔh0*(xj-quanj)),BigFloat(coefH1j-coefΔh1*(xj-quanj)),BigFloat(coefH2j-coefΔh2*(xj-quanj)),BigFloat(coefH3j-coefΔh3*(xj-quanj)),BigFloat(coefH4j-coefΔh4*(xj-quanj)),BigFloat(coefH5j-coefΔh5*(xj-quanj)),BigFloat(coefH6j-coefΔh6*(xj-quanj))]
          =#
          ri=roots(coefi)
          rj=roots(coefj)
          ri2=roots(coefi2)
          rj2=roots(coefj2)
          
          
            #filter roots
            counterRooti=1
            for h in ri
              if abs(imag(h)) < 1.0e-15 && real(h)>0.0
                cacheRootsi[counterRooti]=Float64(real(h))
                counterRooti+=1
              end
            end
            checkRoots(cacheRootsi,coefi,coefΔh0,coefΔh1,coefΔh2,coefΔh3,coefΔh4,coefΔh5,coefΔh6,quani)        
           
            for h in ri2
              if abs(imag(h)) < 1.0e-15 && real(h)>0.0
                cacheRootsi[counterRooti]=Float64(real(h))
                counterRooti+=1
              end
            end
            sort!(cacheRootsi);
  
            counterRootj=1
            for h in rj
              if abs(imag(h)) < 1.0e-15 && real(h)>0.0
                cacheRootsj[counterRootj]=Float64(real(h))
                counterRootj+=1
              end
            end
                 
            checkRoots(cacheRootsj,coefj,coefΔh0,coefΔh1,coefΔh2,coefΔh3,coefΔh4,coefΔh5,coefΔh6,quanj)
            for h in rj2
              if abs(imag(h)) < 1.0e-15 && real(h)>0.0
                cacheRootsj[counterRootj]=Float64(real(h))
                counterRootj+=1
              end
            end
            sort!(cacheRootsj);

      #=     for h in cacheRootsi
            
              h_2=h*h;h_3=h_2*h;h_4=h_3*h;h_5=h_4*h;h_6=h_5*h
              Δ22=4.0+h*coefΔh1+h_2*(coefΔh2)+h_3*(coefΔh3)+h_4*(coefΔh4)+h_5*coefΔh5+h_6*coefΔh6
              if abs(Δ22)!=0.0
                qi=(4.0*xi+h*coefH1+h_2*coefH2+h_3*coefH3+h_4*coefH4+h_5*coefH5+h_6*coefH6)/Δ22
                if abs(qi-xi)>2quani
                  error("your own h not working for you, i")
                end
              end
            
          end
          for h in cacheRootsj
            
              h_2=h*h;h_3=h_2*h;h_4=h_3*h;h_5=h_4*h;h_6=h_5*h
              Δ22=4.0+h*coefΔh1+h_2*(coefΔh2)+h_3*(coefΔh3)+h_4*(coefΔh4)+h_5*coefΔh5+h_6*coefΔh6
              if abs(Δ22)!=0.0
                qj=(4.0*xj+h*coefH1j+h_2*coefH2j+h_3*coefH3j+h_4*coefH4j+h_5*coefH5j+h_6*coefH6j)/Δ22
                if abs(qj-xj)>2quanj
                  error("your own h not working for you, j")
                end
            
              end
          end =#

         
          #construct intervals
          constructIntrval(acceptedi,cacheRootsi)
          constructIntrval(acceptedj,cacheRootsj)
          ki=7;kj=7
          h=-1.0

          while true
                if ki==0 || kj==0
                  error("bug in finding the overlap")
                end
                currentLi=acceptedi[ki][1];currentHi=acceptedi[ki][2]
                currentLj=acceptedj[kj][1];currentHj=acceptedj[kj][2]
                if currentLj<=currentLi<currentHj || currentLi<=currentLj<currentHi#resj[end][1]<=resi[end][1]<=resj[end][2] || resi[end][1]<=resj[end][1]<=resi[end][2] #overlap
                            h=min(currentHi,currentHj )
                            if h==Inf && (currentLj!=0.0 || currentLi!=0.0) # except case both  0 sols  h1=(0,Inf)
                              h=max(currentLj,currentLi#= ,ft-simt,dti =#) #ft-simt in case they are both ver small...elaborate on his later
                            #=   if simt==2.500745083181994e-5
                                @show currentLj,currentLi
                              end =#
                            end
                            if h==Inf # both lower bounds ==0  --> zero sols for both
                                  if aiijj_!=0.0
                                    qi=coefH6/coefΔh6
                                    qj=coefH6j/coefΔh6
                                    #qi=getQfromAsymptote(xi,coefH0,coefH1,coefH2,coefH3,coefH4,coefH5,coefH6,coefΔh0,coefΔh1,coefΔh2,coefΔh3,coefΔh4,coefΔh5,coefΔh6)
                                    #qj=getQfromAsymptote(xj,coefH0j,coefH1j,coefH2j,coefH3j,coefH4j,coefH5j,coefH6j,coefΔh0,coefΔh1,coefΔh2,coefΔh3,coefΔh4,coefΔh5,coefΔh6)
                                    qi1=(aji*uij2-uji2*aij)/aiijj_ #resul of ddx=aq+aq+u=0 ...2eqs 2 unknowns
                                    qj1=(-uij2-aii*qi1)/aij
                                  else
                                    return false
                                  end
                            else #h!=Inf
                                  h_2=h*h;h_3=h_2*h;h_4=h_3*h;h_5=h_4*h;h_6=h_5*h
                                  Δ22=4.0+h*coefΔh1+h_2*(coefΔh2)+h_3*(coefΔh3)+h_4*(coefΔh4)+h_5*coefΔh5+h_6*coefΔh6
                                  Δ1=1.0-h*(aiijj)-h*h*(aiijj_)
                                  if abs(Δ22)!=0.0 && Δ1!=0.0
                                          qi=(4.0*xi+h*coefH1+h_2*coefH2+h_3*coefH3+h_4*coefH4+h_5*coefH5+h_6*coefH6)/Δ22
                                          qj=(4.0*xj+h*coefH1j+h_2*coefH2j+h_3*coefH3j+h_4*coefH4j+h_5*coefH5j+h_6*coefH6j)/Δ22
                                        #=  if 2.500745083181994e-5<=simt<3.00745083181994e-5
                                              @show "simul: ",index,qi,xi,quani,h,simt,j,qj,xj,quanj
                                          end =#
                                  else
                                          println("singularities Δ22=0")#do iterations
                                          maxIter=10000
                                          while (abs(qi - xi) > 2.0*quani || abs(qj - xj) > 2.0*quanj) && (maxIter>0)
                                                  maxIter-=1
                                                  h1 = h * sqrt(quani / abs(qi - xi));
                                                  h2 = h * sqrt(quanj / abs(qj - xj));
                                                  #=  h1 = h * (0.99*1.8*quani / abs(qi - xi));
                                                  h2 = h * (0.99*1.8*quanj / abs(qj - xj)); =#
                                                  h=min(h1,h2)
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
                                          if maxIter==0 || Δ1==0.0 # maxIter=0 :failed to bring q down within a reasonable time. Δ1==0.0 should never happen when maxIter>0...but keep it there for now
                                            return false
                                          end
                                  end
                                  q1parti=aii*qi+aij*qj+uij+h*uij2
                                  q1partj=aji*qi+ajj*qj+uji+h*uji2
                                  qi1=((1-h*ajj)/Δ1)*q1parti+(h*aij/Δ1)*q1partj# Δ1!=0 from iterations
                                  qj1=(h*aji/Δ1)*q1parti+((1-h*aii)/Δ1)*q1partj
                            end
                            break
                else
                          if currentHj==0.0 && currentHi==0.0#empty
                            ki-=1;kj-=1
                          elseif currentHj==0.0 && currentHi!=0.0
                            kj-=1
                          elseif currentHi==0.0 && currentHj!=0.0
                            ki-=1
                          else #both non zero
                                if currentLj<currentLi#remove last seg of resi
                                  ki-=1
                                else #remove last seg of resj
                                  kj-=1
                                end
                          end
                end
    
          end # end while

         

        #checkRoots(cacheRootsi,coefi,coefi2) # for free stress-testing allrealrootintervalnewtonregulafalsi
       # checkRoots(cacheRootsj,coefj,coefj2)  #later check if error in roots makes q-x >2quan...use commented code below

        q[index][0]=(qi)# store back helper vars
        q[j][0]=(qj)
        q[index][1]=(qi1)
        q[j][1]=(qj1)
        trackSimul[1]+=1 # do not have to recomputeNext if qi never changed

          
  end #end if iscycle
    
    
  return iscycle
end






  
  
  
  
  












#= function checkRootsIntNew(cacheRootsi,coefi,coefi2)
    for h in cacheRootsi
      
          zci=abs(coefi[1]*h^6+coefi[2]*h^5+coefi[3]*h^4+coefi[4]*h^3+coefi[5]*h^2+coefi[6]*h+coefi[7])
          zci2=abs(coefi2[1]*h^6+coefi2[2]*h^5+coefi2[3]*h^4+coefi2[4]*h^3+coefi2[5]*h^2+coefi2[6]*h+coefi2[7])
          if zci>1e-6 && zci2>1e-6
            println("error in mliqss_quantizer2: coming from root solver ")
            @show h,coefi,zci
            @show coefi2,zci2
          end
     
    end
end =#

#= function checkRoots(cacheRootsi,coefi,coefi2)
  for h in cacheRootsi
      if h>0 # cache is large and contains 0 0 0 0 0 0
        zci=abs(coefi[7]*h^6+coefi[6]*h^5+coefi[5]*h^4+coefi[4]*h^3+coefi[3]*h^2+coefi[2]*h+coefi[1])
        zci2=abs(coefi2[7]*h^6+coefi2[6]*h^5+coefi2[5]*h^4+coefi2[4]*h^3+coefi2[3]*h^2+coefi2[2]*h+coefi2[1])
        if zci>1e-5 && zci2>1e-5
          println("error in mliqss_quantizer2: coming from root solver ")
          if zci>1e-5  @show h,coefi,zci end
          if zci2>1e-5  @show h,coefi2,zci2 end
        end
      end
   
  end
end =#
function checkRoots(cacheRootsi,coefi,coefΔh0,coefΔh1,coefΔh2,coefΔh3,coefΔh4,coefΔh5,coefΔh6,quani)
  for h in cacheRootsi
      if h>0 # cache is large and contains 0 0 0 0 0 0
        zci=abs(coefi[7]*h^6+coefi[6]*h^5+coefi[5]*h^4+coefi[4]*h^3+coefi[3]*h^2+coefi[2]*h+coefi[1])
        Δ=coefΔh0+h*coefΔh1+coefΔh2*h*h+coefΔh3*h^3+coefΔh4*h^4+coefΔh5*h^5+coefΔh6*h^6
        if (zci/Δ*Δ)>quani
          error("zci not zero makes q changes by more than a quan...root finding not good...Δ=$Δ")
        end

       #=  if zci>1e-5 
          println("error in mliqss_quantizer2: coming from root solver ")
           @show h,coefi,zci 
         
        end =#
      end
   
  end
end

#= function getQfromAsymptote(x::Float64,coefH0::Float64,coefH1::P,coefH2::P,coefH3::P,coefH4::P,coefH5::P,coefH6::P,coefΔh0::Float64,coefΔh1::P,coefΔh2::P,coefΔh3::P,coefΔh4::P,coefΔh5::P,coefΔh6::P) where {P<:Union{BigFloat,Float64}}
  q=0.0
  if coefΔh6==0.0 && coefH6!=0.0
    @show coefΔh6,coefH6
    error("report bug in simul2: could appear when a var becomes a NaN")
   
  elseif coefΔh6==0.0 && coefH6==0.0
        if coefΔh5==0.0 && coefH5!=0.0
          println("simul2: coef5 asymptote should not be outside of +-delta !!!")
        elseif coefΔh5==0.0 && coefH5==0.0
              if coefΔh4==0.0 && coefH4!=0.0
                println("simul2: coef4 asymptote should not be outside of +-delta !!!")
              elseif coefΔh4==0.0 && coefH4==0.0
                  if coefΔh3==0.0 && coefH3!=0.0
                    println("simul2: coef3 asymptote should not be outside of +-delta !!!")
                  elseif coefΔh3==0.0 && coefH3==0.0
                    if coefΔh2==0.0 && coefH2!=0.0
                      println("simul2: coef2 asymptote should not be outside of +-delta !!!")
                    elseif coefΔh2==0.0 && coefH2==0.0
                      if coefΔh1==0.0 && coefH1!=0.0
                        println("simul2: coef1 asymptote should not be outside of +-delta !!!")
                      elseif coefΔh1==0.0 && coefH1==0.0
                        if coefΔh0!=0.0
                          q=coefH0/coefΔh0
                        else
                          q=x
                        end
                      else#coefΔh1!=0.0
                        q=coefH1/coefΔh1
                      end
                    else#coefΔh2!=0.0
                      q=coefH2/coefΔh2
                    end
                  else#coefΔh3!=0.0
                    q=coefH3/coefΔh3
                  end
              else#coefΔh4!=0.0
                q=coefH4/coefΔh4
              end
        else#coefΔh5!=0.0
          q=coefH5/coefΔh5
        end
  else#coefΔh6!=0.0
    q=coefH6/coefΔh6
  end
q
end =#