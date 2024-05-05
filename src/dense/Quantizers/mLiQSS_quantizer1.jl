#analy
function nmisCycle_and_simulUpdate(cacheRootsi::Vector{Float64},cacheRootsj::Vector{Float64},acceptedi::Vector{Vector{Float64}},acceptedj::Vector{Vector{Float64}},aij::Float64,aji::Float64,respp::Ptr{Float64}, pp::Ptr{NTuple{2,Float64}},trackSimul,::Val{1},index::Int,j::Int,dirI::Float64,dti::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
 
  cacheA[1]=0.0;exacteA(q,d,cacheA,index,index,simt)
  aii=cacheA[1]
  cacheA[1]=0.0;exacteA(q,d,cacheA,j,j,simt)
  ajj=cacheA[1]
 #=  exacteA(q,cacheA,index,j)
  aij=cacheA[1]
  exacteA(q,cacheA,j,index)
  aji=cacheA[1] =#


  xi=x[index][0];xj=x[j][0];ẋi=x[index][1];ẋj=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j];quani=quantum[index]
  #qaux[j][1]=qj;
  elapsed = simt - tx[j];x[j][0]= xj+elapsed*ẋj;

  xj=x[j][0]
  tx[j]=simt

 qiminus=qaux[index][1]
   #ujj=ẋj-ajj*qj
    uji=ẋj-ajj*qj-aji*qiminus
    #uii=dxaux[index][1]-aii*qaux[index][1]
    uij=dxaux[index][1]-aii*qiminus-aij*qj
    iscycle=false
    dxj=aji*qi+ajj*qj+uji #only future qi   #emulate fj
   
    

    dxithrow=aii*qi+aij*qj+uij #only future qi
                                                      
  qjplus=xj+sign(dxj)*quanj  #emulate updateQ(j)...

    dxi=aii*qi+aij*qjplus+uij #both future qi & qj   #emulate fi
    #dxj2=ajj*qjplus+aji*qi+uji
    
  
   #=  if abs(dxithrow)<1e-15 && dxithrow!=0.0
      dxithrow=0.0
    end =#
                             
   ########condition:Union 
 #= if abs(dxj)*3<abs(ẋj) || abs(dxj)>3*abs(ẋj) || (dxj*ẋj)<0.0 
    if abs(dxi)>3*abs(ẋi) || abs(dxi)*3<abs(ẋi) ||  (dxi*ẋi)<0.0 
        iscycle=true
    end
  end   =#                           
    ########condition:Union i
  #=   if abs(dxj-ẋj)>(abs(dxj+ẋj)/2)  
      if abs(dxi-dxithrow)>(abs(dxi+dxithrow)/2) 
        iscycle=true
      end
    end =#
    ########condition:Union i union
    if (abs(dxj)*3<abs(ẋj) || abs(dxj)>3*abs(ẋj) || (dxj*ẋj)<0.0)
      cancelCriteria=1e-6*quani
      if abs(dxi)>3*abs(dxithrow) || abs(dxi)*3<abs(dxithrow) ||  (dxi*dxithrow)<0.0  
          iscycle=true
          if abs(dxi)<cancelCriteria && abs(dxithrow)<cancelCriteria #significant change is relative to der values. cancel if all values are small
            iscycle=false
           # println("cancel dxi & dxthr small")
          end
      end
      if  abs(dxi)>10*abs(ẋi) || abs(dxi)*10<abs(ẋi) #||  (dxi*ẋi)<0.0 
        iscycle=true
        if abs(dxi)<cancelCriteria && abs(ẋi)<cancelCriteria #significant change is relative to der values. cancel if all values are small
          iscycle=false
         # println("cancel dxi & ẋi small")
        end
      end
      
      if abs(dxj)<1e-6*quanj && abs(ẋj)<1e-6*quanj #significant change is relative to der values. cancel if all values are small
        iscycle=false
        #println("cancel dxj & ẋj small")
      end

    end   
   

    #=  if 9.087757893221387e-5<=simt<=9.23684654009e-5
      println("simulUpdate $iscycle at simt=$simt index=$index, j=$j")
      @show  ẋi,dxi,dxithrow,ẋj,dxj
      @show qj,xj
     
      
     end  =#
  
 
  if iscycle
    #clear accIntrvals and cache of roots
    for i =1:3# 3 ord1 ,7 ord2
      acceptedi[i][1]=0.0; acceptedi[i][2]=0.0
      acceptedj[i][1]=0.0; acceptedj[i][2]=0.0
    end 
    for i=1:4
      cacheRootsi[i]=0.0
      cacheRootsj[i]=0.0
    end
   
      
       #find positive zeros f=+-Δ
      bi=aii*xi+aij*xj+uij;ci=aij*(aji*xi+uji)-ajj*(aii*xi+uij);αi=-ajj-aii;βi=aii*ajj-aij*aji
      bj=ajj*xj+aji*xi+uji;cj=aji*(aij*xj+uij)-aii*(ajj*xj+uji);αj=-aii-ajj;βj=ajj*aii-aji*aij
    
      coefi=NTuple{3,Float64}((βi*quani-ci,αi*quani-bi,quani))
      coefi2=NTuple{3,Float64}((-βi*quani-ci,-αi*quani-bi,-quani))
      coefj=NTuple{3,Float64}((βj*quanj-cj,αj*quanj-bj,quanj))
      coefj2=NTuple{3,Float64}((-βj*quanj-cj,-αj*quanj-bj,-quanj))

   


      resi1,resi2=quadRootv2(coefi)
      resi3,resi4=quadRootv2(coefi2)
      resj1,resj2=quadRootv2(coefj)
      resj3,resj4=quadRootv2(coefj2)
     #=  unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2)
      allrealrootintervalnewtonregulafalsi(coefi,respp,pp)
      resi1,resi2=unsafe_load(respp,1),unsafe_load(respp,2) 
   
      unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2)
      allrealrootintervalnewtonregulafalsi(coefi2,respp,pp)
      resi3,resi4=unsafe_load(respp,1),unsafe_load(respp,2) 
  
      unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2)
      allrealrootintervalnewtonregulafalsi(coefj,respp,pp)
      resj1,resj2=unsafe_load(respp,1),unsafe_load(respp,2) 
     
      unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2)
      allrealrootintervalnewtonregulafalsi(coefj2,respp,pp)
      resj3,resj4=unsafe_load(respp,1),unsafe_load(respp,2)  =#
      



      #construct intervals
      constructIntrval(acceptedi,resi1,resi2,resi3,resi4)

      constructIntrval(acceptedj,resj1,resj2,resj3,resj4)
  
      #find best H (largest overlap)
      ki=3;kj=3
      while true
            currentLi=acceptedi[ki][1];currentHi=acceptedi[ki][2]
            currentLj=acceptedj[kj][1];currentHj=acceptedj[kj][2]
            if currentLj<=currentLi<currentHj || currentLi<=currentLj<currentHi#resj[end][1]<=resi[end][1]<=resj[end][2] || resi[end][1]<=resj[end][1]<=resi[end][2] #overlap
                        h=min(currentHi,currentHj )
                        #= if simt==6.167911302732222e-9
                          @show h
                        end =#
                        if h==Inf && (currentLj!=0.0 || currentLi!=0.0) # except case both  0 sols  h1=(0,Inf)
                          h=max(currentLj,currentLi#= ,ft-simt,dti =#) #ft-simt in case they are both ver small...elaborate on his later
                        end
                        if h==Inf # both lower bounds ==0  --> zero sols for both
                          qi=getQfromAsymptote(simt,xi,βi,ci,αi,bi)
                          qj=getQfromAsymptote(simt,xj,βj,cj,αj,bj)
                         
                      
                        else
                            Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                            if Δ!=0.0
                              qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
                              qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
                            else
                              qi,qj,h,maxIter=IterationH(h,xi,quani, xj,quanj,aii,ajj,aij,aji,uij,uji)
                              if maxIter==0
                                return false
                              end
                            end
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
          
      q[j][0]=qj
      tq[j]=simt 
      q[index][0]=qi# store back helper vars
      trackSimul[1]+=1 # do not have to recomputeNext if qi never changed
     
  
  end
  return iscycle
end   
 



#iters
# function nmisCycle_and_simulUpdate(acceptedi::Vector{Vector{Float64}},acceptedj::Vector{Vector{Float64}},aij::Float64,aji::Float64,respp::Ptr{Float64}, pp::Ptr{NTuple{2,Float64}},trackSimul,::Val{1},index::Int,j::Int,dirI::Float64,dti::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
#= function nmisCycle_and_simulUpdate(cacheRootsi::Vector{Float64},cacheRootsj::Vector{Float64},acceptedi::Vector{Vector{Float64}},acceptedj::Vector{Vector{Float64}},aij::Float64,aji::Float64,respp::Ptr{Float64}, pp::Ptr{NTuple{2,Float64}},trackSimul,::Val{1},index::Int,j::Int,dirI::Float64,dti::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
 
  cacheA[1]=0.0;exacteA(q,d,cacheA,index,index)
  aii=cacheA[1]
  cacheA[1]=0.0;exacteA(q,d,cacheA,j,j)
  ajj=cacheA[1]
 #=  exacteA(q,cacheA,index,j)
  aij=cacheA[1]
  exacteA(q,cacheA,j,index)
  aji=cacheA[1] =#


  xi=x[index][0];xj=x[j][0];ẋi=x[index][1];ẋj=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j];quani=quantum[index]
  #qaux[j][1]=qj;
  elapsed = simt - tx[j];x[j][0]= xj+elapsed*ẋj;

  xj=x[j][0]
  tx[j]=simt

 qiminus=qaux[index][1]
   #ujj=ẋj-ajj*qj
    uji=ẋj-ajj*qj-aji*qiminus
    #uii=dxaux[index][1]-aii*qaux[index][1]
    uij=dxaux[index][1]-aii*qiminus-aij*qj
    iscycle=false
    dxj=aji*qi+ajj*qj+uji #only future qi   #emulate fj
   
    

    dxithrow=aii*qi+aij*qj+uij #only future qi
                                                      
  qjplus=xj+sign(dxj)*quanj  #emulate updateQ(j)...

    dxi=aii*qi+aij*qjplus+uij #both future qi & qj   #emulate fi
    #dxj2=ajj*qjplus+aji*qi+uji
    
  
   #=  if abs(dxithrow)<1e-15 && dxithrow!=0.0
      dxithrow=0.0
    end =#
                             
   ########condition:Union 
 #= if abs(dxj)*3<abs(ẋj) || abs(dxj)>3*abs(ẋj) || (dxj*ẋj)<0.0 
    if abs(dxi)>3*abs(ẋi) || abs(dxi)*3<abs(ẋi) ||  (dxi*ẋi)<0.0 
        iscycle=true
    end
  end   =#                           
    ########condition:Union i
  #=   if abs(dxj-ẋj)>(abs(dxj+ẋj)/2)  
      if abs(dxi-dxithrow)>(abs(dxi+dxithrow)/2) 
        iscycle=true
      end
    end =#
    ########condition:Union i union
    if (abs(dxj)*3<abs(ẋj) || abs(dxj)>3*abs(ẋj) || (dxj*ẋj)<0.0)
      cancelCriteria=1e-3*quani
      if abs(dxi)>3*abs(dxithrow) || abs(dxi)*3<abs(dxithrow) ||  (dxi*dxithrow)<0.0  
          iscycle=true
          if abs(dxi)<cancelCriteria && abs(dxithrow)<cancelCriteria #significant change is relative to der values. cancel if all values are small
            iscycle=false
           # println("cancel dxi & dxthr small")
          end
      end
      if  abs(dxi)>10*abs(ẋi) || abs(dxi)*10<abs(ẋi) #||  (dxi*ẋi)<0.0 
        iscycle=true
        if abs(dxi)<cancelCriteria && abs(ẋi)<cancelCriteria #significant change is relative to der values. cancel if all values are small
          iscycle=false
         # println("cancel dxi & ẋi small")
        end
      end
      
      if abs(dxj)<1e-6*quanj && abs(ẋj)<1e-6*quanj #significant change is relative to der values. cancel if all values are small
        iscycle=false
        #println("cancel dxj & ẋj small")
      end

    end   
   

    #=  if 9.087757893221387e-5<=simt<=9.23684654009e-5
      println("simulUpdate $iscycle at simt=$simt index=$index, j=$j")
      @show  ẋi,dxi,dxithrow,ẋj,dxj
      @show qj,xj
     
      
     end  =#
  
 
   if iscycle
   
    h=ft-simt
    # h=firstguessH
     Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
     qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
     qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
   if (abs(qi - xi) > 1*quani || abs(qj - xj) > 1*quanj) #checking qi-xi is not needed since firstguess just made it less than delta
     h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));
     h=min(h1,h2)
     h_two=h
     Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
     #= if Δ==0
       Δ=1e-12
     end =#
     qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
     qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
   end
   maxIter=1000
   while (abs(qi - xi) > 1*quani || abs(qj - xj) > 1*quanj) && (maxIter>0)
       maxIter-=1
       h1 = h * (0.99*quani / abs(qi - xi));
      #=  Δtemp=(1-h1*aii)*(1-h1*ajj)-h1*h1*aij*aji
       qitemp = ((1-h1*ajj)*(xi+h1*uij)+h1*aij*(xj+h1*uji))/Δtemp =#
       h2 = h * (0.99*quanj / abs(qj - xj));
       #= Δtemp=(1-h2*aii)*(1-h2*ajj)-h2*h2*aij*aji
       qjtemp = ((1-h2*ajj)*(xi+h2*uij)+h2*aij*(xj+h2*uji))/Δtemp
       if abs(qitemp - xi) > 1*quani || abs(qjtemp - xj) > 1*quanj
         println("regle de croix did not work")
       end =#
       h=min(h1,h2)
       h_three=h
       Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
       #= if Δ==0
         Δ=1e-12
         println("delta liqss1 simulupdate==0")
       end =#
       if Δ!=0
       qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
       qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
       end
      
       if maxIter < 1 println("maxiter of updateQ      = ",maxIter) end
   end 
   if maxIter < 1  
    return false
   end
  
   q[index][0]=qi# store back helper vars
   trackSimul[1]+=1 # do not have to recomputeNext if qi never changed
   q[j][0]=qj
   #= push!(simulqxiVals,abs(qi-xi))
   push!(simulqxjVals,abs(qj-xj))
   push!(simuldeltaiVals,quani)
   push!(simuldeltajVals,quanj) =#
   tq[j]=simt 
 end #end second dependecy check

  return iscycle
end   =# 




#= function nmisCycle_and_simulUpdate(acceptedi::Vector{Vector{Float64}},acceptedj::Vector{Vector{Float64}},aij::Float64,aji::Float64,respp::Ptr{Float64}, pp::Ptr{NTuple{2,Float64}},trackSimul,::Val{1},index::Int,j::Int,dirI::Float64,dti::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
  exacteA(q,d,cacheA,index,index)
  aii=cacheA[1]
  exacteA(q,d,cacheA,j,j)
  ajj=cacheA[1]
  exacteA(q,d,cacheA,index,j)
  aij=cacheA[1]
  exacteA(q,d,cacheA,j,index)
  aji=cacheA[1]

  xi=x[index][0];xj=x[j][0];ẋi=x[index][1];ẋj=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j];quani=quantum[index]
  qaux[j][1]=qj;#olddx[j][1]=ẋj
     
  elapsed = simt - tx[j];x[j][0]= xj+elapsed*ẋj;
  xj=x[j][0]
  tx[j]=simt

 
   #ujj=ẋj-ajj*qj
    uji=ẋj-ajj*qj-aji*qaux[index][1]
    #uii=dxaux[index][1]-aii*qaux[index][1]
    uij=dxaux[index][1]-aii*qaux[index][1]-aij*qj#qaux[j][1]
    iscycle=false
    dxj=aji*qi+ajj*qaux[j][1]+uji #only future qi
    qjplus=xj+sign(dxj)*quanj

      dxi=aii*qi+aij*qjplus+uij #both future qi & qj
      dxi2=aii*qi+aij*qj+uij #only future qi
      dxj2=aji*qi+ajj*qjplus+uji#both future qi & qj
 
        
   #=   if simt>=0.22816661756287676
    @show simt,index
    @show ẋi,dxi
    @show ẋj,dxj
    @show xj,qj,qjplus,quanj
    @show xi,qaux[index][1],qi,quani
    @show firstguessH
   end =#
     if  #= abs(dxj-ẋj)>abs(dxj+ẋj)/1.8 =#((abs(dxj)*3<abs(ẋj) || abs(dxj)>3*abs(ẋj))#= && abs(ẋj * dxj) > 1e-3  =# )|| (dxj*ẋj)<=0.0#= && abs(dxj-ẋj)>abs(dxj+ẋj)/20 =#
    #  if (dxj*ẋj)<=0.0
      # if abs(a[j][index]*2*quantum[index])>abs(x[j][1])
    ############################################################################
   
     if #= abs(dxi-ẋi)>abs(dxi+ẋi)/1.8  =#((abs(dxi)*3<abs(ẋi) || abs(dxi)>3*abs(ẋi))#= && abs(ẋi * dxi)> 1e-3 =#  )|| (dxi*ẋi)<=0.0 #= && abs(dxi-ẋi)>abs(dxi+ẋi)/20 =#
    #if (dxi*ẋi)<=0.0 
    iscycle=true
   
      
        h_two=-1.0;h_three=-1.0
      
     
          h=ft-simt
         # h=firstguessH
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
          qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
        if (abs(qi - xi) > 1*quani || abs(qj - xj) > 1*quanj) #checking qi-xi is not needed since firstguess just made it less than delta
          h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));
          h=min(h1,h2)
          h_two=h
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          if Δ==0
            Δ=1e-12
          end
          qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
          qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
        end
        maxIter=1000
        while (abs(qi - xi) > 1*quani || abs(qj - xj) > 1*quanj) && (maxIter>0)
            maxIter-=1
            h1 = h * (0.99*quani / abs(qi - xi));
           #=  Δtemp=(1-h1*aii)*(1-h1*ajj)-h1*h1*aij*aji
            qitemp = ((1-h1*ajj)*(xi+h1*uij)+h1*aij*(xj+h1*uji))/Δtemp =#
            h2 = h * (0.99*quanj / abs(qj - xj));
            #= Δtemp=(1-h2*aii)*(1-h2*ajj)-h2*h2*aij*aji
            qjtemp = ((1-h2*ajj)*(xi+h2*uij)+h2*aij*(xj+h2*uji))/Δtemp
            if abs(qitemp - xi) > 1*quani || abs(qjtemp - xj) > 1*quanj
              println("regle de croix did not work")
            end =#
            h=min(h1,h2)
            h_three=h
            Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
            if Δ==0
              Δ=1e-12
              println("delta liqss1 simulupdate==0")
            end
            qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
            qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
           
            if maxIter < 1 println("maxiter of updateQ      = ",maxIter) end
        end 
       # if maxIter < 950 println("maxiter      = ",maxIter) end
       # @show simt,index,h
         #= push!(simulHTimes,simt)
       push!(simulHVals,h) =#
    #   @show h,qi,qj
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
        #= push!(simulqxiVals,abs(qi-xi))
        push!(simulqxjVals,abs(qj-xj))
        push!(simuldeltaiVals,quani)
        push!(simuldeltajVals,quanj) =#
        tq[j]=simt 
      end #end second dependecy check
   end # end outer dependency check
   return iscycle
end   =# 


function getQfromAsymptote(simt,x::Float64,β::P,c::P,α::P,b::P) where {P<:Union{BigFloat,Float64}}
  q=0.0
  if β==0.0 && c!=0.0
    println("report bug: β==0 && c!=0.0: this leads to asym=Inf while this code came from asym<Delta at $simt")
  elseif β==0.0 && c==0.0
    if α==0.0 && b!=0.0
      println("report bug: α==0 && b!=0.0 this leads to asym=Inf while this code came from asym<Delta at $simt")
    elseif α==0.0 && b==0.0
      q=x
    else
      q=x+b/α
    end
  else
       q=x+c/β #asymptote...even though not guaranteed to be best
  end
q
end

function iterationH(h::Float64,xi::Float64,quani::Float64, xj::Float64,quanj::Float64,aii::Float64,ajj::Float64,aij::Float64,aji::Float64,uij::Float64,uji::Float64)
  qi,qj=0.0,0.0
  maxIter=1000
  while (abs(qi - xi) > 1*quani || abs(qj - xj) > 1*quanj) && (maxIter>0)
      maxIter-=1
      h1 = h * (0.99*quani / abs(qi - xi));
      h2 = h * (0.99*quanj / abs(qj - xj));
      h=min(h1,h2)
      h_three=h
      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
       if Δ!=0
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
        qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
      end
     
     # if maxIter < 1 println("maxiter of updateQ      = ",maxIter) end
  end 
 
  qi,qj,h,maxIter
end