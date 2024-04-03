
 



#iters


function nmisCycle_and_simulUpdate(cacheRootsi::Vector{Float64},cacheRootsj::Vector{Float64},acceptedi::Vector{Vector{Float64}},acceptedj::Vector{Vector{Float64}},aij::Float64,aji::Float64,respp::Ptr{Float64}, pp::Ptr{NTuple{2,Float64}},trackSimul,::Val{1},index::Int,j::Int,dirI::Float64,dti::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exactA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
 
  cacheA[1]=0.0;exactA(q,d,cacheA,index,index,simt)
  aii=cacheA[1]
  cacheA[1]=0.0;exactA(q,d,cacheA,j,j,simt)
  ajj=cacheA[1]
 #=  exactA(q,cacheA,index,j)
  aij=cacheA[1]
  exactA(q,cacheA,j,index)
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
   
    h=ft-simt
    # h=firstguessH
     Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
     qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
     qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
   if (abs(qi - xi) > 2*quani || abs(qj - xj) > 2*quanj) #checking qi-xi is not needed since firstguess just made it less than delta
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
   while (abs(qi - xi) > 2*quani || abs(qj - xj) > 2*quanj) && (maxIter>0)
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
end    


