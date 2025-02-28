
"""
    nmisCycle_and_simulUpdate(aij::Float64,aji::Float64,trackSimul,::Val{1},Val(M),index::Int,j::Int,dirI::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exactA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64,f::F) where {M,F}

Performs a simultaneous update of two quantized variables `qi` and `qj` if cycle conditions are met.

# Arguments
- `aij::Float64`: linear approximation Coefficient for the jacobian entry between variable i and variable j.
- `aji::Float64`: linear approximation Coefficient for the jacobian entry between variable j and variable i.
- `trackSimul`: A tracking object for the simulataneous update.
- `::Val{1}`: A type parameter indicating the method order is order 1.
- `::Val{M}`: A type parameter indicating the detection mechanism.
- `index::Int`: The index of the current variable.
- `j::Int`: The index of the interacting variable.
- `dirI::Float64`: Direction of variable i.
- `x::Vector{Taylor0}`: State vector of Taylor series coefficients for the variables.
- `q::Vector{Taylor0}`: Quantized state vector of Taylor series coefficients for the variables.
- `quantum::Vector{Float64}`: Quantum levels for the variables.
- `exactA::Function`: Function to compute the exact value of a jacobian entry.
- `d::Vector{Float64}`: discrete variables.
- `cacheA::MVector{1,Float64}`: Cache for jacobian entry computation.
- `dxaux::Vector{MVector{1,Float64}}`: Auxiliary vector for old x values.
- `qaux::Vector{MVector{1,Float64}}`: Auxiliary vector for old quantized values.
- `tx::Vector{Float64}`: Time vector for state updates.
- `tq::Vector{Float64}`: Time vector for quantized state updates.
- `simt::Float64`: Current simulation time.
- `ft::Float64`: Final time for the simulation.

# Returns
- None. The function performs in-place updates on the quantized state vectors.
"""
function nmisCycle_and_simulUpdate(aij::Float64,aji::Float64,trackSimul,::Val{1},::Val{M},index::Int,j::Int,dirI::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exactA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64,f::F) where {M,F}
  cacheA[1]=0.0;exactA(q,d,cacheA,index,index, simt,f)
  aii=cacheA[1]
  cacheA[1]=0.0;exactA(q,d,cacheA,j,j, simt,f)
  ajj=cacheA[1]
 xi = x[index][0]
  xj = x[j][0]
  ẋi = x[index][1]
  ẋj = x[j][1]
  qi = q[index][0]
  qj = q[j][0]
  quanj = quantum[j]
  quani = quantum[index]
  xj = x[j][0]
  qiminus = qaux[index][1]
  uji = ẋj - ajj * qj - aji * qiminus
  uij = dxaux[index][1] - aii * qiminus - aij * qj
  iscycle = false
  dxj = aji * qi + ajj * qj + uji #only future qi   #emulate fj
  dxP = aii * qi + aij * qj + uij #only future qi                                          
  qjplus = xj + sign(dxj) * quanj  #emulate updateQ(j)...
  dxi = aii * qi + aij * qjplus + uij #both future qi & qj   #emulate fi

  iscycle=detect1(Val(M),ẋi,dxP,dxi,ẋj,dxj)


  if iscycle 
      h = ft-simt
      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
      qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
      qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
      if (abs(qi - xi) > 2.0*quani || abs(qj - xj) > 2.0*quanj) 
        h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
        qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
      end
      maxIter=2
      while (abs(qi - xi) >2.0* quani || abs(qj - xj) >2.0*quanj) && (maxIter>0)
        maxIter-=1
        if maxIter < 1
          return false
        end
        h1 = h * sqrt(quani / abs(qi - xi));
        h2 = h * sqrt(quanj / abs(qj - xj));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
        qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
       
      end 
      #trackSimul[1]+=1
      q[index][0]=qi# store back helper vars
      q[j][0]=qj
      tq[j]=simt 
 end # end outer dependency check
 return iscycle
end 


#tripleUpdates
#= function nmisCycle_and_simulUpdate(aij::Float64,aji::Float64,trackSimul,A,I,U,X,SD,jac,taylorOpsCache,t,CS,ff,clF,nextStateTime,::Val{1},::Val{10},index::Int,j::Int,dirI::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exactA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64,f::F) where {F}
  cacheA[1]=0.0;exactA(q,d,cacheA,index,index, simt,f)
  aii=cacheA[1]
  cacheA[1]=0.0;exactA(q,d,cacheA,j,j, simt,f)
  ajj=cacheA[1]
  xi = x[index][0]
  xj = x[j][0]
  ẋi = x[index][1]
  ẋj = x[j][1]
  qi = q[index][0]
  qj = q[j][0]
  quanj = quantum[j]
  quani = quantum[index]
  xj = x[j][0]
  qiminus = qaux[index][1]
  uji = ẋj - ajj * qj - aji * qiminus
  uij = dxaux[index][1] - aii * qiminus - aij * qj
  iscycle = false
  dxj = aji * qi + ajj * qj + uji #only future qi   #emulate fj
  dxP = aii * qi + aij * qj + uij #only future qi                                          
  qjplus = xj + sign(dxj) * quanj  #emulate updateQ(j)...
  dxi = aii * qi + aij * qjplus + uij #both future qi & qj   #emulate fi

  iscycle=detect1(Val(1),ẋi,dxP,dxi,ẋj,dxj)

  istriplecycle = false
  if iscycle 
        for k in SD(index)
          if trackSimul[1] != k && j!= k
            for b in (jac(k))    # update Qb: to be used to calculate exacte Akb
              elapsedq = simt - tq[b]
              if elapsedq > 0
                integrateState(Val(0), q[b], elapsedq)
                tq[b] = simt
              end
            end
            cacheA[1] = 0.0
            exactA(q, d, cacheA, index, k, simt,f)
            aik = cacheA[1]# 
            cacheA[1] = 0.0
            exactA(q, d, cacheA, k, index, simt,f)
            aki = cacheA[1]
            if k != index && aik * aki != 0.0
              cacheA[1]=0.0;exactA(q,d,cacheA,k,k, simt,f)
              akk=cacheA[1]
              xk = x[k][0]
              ẋk = x[k][1]
              qk = q[k][0]
              quank = quantum[k]
              xk = x[k][0]
              uki = ẋk - akk * qk - aki * qiminus
              uik = dxaux[index][1] - aii * qiminus - aik * qk
              
              dxk = aki * qi + akk * qk + uki #only future qi   #emulate fk
              dxP = aii * qi + aik * qk + uik #only future qi                                          
              qkplus = xk + sign(dxk) * quank  #emulate updateQ(k)...
              dxi = aii * qi + aik * qkplus + uik #both future qi & qk   #emulate fi
         
              istriplecycle=detect1(Val(1),ẋi,dxP,dxi,ẋk,dxk)
              if istriplecycle
                exactA(q, d, cacheA, j, k, simt,f)
                ajk = cacheA[1]# 
                cacheA[1] = 0.0
                exactA(q, d, cacheA, k, j, simt,f)
                akj = cacheA[1]
                A[1,1]=aii;A[1,2]=aij;A[1,3]=aki
                A[2,1]=aji;A[2,2]=ajj;A[2,3]=ajk
                A[3,1]=aki;A[3,2]=akj;A[3,3]=akk
                h = ft-simt
                X[1]=xi;X[2]=xj;X[3]=xk
                uijk=uij-aik*qk
                ujik=uji-ajk*qk
                ukij=uki-akj*qj
                U[1]=uijk;U[2]=ujik;U[3]=ukij
                Q=inv(I-h*A)*(X+h*U)
                qi=Q[1];qj=Q[2];qk=Q[3]
                if (abs(qi - xi) > 2.0*quani || abs(qj - xj) > 2.0*quanj || abs(qk - xk) > 2.0*quank) 
                  h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));h3 = (abs(quank / ẋk));
                  h=min(h1,h2,h3)
                  Q=inv(I-h*A)*(X+h*U)
                  qi=Q[1];qj=Q[2];qk=Q[3]
                end
                maxIter=30
                while (abs(qi - xi) >2.0* quani || abs(qj - xj) >2.0*quanj || abs(qk - xk) >2.0*quank) && (maxIter>0)
                  maxIter-=1
                  if maxIter < 1
                    istriplecycle=false
                  end
                  h1 = h * sqrt(quani / abs(qi - xi));
                  h2 = h * sqrt(quanj / abs(qj - xj));
                  h3 = h * sqrt(quank / abs(qk - xk));
                  h=min(h1,h2,h3)
                 # if h>1e-2
                    Q=inv(I-h*A)*(X+h*U)
                    qi=Q[1];qj=Q[2];qk=Q[3]
                 # end
         
                end
                if h<1e-12
                  istriplecycle=false
                end
                if istriplecycle
                  q[index][0]=qi# store back helper vars
                  q[j][0]=qj
                  q[k][0]=qk
                  tq[j]=simt
                  tq[k]=simt
                  trackSimul[1] = k
                    for c in SD(k)  #k influences c
                      if c != index && c != k && c != j
                        elapsedx = simt - tx[c]
                        x[c].coeffs[1] = x[c](elapsedx)
                        tx[c] = simt
                        elapsedq = simt - tq[c]
                        if elapsedq > 0 integrateState(Val(0), q[c], elapsedq) ;tq[c] = simt end
                        for b in (jac(c))
                          elapsedq = simt - tq[b]
                          if elapsedq > 0
                            integrateState(Val(0), q[b], elapsedq)
                            tq[b] = simt
                          end
                        end
                        clearCache(taylorOpsCache, Val(CS), Val(1))
                        ff(c, q, t, d,taylorOpsCache,clF)
                        computeDerivative(Val(1), x[c], taylorOpsCache[1])
                        Liqss_reComputeNextTime(Val(1), c, simt, nextStateTime, x, q, quantum)
                      end#end if c!=0
                    end#end for c depend on k      







                  
                    #break # only one triple update is allowed for now
                    return true
                  end
              end# end if istriplecycle
            end # end if k != index && aik * aki != 0.0
          end# end if trackSimul[1] != k
        end# end for k in SD(index)

        if !istriplecycle
        
              h = ft-simt
              Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
              qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
              qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
              if (abs(qi - xi) > 2.0*quani || abs(qj - xj) > 2.0*quanj) 
                h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));
                h=min(h1,h2)
                Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                if Δ==0
                  Δ=1e-12
                end
                qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
                qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
              end
              maxIter=2
              while (abs(qi - xi) >2.0* quani || abs(qj - xj) >2.0*quanj) && (maxIter>0)
                maxIter-=1
                if maxIter < 1
                  return false
                end
                h1 = h * sqrt(quani / abs(qi - xi));
                h2 = h * sqrt(quanj / abs(qj - xj));
                h=min(h1,h2)
                Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                if Δ==0
                  Δ=1e-12
                end
                qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
                qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
              
              end 
              q[index][0]=qi# store back helper vars
              q[j][0]=qj
              tq[j]=simt 
        end
      #trackSimul[1]+=1

 end # end outer dependency check
 return iscycle
end  =#

function detect1(::Val{0},ẋi,dxP,dxi,ẋj,dxj)
 
  return false

end
function detect1(::Val{1},ẋi,dxP,dxi,ẋj,dxj)
if (dxj*ẋj)<0.0 && (dxi*ẋi)<0.0
    return true
else  
  return false
end
end

function detect1(::Val{2},ẋi,dxP,dxi,ẋj,dxj)
if (dxj*ẋj)<0.0 && (dxi*dxP)<0.0
    return true
else  
  #@show dxj,ẋj,dxi,dxP
  return false
end
end

function detect1(::Val{3},ẋi,dxP,dxi,ẋj,dxj)
if abs(dxj-ẋj)>abs(dxj+ẋj)/2 && abs(dxi-ẋi)>abs(dxi+ẋi)/2
 
    return true
 
else  
  return false
end
end

function detect1(::Val{4},ẋi,dxP,dxi,ẋj,dxj)
if abs(dxj-ẋj)>abs(dxj+ẋj)/2 && abs(dxi-dxP)>abs(dxi+dxP)/2
 
    return true
 
else  
  return false
end
end

function detect1(::Val{5},ẋi,dxP,dxi,ẋj,dxj)
  if ((abs(dxj)>3*abs(ẋj)||3*abs(dxj)<abs(ẋj)) && dxi*ẋi<0.0) || ((abs(dxi)>3*abs(ẋi)||3*abs(dxi)<abs(ẋi)) && dxj*ẋj<0.0) || ((dxj*ẋj)<0.0 && (dxi*ẋi)<0.0)
   
      return true
   
  else  
    return false
  end
end
function detect1(::Val{6},ẋi,dxP,dxi,ẋj,dxj)
  if ((abs(dxj)>3*abs(ẋj)||3*abs(dxj)<abs(ẋj)) && dxi*dxP<0.0) || ((abs(dxi)>3*abs(dxP)||3*abs(dxi)<abs(dxP)) && dxj*ẋj<0.0) || ((dxj*ẋj)<0.0 && (dxi*dxP)<0.0)
    
      return true
    
  else  
    return false
  end
end
#cases to remove: one value null and other very small

function detect1(::Val{7},ẋi,dxP,dxi,ẋj,dxj)
  if (dxj*ẋj)<=0.0 && (dxi*ẋi)<=0.0
    if (dxj)==0.0 && abs(ẋj)<1e-9
      return false
    end
    if (ẋj)==0.0 && abs(dxj)<1e-9
      return false
    end
    if (dxi)==0.0 && abs(ẋi)<1e-9
      return false
    end
    if (ẋi)==0.0 && abs(dxi)<1e-9
      return false
    end
      return true
    
  else  
    return false
  end
  end
function detect1(::Val{8},ẋi,dxP,dxi,ẋj,dxj)
  if (dxj*ẋj)<=0.0 && (dxi*dxP)<=0.0
    if (dxj)==0.0 && abs(ẋj)<1e-9
      return false
    end
    if (ẋj)==0.0 && abs(dxj)<1e-9
      return false
    end
    if (dxi)==0.0 && abs(dxP)<1e-9
      return false
    end
    if (dxP)==0.0 && abs(dxi)<1e-9
      return false
    end
      return true
    
  else  
    #@show dxj,ẋj,dxi,dxP
    return false
  end
end
  
function detect1(::Val{9},ẋi,dxP,dxi,ẋj,dxj)
  if abs(dxj-ẋj)>abs(dxj+ẋj)/2 && abs(dxi-ẋi)>abs(dxi+ẋi)/2
    if (dxj)==0.0 && abs(ẋj)<1e-9
      return false
    end
    if (ẋj)==0.0 && abs(dxj)<1e-9
      return false
    end
    if (dxi)==0.0 && abs(ẋi)<1e-9
      return false
    end
    if (ẋi)==0.0 && abs(dxi)<1e-9
      return false
    end
      return true
    
  else  
    return false
  end
end

function detect1(::Val{10},ẋi,dxP,dxi,ẋj,dxj)
  if abs(dxj-ẋj)>abs(dxj+ẋj)/2 && abs(dxi-dxP)>abs(dxi+dxP)/2
    if (dxj)==0.0 && abs(ẋj)<1e-9
      return false
    end
    if (ẋj)==0.0 && abs(dxj)<1e-9
      return false
    end
    if (dxi)==0.0 && abs(dxP)<1e-9
      return false
    end
    if (dxP)==0.0 && abs(dxi)<1e-9
      return false
    end
      return true
    
  else  
    return false
  end
end

function detect1(::Val{11},ẋi,dxP,dxi,ẋj,dxj)
  if (dxj*ẋj)<=0.0 && (dxi*ẋi)<=0.0
    
      return true
    
  else  
    return false
  end
end

function detect1(::Val{12},ẋi,dxP,dxi,ẋj,dxj)
  if  (dxi*dxP)<0.0
    
      return true
    
  else  
    return false
  end
end

function detect1(::Val{13},ẋi,dxP,dxi,ẋj,dxj)
    if  (dxi*dxP)<=0.0
      
        return true
      
    else  
      return false
    end
end

        

#= function detect1(::Val{14},ẋi,dxP,dxi,ẋj,dxj)
  iscycle=false 
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
  return iscycle

end =#