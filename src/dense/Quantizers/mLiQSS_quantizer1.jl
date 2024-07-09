#analy
function nmisCycle_and_simulUpdate(cacheRootsi::Vector{Float64}, cacheRootsj::Vector{Float64}, acceptedi::Vector{Vector{Float64}}, acceptedj::Vector{Vector{Float64}}, aij::Float64, aji::Float64, respp::Ptr{Float64}, pp::Ptr{NTuple{2,Float64}}, trackSimul, ::Val{1}, index::Int, j::Int, dirI::Float64, dti::Float64, x::Vector{Taylor0}, q::Vector{Taylor0}, quantum::Vector{Float64}, exacteA::Function, d::Vector{Float64}, cacheA::MVector{1,Float64}, dxaux::Vector{MVector{1,Float64}}, qaux::Vector{MVector{1,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64)
  cacheA[1] = 0.0
  exacteA(q, d, cacheA, index, index, simt)
  aii = cacheA[1]
  cacheA[1] = 0.0
  exacteA(q, d, cacheA, j, j, simt)
  ajj = cacheA[1]
  xi = x[index][0]
  xj = x[j][0]
  ẋi = x[index][1]
  ẋj = x[j][1]
  qi = q[index][0]
  qj = q[j][0]
  quanj = quantum[j]
  quani = quantum[index]
  elapsed = simt - tx[j]
  x[j][0] = xj + elapsed * ẋj
  xj = x[j][0]
  tx[j] = simt
  qiminus = qaux[index][1]
  uji = ẋj - ajj * qj - aji * qiminus
  uij = dxaux[index][1] - aii * qiminus - aij * qj
  iscycle = false
  dxj = aji * qi + ajj * qj + uji #only future qi   #emulate fj
  dxithrow = aii * qi + aij * qj + uij #only future qi                                          
  qjplus = xj + sign(dxj) * quanj  #emulate updateQ(j)...
  dxi = aii * qi + aij * qjplus + uij #both future qi & qj   #emulate fi
  ########condition:Union 
  #= if abs(dxj)*3<abs(ẋj) || abs(dxj)>3*abs(ẋj) || (dxj*ẋj)<0.0 
     if abs(dxi)>3*abs(ẋi) || abs(dxi)*3<abs(ẋi) ||  (dxi*ẋi)<0.0 
         iscycle=true
     end
   end   =#
  ########condition:Union i union
  if (abs(dxj) * 3 < abs(ẋj) || abs(dxj) > 3 * abs(ẋj) || (dxj * ẋj) < 0.0)
    cancelCriteria = 1e-6 * quani
    if abs(dxi) > 3 * abs(dxithrow) || abs(dxi) * 3 < abs(dxithrow) || (dxi * dxithrow) < 0.0
      iscycle = true
      if abs(dxi) < cancelCriteria && abs(dxithrow) < cancelCriteria #significant change is relative to der values. cancel if all values are small
        iscycle = false
      end
    end
    if abs(dxi) > 10 * abs(ẋi) || abs(dxi) * 10 < abs(ẋi) #||  (dxi*ẋi)<0.0 
      iscycle = true
      if abs(dxi) < cancelCriteria && abs(ẋi) < cancelCriteria #significant change is relative to der values. cancel if all values are small
        iscycle = false
      end
    end
    if abs(dxj) < 1e-6 * quanj && abs(ẋj) < 1e-6 * quanj #significant change is relative to der values. cancel if all values are small
      iscycle = false
    end
  end
  if iscycle
    #clear accIntrvals and cache of roots
    for i = 1:3# 3 ord1 ,7 ord2
      acceptedi[i][1] = 0.0
      acceptedi[i][2] = 0.0
      acceptedj[i][1] = 0.0
      acceptedj[i][2] = 0.0
    end
    for i = 1:4
      cacheRootsi[i] = 0.0
      cacheRootsj[i] = 0.0
    end
    #find positive zeros f=+-Δ
    bi = aii * xi + aij * xj + uij
    ci = aij * (aji * xi + uji) - ajj * (aii * xi + uij)
    αi = -ajj - aii
    βi = aii * ajj - aij * aji
    bj = ajj * xj + aji * xi + uji
    cj = aji * (aij * xj + uij) - aii * (ajj * xj + uji)
    αj = -aii - ajj
    βj = ajj * aii - aji * aij
    coefi = NTuple{3,Float64}((βi * quani - ci, αi * quani - bi, quani))
    coefi2 = NTuple{3,Float64}((-βi * quani - ci, -αi * quani - bi, -quani))
    coefj = NTuple{3,Float64}((βj * quanj - cj, αj * quanj - bj, quanj))
    coefj2 = NTuple{3,Float64}((-βj * quanj - cj, -αj * quanj - bj, -quanj))
    resi1, resi2 = quadRootv2(coefi)
    resi3, resi4 = quadRootv2(coefi2)
    resj1, resj2 = quadRootv2(coefj)
    resj3, resj4 = quadRootv2(coefj2)
    #construct intervals
    constructIntrval(acceptedi, resi1, resi2, resi3, resi4)
    constructIntrval(acceptedj, resj1, resj2, resj3, resj4)
    #find best H (largest overlap)
    ki = 3
    kj = 3
    while true
      currentLi = acceptedi[ki][1]
      currentHi = acceptedi[ki][2]
      currentLj = acceptedj[kj][1]
      currentHj = acceptedj[kj][2]
      if currentLj <= currentLi < currentHj || currentLi <= currentLj < currentHi#resj[end][1]<=resi[end][1]<=resj[end][2] || resi[end][1]<=resj[end][1]<=resi[end][2] #overlap
        h = min(currentHi, currentHj)
        if h == Inf && (currentLj != 0.0 || currentLi != 0.0) # except case both  0 sols  h1=(0,Inf)
          h = max(currentLj, currentLi) #ft-simt in case they are both ver small...elaborate on his later
        end
        if h == Inf # both lower bounds ==0  --> zero sols for both
          qi = getQfromAsymptote(simt, xi, βi, ci, αi, bi)
          qj = getQfromAsymptote(simt, xj, βj, cj, αj, bj)
        else
          Δ = (1 - h * aii) * (1 - h * ajj) - h * h * aij * aji
          if Δ != 0.0
            qi = ((1 - h * ajj) * (xi + h * uij) + h * aij * (xj + h * uji)) / Δ
            qj = ((1 - h * aii) * (xj + h * uji) + h * aji * (xi + h * uij)) / Δ
          else
            qi, qj, h, maxIter = IterationH(h, xi, quani, xj, quanj, aii, ajj, aij, aji, uij, uji)
            if maxIter == 0
              return false
            end
          end
        end
        break
      else
        if currentHj == 0.0 && currentHi == 0.0#empty
          ki -= 1
          kj -= 1
        elseif currentHj == 0.0 && currentHi != 0.0
          kj -= 1
        elseif currentHi == 0.0 && currentHj != 0.0
          ki -= 1
        else #both non zero
          if currentLj < currentLi#remove last seg of resi
            ki -= 1
          else #remove last seg of resj
            kj -= 1
          end
        end
      end
    end # end while       
    q[j][0] = qj
    tq[j] = simt
    q[index][0] = qi# store back helper vars
    trackSimul[1] += 1 # do not have to recomputeNext if qi never changed
  end
  return iscycle
end
function getQfromAsymptote(simt, x::Float64, β::P, c::P, α::P, b::P) where {P<:Union{BigFloat,Float64}} 
  q = 0.0
  if β == 0.0 && c != 0.0
    println("report bug: β==0 && c!=0.0: this leads to asym=Inf while this code came from asym<Delta at $simt")
  elseif β == 0.0 && c == 0.0
    if α == 0.0 && b != 0.0
      println("report bug: α==0 && b!=0.0 this leads to asym=Inf while this code came from asym<Delta at $simt")
    elseif α == 0.0 && b == 0.0
      q = x
    else
      q = x + b / α
    end
  else
    q = x + c / β #asymptote...even though not guaranteed to be best
  end
  q
end
function iterationH(h::Float64, xi::Float64, quani::Float64, xj::Float64, quanj::Float64, aii::Float64, ajj::Float64, aij::Float64, aji::Float64, uij::Float64, uji::Float64)
  qi, qj = 0.0, 0.0
  maxIter = 1000
  while (abs(qi - xi) > 1 * quani || abs(qj - xj) > 1 * quanj) && (maxIter > 0)
    maxIter -= 1
    h1 = h * (0.99 * quani / abs(qi - xi))
    h2 = h * (0.99 * quanj / abs(qj - xj))
    h = min(h1, h2)
    h_three = h
    Δ = (1 - h * aii) * (1 - h * ajj) - h * h * aij * aji
    if Δ != 0
      qi = ((1 - h * ajj) * (xi + h * uij) + h * aij * (xj + h * uji)) / Δ
      qj = ((1 - h * aii) * (xj + h * uji) + h * aji * (xi + h * uij)) / Δ
    end
    # if maxIter < 1 println("maxiter of updateQ      = ",maxIter) end
  end
  qi, qj, h, maxIter
end