#using StaticArrays
function minPosRoot(coeff::SVector{2,Float64}, ::Val{1}) # coming from val(1) means coef has x and derx only...size 2 #call from reComputeNextTime
  mpr = -1
  if coeff[2] == 0
    mpr = Inf
  else
    mpr = -coeff[1] / coeff[2]
  end
  if mpr < 0
    mpr = Inf
  end
  # println("mpr inside minPosRoot in utils= ",mpr)
  return mpr
end
function minPosRoot(coeff::Taylor0, ::Val{1}) # coming from val(1) means coef has x and derx only...size 2 call from compute event order1
  mpr = -1
  if coeff[1] == 0
    mpr = Inf
  else
    mpr = -coeff[0] / coeff[1]
  end
  if mpr < 0
    mpr = Inf
  end
  return mpr
end
function minPosRootv1(coeff::NTuple{3,Float64}) #called from updateQ order2
  a = coeff[1]
  b = coeff[2]
  c = coeff[3]
  mpr = -1 #coef1=c, coef2=b, coef3=a
  if a == 0 || 10000 * abs(a) < abs(b)# coef3 is the coef of t^2
    if b == 0 # julia allows to divide by zero and result is Inf ...so no need to check
      mpr = Inf
    else
      mpr = -c / b
    end
    if mpr < 0
      mpr = Inf
    end
  else
    #double disc;
    disc = b * b - 4 * a * c#b^2-4ac
    if disc < 0 # no real roots
      mpr = Inf
    else
      #double sd, r1;
      sd = sqrt(disc)
      r1 = (-b + sd) / (2 * a)
      if r1 > 0
        mpr = r1
      else
        mpr = Inf
      end
      r1 = (-b - sd) / (2 * a)
      if ((r1 > 0) && (r1 < mpr))
        mpr = r1
      end
    end
  end
  return mpr
end
function minPosRoot(coeff::SVector{3,Float64}, ::Val{2}) # credit goes to github.com/CIFASIS/qss-solver # call from reComputeNextTime order2
  mpr = -1
  a = coeff[3]
  b = coeff[2]
  c = coeff[1]  # a is coeff 3 because in taylor representation 1 is var 2 is der 3 is derder
  if a == 0 || (10000 * abs(a)) < abs(b)# coef3 is the coef of t^2
    if b == 0
      mpr = Inf
    else
      mpr = -c / b
    end
    if mpr < 0
      mpr = Inf
    end
  else
    #double disc;
    disc = b * b - 4 * a * c#b^2-4ac
    if disc < 0 # no real roots
      mpr = Inf
    else
      #double sd, r1;
      sd = sqrt(disc)
      r1 = (-b + sd) / (2 * a)
      if r1 > 0
        mpr = r1
      else
        mpr = Inf
      end
      r1 = (-b - sd) / (2 * a)
      if ((r1 > 0) && (r1 < mpr))
        mpr = r1
      end
    end
  end
  return mpr
end
function minPosRoot(coeff::Taylor0, ::Val{2}) # credit goes to github.com/CIFASIS/qss-solver # call from compute event(::Val{2})
  mpr = -1
  a = coeff[2]
  b = coeff[1]
  c = coeff[0]  # a is coeff 3 because in taylor representation 1 is var 2 is der 3 is derder
  if a == 0 || (100000 * abs(a)) < abs(b)# coef3 is the coef of t^2
    if b == 0
      mpr = Inf
    else
      mpr = -c / b
    end
    if mpr < 0
      mpr = Inf
    end
  else
    #double disc;
    disc = b * b - 4 * a * c#b^2-4ac
    if disc < 0 # no real roots
      mpr = Inf
    else
      #double sd, r1;
      sd = sqrt(disc)
      r1 = (-b + sd) / (2 * a)
      if r1 > 0
        mpr = r1
      else
        mpr = Inf
      end
      r1 = (-b - sd) / (2 * a)
      if ((r1 > 0) && (r1 < mpr))
        mpr = r1
      end
    end
  end
  if DEBUG && mpr != Inf
    sl = c + mpr * b + a * mpr * mpr
    @show sl
  end
  return mpr
end
function quadRootv2(coeff::NTuple{3,Float64}) # call from mliqss1/simultupdate1
  mpr = (-1.0, -1.0) #size 2 to use mpr[2] in quantizer
  a = coeff[1]
  b = coeff[2]
  c = coeff[3]
  if a == 0.0 || (1e7 * abs(a)) < abs(b)# coef3 is the coef of t^2
    if b != 0.0
      if 0 < -c / b < 1e7  # neglecting a small 'a' then having a large h would cause an error  'a*h*h' because large
        mpr = (-c / b, -1.0)
      end
    end
  elseif b == 0.0
    if -c / a > 0
      mpr = (sqrt(-c / a), -1.0)
    end
  elseif c == 0.0
    mpr = (-1.0, -b / a)
  else
    #double disc;
    Δ = 1.0 - 4.0 * c * a / (b * b)
    if Δ > 0.0
      sq = sqrt(Δ)
      r1 = -0.5 * (1.0 + sq) * b / a
      r2 = -0.5 * (1.0 - sq) * b / a
      mpr = (r1, r2)
    elseif Δ == 0.0
      r1 = -0.5 * b / a
      mpr = (r1, r1 - 1e-12)
    end
  end
  return mpr
end
function constructIntrval2(cache::Vector{Vector{Float64}}, t1::Float64, t2::Float64)
  h1 = min(t1, t2)
  h4 = max(t1, t2)
  # res=((0.0,h1),(h4,Inf))
  cache[1][1] = 0.0
  cache[1][2] = h1
  cache[2][1] = h4
  cache[2][2] = Inf
  return nothing
end
#constructIntrval(tup2::Tuple{},tup::NTuple{2, Float64})=constructIntrval(tup,tup2)
#merge has 3 elements
function constructIntrval3(cache::Vector{Vector{Float64}}, t1::Float64, t2::Float64, t3::Float64)
  #t1=tup[1];t2=tup[2];t3=tup2[1]
  t12 = min(t1, t2)
  t21 = max(t1, t2)
  if t3 < t12
    # res=((0.0,t3),(t12,t21))
    cache[1][1] = 0.0
    cache[1][2] = t3
    cache[2][1] = t12
    cache[2][2] = t21
  else
    if t3 < t21
      # res=((0.0,t12),(t3,t21))
      cache[1][1] = 0.0
      cache[1][2] = t12
      cache[2][1] = t3
      cache[2][2] = t21
    else
      #res=((0.0,t12),(t21,t3))
      cache[1][1] = 0.0
      cache[1][2] = t12
      cache[2][1] = t21
      cache[2][2] = t3
    end
  end
  return nothing
end
#merge has 4 elements
function constructIntrval4(cache::Vector{Vector{Float64}}, t1::Float64, t2::Float64, t3::Float64, t4::Float64,)
  # t1=tup[1];t2=tup[2];t3=tup2[1];t4=tup2[2]
  h1 = min(t1, t2, t3, t4)
  h4 = max(t1, t2, t3, t4)
  if t1 == h1
    if t2 == h4
      # res=((0,t1),(min(t3,t4),max(t3,t4)),(t2,Inf))
      cache[1][1] = 0.0
      cache[1][2] = t1
      cache[2][1] = min(t3, t4)
      cache[2][2] = max(t3, t4)
      cache[3][1] = t2
      cache[3][2] = Inf
    elseif t3 == h4
      #res=((0,t1),(min(t2,t4),max(t2,t4)),(t3,Inf))
      cache[1][1] = 0.0
      cache[1][2] = t1
      cache[2][1] = min(t2, t4)
      cache[2][2] = max(t2, t4)
      cache[3][1] = t3
      cache[3][2] = Inf
    else
      # res=((0,t1),(min(t2,t3),max(t2,t3)),(t4,Inf))
      cache[1][1] = 0.0
      cache[1][2] = t1
      cache[2][1] = min(t3, t2)
      cache[2][2] = max(t3, t2)
      cache[3][1] = t4
      cache[3][2] = Inf
    end
  elseif t2 == h1
    if t1 == h4
      #res=((0,t2),(min(t3,t4),max(t3,t4)),(t1,Inf))
      cache[1][1] = 0.0
      cache[1][2] = t2
      cache[2][1] = min(t3, t4)
      cache[2][2] = max(t3, t4)
      cache[3][1] = t1
      cache[3][2] = Inf
    elseif t3 == h4
      #res=((0,t2),(min(t1,t4),max(t1,t4)),(t3,Inf))
      cache[1][1] = 0.0
      cache[1][2] = t2
      cache[2][1] = min(t1, t4)
      cache[2][2] = max(t1, t4)
      cache[3][1] = t3
      cache[3][2] = Inf
    else
      # res=((0,t2),(min(t1,t3),max(t1,t3)),(t4,Inf))
      cache[1][1] = 0.0
      cache[1][2] = t2
      cache[2][1] = min(t3, t1)
      cache[2][2] = max(t3, t1)
      cache[3][1] = t4
      cache[3][2] = Inf
    end
  elseif t3 == h1
    if t1 == h4
      #res=((0,t3),(min(t2,t4),max(t2,t4)),(t1,Inf))
      cache[1][1] = 0.0
      cache[1][2] = t3
      cache[2][1] = min(t2, t4)
      cache[2][2] = max(t2, t4)
      cache[3][1] = t1
      cache[3][2] = Inf
    elseif t2 == h4
      #res=((0,t3),(min(t1,t4),max(t1,t4)),(t2,Inf))
      cache[1][1] = 0.0
      cache[1][2] = t3
      cache[2][1] = min(t1, t4)
      cache[2][2] = max(t1, t4)
      cache[3][1] = t2
      cache[3][2] = Inf
    else
      #res=((0,t3),(min(t1,t2),max(t1,t2)),(t4,Inf))
      cache[1][1] = 0.0
      cache[1][2] = t3
      cache[2][1] = min(t1, t2)
      cache[2][2] = max(t1, t2)
      cache[3][1] = t4
      cache[3][2] = Inf
    end
  else
    if t1 == h4
      cache[1][1] = 0.0
      cache[1][2] = t4
      cache[2][1] = min(t3, t2)
      cache[2][2] = max(t3, t2)
      cache[3][1] = t1
      cache[3][2] = Inf
    elseif t2 == h4
      cache[1][1] = 0.0
      cache[1][2] = t4
      cache[2][1] = min(t3, t1)
      cache[2][2] = max(t3, t1)
      cache[3][1] = t2
      cache[3][2] = Inf
    else
      cache[1][1] = 0.0
      cache[1][2] = t4
      cache[2][1] = min(t1, t2)
      cache[2][2] = max(t1, t2)
      cache[3][1] = t3
      cache[3][2] = Inf
    end
  end
  return nothing
end
@inline function constructIntrval(cache::Vector{Vector{Float64}}, res1::Float64, res2::Float64, res3::Float64, res4::Float64)
  if res1 <= 0.0 && res2 <= 0.0 && res3 <= 0.0 && res4 <= 0.0
    cache[1][1] = 0.0
    cache[1][2] = Inf
  elseif res1 > 0.0 && res2 <= 0.0 && res3 <= 0.0 && res4 <= 0.0
    cache[1][1] = 0.0
    cache[1][2] = res1
  elseif res2 > 0.0 && res1 <= 0.0 && res3 <= 0.0 && res4 <= 0.0
    cache[1][1] = 0.0
    cache[1][2] = res2
  elseif res3 > 0.0 && res2 <= 0.0 && res1 <= 0.0 && res4 <= 0.0
    cache[1][1] = 0.0
    cache[1][2] = res3
  elseif res4 > 0.0 && res2 <= 0.0 && res3 <= 0.0 && res1 <= 0.0
    cache[1][1] = 0.0
    cache[1][2] = res4
  elseif res1 > 0.0 && res2 > 0.0 && res3 <= 0.0 && res4 <= 0.0
    constructIntrval2(cache, res1, res2)
  elseif res1 > 0.0 && res3 > 0.0 && res2 <= 0.0 && res4 <= 0.0
    constructIntrval2(cache, res1, res3)
  elseif res1 > 0.0 && res4 > 0.0 && res2 <= 0.0 && res4 <= 0.0
    constructIntrval2(cache, res1, res4)
  elseif res2 > 0.0 && res3 > 0.0 && res1 <= 0.0 && res4 <= 0.0
    constructIntrval2(cache, res2, res3)
  elseif res2 > 0.0 && res4 > 0.0 && res1 <= 0.0 && res3 <= 0.0
    constructIntrval2(cache, res2, res4)
  elseif res3 > 0.0 && res4 > 0.0 && res1 <= 0.0 && res2 <= 0.0
    constructIntrval2(cache, res3, res4)
  elseif res1 > 0.0 && res2 > 0.0 && res3 > 0.0 && res4 <= 0.0
    constructIntrval3(cache, res1, res2, res3)
  elseif res1 > 0.0 && res2 > 0.0 && res4 > 0.0 && res3 <= 0.0
    constructIntrval3(cache, res1, res2, res4)
  elseif res1 > 0.0 && res3 > 0.0 && res4 > 0.0 && res2 <= 0.0
    constructIntrval3(cache, res1, res3, res4)
  elseif res1 <= 0.0 && res3 > 0.0 && res4 > 0.0 && res2 > 0.0
    constructIntrval3(cache, res2, res3, res4)
  else
    constructIntrval4(cache, res1, res2, res3, res4)
  end
  return nothing
end
