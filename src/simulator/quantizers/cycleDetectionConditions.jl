# Order 1
function detect1(::Val{0},ẋi,dxP,dxi,ẋj,dxj)
  return false
end

"""
    detect1(::Val{1}, ẋi, dxP, dxi, ẋj, dxj)

cycle detection in order 1 that uses mechanism type 1. `ẋi` is used and only change of direction is checked.

# Arguments
- `::Val{1}`: Type parameter indicating the cycle detection mechanism type 1.
- `ẋi`: Derivative of the first variable.
- `dxP`: derivative the first variable used to update the quantized variable `q`.
- `dxi`: prediction of the derivative of the first variable.
- `ẋj`: Derivative of the second variable.
- `dxj`: prediction of the derivative of the second variable.

# Returns
- A boolean indicating whether a cycle is present.
```julia
return  (dxj*ẋj)<0.0 && (dxi*ẋi)<0.0
```

# Notes
This function is used internally by the integrator to detect the presence of a cycle.
"""
function detect1(::Val{1},ẋi,dxP,dxi,ẋj,dxj)
  return  (dxj*ẋj)<0.0 && (dxi*ẋi)<0.0
end


"""
    detect1(::Val{2}, ẋi, dxP, dxi, ẋj, dxj)

cycle detection in order 1 that uses mechanism type 2. `dxP` is used and only change of direction is checked.
# Returns
- A boolean indicating whether a cycle is present.
```julia
return  (dxj*ẋj)<0.0 && (dxi*dxP)<0.0
```
"""
function detect1(::Val{2},ẋi,dxP,dxi,ẋj,dxj)
  return (dxj*ẋj)<0.0 && (dxi*dxP)<0.0
end

"""
    detect1(::Val{3}, ẋi, dxP, dxi, ẋj, dxj)

cycle detection in order 1 that uses mechanism type 3. `ẋi` is used and significant changes in the derivatives are also checked.
# Returns
- A boolean indicating whether a cycle is present.
```julia
return  abs(dxj-ẋj)>abs(dxj+ẋj)/2 && abs(dxi-ẋi)>abs(dxi+ẋi)/2
```
"""
function detect1(::Val{3},ẋi,dxP,dxi,ẋj,dxj)
  return abs(dxj-ẋj)>abs(dxj+ẋj)/2 && abs(dxi-ẋi)>abs(dxi+ẋi)/2
end

"""
    detect1(::Val{4}, ẋi, dxP, dxi, ẋj, dxj)

cycle detection in order 1 that uses mechanism type 3. `dxP` is used and significant changes in the derivatives are also checked.
# Returns
- A boolean indicating whether a cycle is present.
```julia
return  abs(dxj-ẋj)>abs(dxj+ẋj)/2 && abs(dxi-dxP)>abs(dxi+dxP)/2
```
"""
function detect1(::Val{4},ẋi,dxP,dxi,ẋj,dxj)
  return abs(dxj-ẋj)>abs(dxj+ẋj)/2 && abs(dxi-dxP)>abs(dxi+dxP)/2
end



function detect1(::Val{5},ẋi,dxP,dxi,ẋj,dxj)
  return ((abs(dxj)>3*abs(ẋj)||3*abs(dxj)<abs(ẋj)) && dxi*ẋi<0.0) || ((abs(dxi)>3*abs(ẋi)||3*abs(dxi)<abs(ẋi)) && dxj*ẋj<0.0) || ((dxj*ẋj)<0.0 && (dxi*ẋi)<0.0)
end
function detect1(::Val{6},ẋi,dxP,dxi,ẋj,dxj)
  return ((abs(dxj)>3*abs(ẋj)||3*abs(dxj)<abs(ẋj)) && dxi*dxP<0.0) ||((abs(dxi)>3*abs(dxP)||3*abs(dxi)<abs(dxP)) && dxj*ẋj<0.0) || ((dxj*ẋj)<0.0 && (dxi*dxP)<0.0)
end
#cases to remove: one value null and other very small
#= 
function detect1(::Val{7},ẋi,dxP,dxi,ẋj,dxj)
  if (dxj*ẋj)<=0.0 && (dxi*ẋi)<=0.0
    if (dxj)==0.0 && abs(ẋj)<1e-6
      return false
    end
    if (ẋj)==0.0 && abs(dxj)<1e-6
      return false
    end
    if (dxi)==0.0 && abs(ẋi)<1e-6
      return false
    end
    if (ẋi)==0.0 && abs(dxi)<1e-6
      return false
    end
      return true
  else  
    return false
  end
  end
function detect1(::Val{8},ẋi,dxP,dxi,ẋj,dxj)
  if (dxj*ẋj)<=0.0 && (dxi*dxP)<=0.0
    if (dxj)==0.0 && abs(ẋj)<1e-6
      return false
    end
    if (ẋj)==0.0 && abs(dxj)<1e-6
      return false
    end
    if (dxi)==0.0 && abs(dxP)<1e-6
      return false
    end
    if (dxP)==0.0 && abs(dxi)<1e-6
      return false
    end
    return true
  else  
    return false
  end
end
  
function detect1(::Val{9},ẋi,dxP,dxi,ẋj,dxj)
  if abs(dxj-ẋj)>abs(dxj+ẋj)/2 && abs(dxi-ẋi)>abs(dxi+ẋi)/2
    if (dxj)==0.0 && abs(ẋj)<1e-6
      return false
    end
    if (ẋj)==0.0 && abs(dxj)<1e-6
      return false
    end
    if (dxi)==0.0 && abs(ẋi)<1e-6
      return false
    end
    if (ẋi)==0.0 && abs(dxi)<1e-6
      return false
    end
    return true
  else  
    return false
  end
end

function detect1(::Val{10},ẋi,dxP,dxi,ẋj,dxj)
  if abs(dxj-ẋj)>abs(dxj+ẋj)/2 && abs(dxi-dxP)>abs(dxi+dxP)/2
    if (dxj)==0.0 && abs(ẋj)<1e-6
      return false
    end
    if (ẋj)==0.0 && abs(dxj)<1e-6
      return false
    end
    if (dxi)==0.0 && abs(dxP)<1e-6
      return false
    end
    if (dxP)==0.0 && abs(dxi)<1e-6
      return false
    end
    return true
  else  
    return false
  end
end

function detect1(::Val{11},ẋi,dxP,dxi,ẋj,dxj)
  return (dxj*ẋj)<=0.0 && (dxi*ẋi)<=0.0
end
function detect1(::Val{12},ẋi,dxP,dxi,ẋj,dxj)
  return (dxj*ẋj)<=0.0 && (dxi*dxP)<=0.0
end
 =#


function detect2(::Val{0},xi1,dxi,dxithrow,xi2,ddxi,ddxithrow,βidir,βidth,xj1,dxj,xj2,ddxj,dqjplus,recentjDir,dirI)
  return false
end

"""
    detect2(::Val{1},xi1,dxi,dxithrow,xi2,ddxi,ddxithrow,βidir,βidth,xj1,dxj,xj2,ddxj,dqjplus,recentjDir,dirI)

cycle detection in order 2 that uses mechanism type 1. `dxithrow` and `ddxithrow` are used. changes in the derivatives are checked.

# Returns
- A boolean indicating whether a cycle is present.
```julia
  if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  
    return (abs(dxi-dxithrow)>(abs(dxi+dxithrow)/2) || abs(ddxi-ddxithrow)>(abs(ddxi+ddxithrow)/2)) 
  else  
    return false
  end
```
"""
function detect2(::Val{1},xi1,dxi,dxithrow,xi2,ddxi,ddxithrow,βidir,βidth,xj1,dxj,xj2,ddxj,dqjplus,recentjDir,dirI)
  if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  
    return (abs(dxi-dxithrow)>(abs(dxi+dxithrow)/2) || abs(ddxi-ddxithrow)>(abs(ddxi+ddxithrow)/2)) 
  else  
    return false
  end
end

"""
    detect2(::Val{2},xi1,dxi,dxithrow,xi2,ddxi,ddxithrow,βidir,βidth,xj1,dxj,xj2,ddxj,dqjplus,recentjDir,dirI)

cycle detection in order 2 that uses mechanism type 2. `xi1` and `xi2` are used. changes in the derivatives are checked.

# Returns
- A boolean indicating whether a cycle is present.
```julia
 if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  
    return (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2))  
  else  
    return false
  end
```
"""
function detect2(::Val{2},xi1,dxi,dxithrow,xi2,ddxi,ddxithrow,βidir,βidth,xj1,dxj,xj2,ddxj,dqjplus,recentjDir,dirI)
  if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  
    return (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2))  
  else  
    return false
  end
end



"""
    detect2(::Val{3},xi1,dxi,dxithrow,xi2,ddxi,ddxithrow,βidir,βidth,xj1,dxj,xj2,ddxj,dqjplus,recentjDir,dirI)

cycle detection in order 2 that uses mechanism type 3. `xi1` and `xi2` are used. changes in the derivatives are checked, changes in direction are also checked.

# Returns
- A boolean indicating whether a cycle is present.
```julia
  if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  || dqjplus*recentjDir<0.0 
    return (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)) || βidir*dirI<0.0
  else  
    return false
  end
```
"""
function detect2(::Val{3},xi1,dxi,dxithrow,xi2,ddxi,ddxithrow,βidir,βidth,xj1,dxj,xj2,ddxj,dqjplus,recentjDir,dirI)
  if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  || dqjplus*recentjDir<0.0 #(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
    return (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)) || βidir*dirI<0.0
  else  
    return false
  end
end

"""
    detect2(::Val{4},xi1,dxi,dxithrow,xi2,ddxi,ddxithrow,βidir,βidth,xj1,dxj,xj2,ddxj,dqjplus,recentjDir,dirI)

cycle detection in order 2 that uses mechanism type 3. `dxithrow` and `ddxithrow` are used. changes in the derivatives are checked, changes in direction are also checked.

# Returns
- A boolean indicating whether a cycle is present.
```julia
  if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  || dqjplus*recentjDir<0.0 
    return (abs(dxi-dxithrow)>(abs(dxi+dxithrow)/2) || abs(ddxi-ddxithrow)>(abs(ddxi+ddxithrow)/2)) || βidir*βidth<0.0
  else  
    return false
  end
```
"""
function detect2(::Val{4},xi1,dxi,dxithrow,xi2,ddxi,ddxithrow,βidir,βidth,xj1,dxj,xj2,ddxj,dqjplus,recentjDir,dirI) 
  if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  || dqjplus*recentjDir<0.0 #(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
    return (abs(dxi-dxithrow)>(abs(dxi+dxithrow)/2) || abs(ddxi-ddxithrow)>(abs(ddxi+ddxithrow)/2)) || βidir*βidth<0.0
  else  
    return false
  end
end


function detect2(::Val{5},xi1,dxi,dxithrow,xi2,ddxi,ddxithrow,βidir,βidth,xj1,dxj,xj2,ddxj,dqjplus,recentjDir,dirI)
  return dqjplus*recentjDir<=0.0 && βidir*dirI<=0.0
end

function detect2(::Val{6},xi1,dxi,dxithrow,xi2,ddxi,ddxithrow,βidir,βidth,xj1,dxj,xj2,ddxj,dqjplus,recentjDir,dirI)
    return dqjplus*recentjDir<=0.0 && βidth*βidir<=0.0
end