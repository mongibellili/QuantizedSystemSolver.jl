

"""
    minPosRoot(c::Float64,b::Float64, ::Val{1})

Finds the minimum positive root for a linear equation represented by the coefficients.

# Returns
- The minimum positive root of the linear equation. 
"""
function minPosRoot(c::Float64,b::Float64, ::Val{1}) #call from reComputeNextTime
  mpr = -1
  if b == 0
    mpr = Inf
  else
    mpr = -c / b
  end
  if mpr < 0
    mpr = Inf
  end
  # println("mpr inside minPosRoot in utils= ",mpr) 
  return mpr
end

"""
    minPosRoot(coeff::Taylor0, ::Val{1})

Finds the minimum positive root for a linear equation represented by the Taylor series coefficients.

# Arguments
- `coeff::Taylor0`: The Taylor series coefficients of the linear equation.
- `::Val{1}`: A type parameter indicating the order.

# Returns
- The minimum positive root of the linear equation.
"""
function minPosRoot(coeff::Taylor0, ::Val{1}) # call from compute event order1
  minPosRoot(coeff[0],coeff[1], Val(1))
end

"""
    minPosRoot(c::Float64,b::Float64,a::Float64,::Val{2})

Finds the minimum positive root for a quadratic equation represented by the coefficients.

# Returns
- The minimum positive root of the quadratic equation.
"""
function minPosRoot(c::Float64,b::Float64,a::Float64,::Val{2}) #called from updateQ order2 and reComputeNextTime
  mpr = -1
  if a == 0 || (1000000 * abs(a)) < abs(b)# coef3 is the coef of t^2
    if b == 0
      mpr = Inf
    else
      if a==0
        mpr = -c / b
      else
        if  0 <-c / b < 1e7 # neglecting a small 'a' then having a large h would cause an error  'a*h*h' because large
          mpr = -c / b
        end
      end
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

"""
    minPosRoot(coeff::Taylor0, ::Val{2})

Finds the minimum positive root for a quadratic equation represented by the Taylor series coefficients.

# Arguments
- `coeff::Taylor0`: The Taylor series coefficients of the quadratic equation.
- `::Val{2}`: A type parameter indicating the order.

# Returns
- The minimum positive root of the quadratic equation.
"""
function minPosRoot(coeff::Taylor0, ::Val{2}) #  call from compute event(::Val{2})
  minPosRoot(coeff[0],coeff[1],coeff[2], Val(2))
end


