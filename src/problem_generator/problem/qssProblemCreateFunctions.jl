
"""
    createExactJacFun(otherCode::Expr,Exactjac:: Dict{Expr,Union{Float64,Int,Symbol,Expr}},funName::Symbol,f) 

 constructs the exact jacobian entries as a function from the existing dictionary Exactjac resulted from extractJacDep and extractJacDepLoop functions. \n

  From the input dictionary exactJac, we see that the key is always a expression that hold a tuple of size 2: the first element is going to be the index 'i' and the second element if the index 'j'. The value corresponding to this key, which is the exact jacobian entry, is put in a cache. i.e the function maps the keys of the dictionary to their values using an 'if-statement.. This approach does not depend on the size of the problem.
    # Example:   
```jldoctest
using QuantizedSystemSolver
exacteJacExpr=Dict{Expr,Union{Float64,Int,Symbol,Expr}}(:((1, 1)) => :(-2.0 * (q[2])[0]), :((1, 2)) => :(1 - 2.0 * (q[1])[0]),:(((2, 9), i - 1)) => :((q[i])[0]), :(((2, 9), i)) => :((q[i - 1])[0]),:((10, 10)) => -1, :((10, 1)) => 1);
exactJac=QuantizedSystemSolver.createExactJacFun(:(),exacteJacExpr,:f,0);
exactJac

# output

:(function exactJacf(q::Vector{Taylor0}, p::Vector{Float64}, cache::AbstractVector{Float64}, i::Int, j::Int, t::Float64, f_)
      (if i == 0
              return nothing
          elseif i == 1 && j == 1
              cache[1] = -2.0 * (q[2])[0]
              return nothing
          elseif i == 10 && j == 10
              cache[1] = -1
              return nothing
          elseif 2 <= i <= 9 && j == i - 1
              cache[1] = (q[i])[0]
              return nothing
          elseif i == 1 && j == 2
              cache[1] = 1 - 2.0 * (q[1])[0]
              return nothing
          elseif 2 <= i <= 9 && j == i
              cache[1] = (q[i - 1])[0]
              return nothing
          elseif i == 10 && j == 1
              cache[1] = 1
              return nothing
          end,)
  end)
```
"""
function createExactJacFun(passed_otherCode::Expr,Exactjac:: Dict{Expr,Union{Float64,Int,Symbol,Expr}},funName::Symbol,f) 
  ss="if i==0 return nothing\n"
  for dictElement in Exactjac
      if dictElement[1].args[1] isa Int
          ss*="elseif i==$(dictElement[1].args[1]) && j==$(dictElement[1].args[2]) \n"
          ss*="cache[1]=$(dictElement[2]) \n"
          ss*=" return nothing \n"
      elseif dictElement[1].args[1] isa Expr
          ss*="elseif $(dictElement[1].args[1].args[1])<=i<=$(dictElement[1].args[1].args[2]) && j==$(dictElement[1].args[2]) \n"
          ss*="cache[1]=$(dictElement[2]) \n"
          ss*=" return nothing \n"
      end     
  end 
  ss*=" end \n"        
  myex1=Meta.parse(ss)
  Base.remove_linenums!(myex1)

  otherCode=copy(passed_otherCode)# to avoid changing the original
  push!(otherCode.args,myex1)
  
  Base.remove_linenums!(otherCode)
  def1=Dict{Symbol,Any}() 
  def1[:head] = :function
  def1[:name] = Symbol(:exactJac,funName)  
  #def1[:args] = [:(q::Vector{Taylor0}),:(p::Vector{Float64}),:(cache::AbstractVector{Float64}),:(i::Int),:(j::Int),:(t::Float64),:(f_)]
  #def1[:args] = [:(q::Vector{Taylor0}),:(p::Vector{Any}),:(cache::AbstractVector{Float64}),:(i::Int),:(j::Int),:(t::Float64),:(f_)]
  def1[:args] = [:(q::Vector{Taylor0}),:(p),:(cache::AbstractVector{Float64}),:(i::Int),:(j::Int),:(t::Float64),:(f_)]
  def1[:body] = otherCode
  functioncode1=combinedef(def1)
end
 


function createEqFun(passed_otherCode::Expr,equs::Dict{Union{Int,Expr},Union{Int,Symbol,Expr}},zcequs::Vector{Expr},eventequs::Vector{Expr},fname::Symbol,f)# where{F}
  #allEpxpr=Expr(:block)
  ##############diffEqua###############
  s="if i==0 return nothing\n"  # :i is the mute var
  for elmt in equs
      Base.remove_linenums!(elmt[1])
      Base.remove_linenums!(elmt[2])
      if elmt[1] isa Int
          s*="elseif i==$(elmt[1]) $(elmt[2]) ;return nothing\n"
      end
      if elmt[1] isa Expr
          s*="elseif $(elmt[1].args[1])<=i<=$(elmt[1].args[2]) $(elmt[2]) ;return nothing\n"
      end
  end
  s*=" end "
  myex1=Meta.parse(s)


  
  #push!(otherCode.args,myex1) 
  otherCodeCopy=copy(passed_otherCode)# to avoid changing the original
  push!(otherCodeCopy.args,myex1)
   ##############ZCF###################
  if length(zcequs)>0
      s="if zc==1  $(zcequs[1]) ;return nothing"
      for i=2:length(zcequs)
          s*= " elseif zc==$i $(zcequs[i]) ;return nothing"
      end
      s*= " end "
      myex2=Meta.parse(s) 
      push!(otherCodeCopy.args,myex2)
  end
   #############events#################
  if length(eventequs)>0
      s= "if ev==1  $(eventequs[1]) ;return nothing"
      for i=2:length(eventequs)
          s*= " elseif ev==$i $(eventequs[i]) ;return nothing"
      end
      s*= " end "
      myex3=Meta.parse(s)
      push!(otherCodeCopy.args,myex3)
  end
 
  Base.remove_linenums!(otherCodeCopy)
  def=Dict{Symbol,Any}()
  def[:head] = :function
  def[:name] = fname  
  #def[:args] = [:(i::Int),:(zc::Int),:(ev::Int),:(q::Vector{Taylor0}),:(p::Vector{Float64}), :(t::Taylor0),:(cache::Vector{Taylor0}),:(f_)]  
  #def[:args] = [:(i::Int),:(zc::Int),:(ev::Int),:(q::Vector{Taylor0}),:(p::Vector{Any}), :(t::Taylor0),:(cache::Vector{Taylor0}),:(f_)]  
  def[:args] = [:(i::Int),:(zc::Int),:(ev::Int),:(q::Vector{Taylor0}),:(p), :(t::Taylor0),:(cache::Vector{Taylor0}),:(f_)]  
  def[:body] = otherCodeCopy 
  functioncode=combinedef(def)
 # @show functioncode

end



