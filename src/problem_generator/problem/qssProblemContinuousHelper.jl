
"""
    createContEqFun(otherCode::Expr,equs::Dict{Union{Int,Expr},Union{Int,Symbol,Expr}},fname::Symbol,f)

 constructs one function from all differential equations in the problem, which are transformed and stored in a dictionary in the NLodeProblemFunc function.\n
 
  The function maps the keys of the dictionary to their values using an 'if-statement.
    # Example:   
```jldoctest
using QuantizedSystemSolver
equs = Dict{Union{Int,Expr},Union{Int,Symbol,Expr}}(10 => :(subT(q[1], q[10], cache[1])), :((2, 9)) => :(mulT(q[i], q[i - 1], cache[1])), 1 => :(subT(q[2], mulTT(2.0, q[1], q[2], cache[2], cache[3]), cache[1])));
diffEqfun=QuantizedSystemSolver.createContEqFun(:(),equs,:f,0); 
diffEqfun 

# output

:(function f(i::Int, q::Vector{Taylor0}, p::Vector{Float64}, t::Taylor0, cache::Vector{Taylor0}, f_)    
      (if i == 0
              return nothing
          elseif i == 10
              subT(q[1], q[10], cache[1])
              return nothing
          elseif 2 <= i <= 9
              mulT(q[i], q[i - 1], cache[1])
              return nothing
          elseif i == 1
              subT(q[2], mulTT(2.0, q[1], q[2], cache[2], cache[3]), cache[1])
              return nothing
          end,)
  end)
```

"""
function createContEqFun(otherCode::Expr,equs::Dict{Union{Int,Expr},Union{Int,Symbol,Expr}},fname::Symbol,f) 
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
  Base.remove_linenums!(myex1)
  otherCodeCopy=copy(otherCode)# to avoid changing the original
  push!(otherCodeCopy.args,myex1)
  Base.remove_linenums!(otherCodeCopy)
  def=Dict{Symbol,Any}()
  def[:head] = :function
  def[:name] = fname   
  def[:args] = [:(i::Int),:(q::Vector{Taylor0}),:(p::Vector{Float64}),:(t::Taylor0),:(cache::Vector{Taylor0}),:(f_)]
  def[:body] = otherCodeCopy#push!(myex1.args,funcEx.args[1])  
  functioncode=combinedef(def)
end






