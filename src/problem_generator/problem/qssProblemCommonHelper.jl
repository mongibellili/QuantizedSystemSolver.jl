"""
    changeExprToFirstValue(ex::Expr)

 changes an expression in the form u[1] to an expression in the form u[1][0] inside exact jacobian expressions and inside events, because linear coefficients (a_{ii}) do not have derivatives, and updates in events affect the value of a variable directly and there is no need to update its higher derivatives. It is called by the [`restoreRef`](@ref) function for jacobian expressions, and called by the [`handleEvents`](@ref) function for events.
# Example:
```jldoctest
using QuantizedSystemSolver 
ex=:(q[i - 1]) 
newEx=QuantizedSystemSolver.changeExprToFirstValue(ex)

# output

:((q[i - 1])[0])
```
"""
function changeExprToFirstValue(ex::Expr)
  newEx=postwalk(ex) do a  
      if a isa Expr && a.head == :ref && a.args[1]==:q  
          outerRef=Expr(:ref)
          push!(outerRef.args,a)
          push!(outerRef.args,:(0))
          a=outerRef
      end
      return a
  end
  newEx
end

"""
    symbolFromRef(el::Symbol,refEx::Union{Int64,Expr,Symbol})  

 gets a symbol qi, qiplusNumber, qiminusNumber, or qitimesNumber from a symbol i or expressions like i+Number, i+Number, i+Number.
 It is called by the [`changeVarNames_params`](@ref), the [`extractJacDepNormal`](@ref) and the [`extractJacDepLoop`](@ref) functions
  # Example:
```jldoctest
using QuantizedSystemSolver
ex=:(i - 1) 
newEx=QuantizedSystemSolver.symbolFromRef(:q,ex)
(ex,newEx)

# output

(:(i - 1), :qiminus1)
```
"""
function symbolFromRef(el::Symbol,refEx::Union{Int64,Expr,Symbol})  
  if refEx isa Expr 
    if refEx.args[1]==:+
      refEx=Symbol(el,(refEx.args[2]), "plus",(refEx.args[3]))
    elseif refEx.args[1]==:-
      refEx=Symbol(el,(refEx.args[2]), "minus",(refEx.args[3]))
    elseif refEx.args[1]==:*
      refEx=Symbol(el,(refEx.args[2]), "times",(refEx.args[3]))
    elseif refEx.args[1]==:/
      Error("parameter_index division is not implemented Yet")
    else
      error("the provided parameter_index $(refEx) is not implemented Yet in symbolFromRef")
    end
  else
    refEx=Symbol(el,(refEx))
  end
  return refEx
end


"""
    restoreRef(coefExpr,symDict)

This function is the opposite of symbolFromRef. After using the symbols in symbolic differentiation, it gets back expressions like p[i+Number] and q[i+Number][0] from symbols diplusNumber and qiplusNumber. Adding a zero to q variables is beacause q is a taylor variable while p is a vector.\n 
  # arguments:
- `coefExpr::Expr`: the expression to be changed
- `symDict::Dict{Symbol,Expr}`: the dictionary to store the translation of symbols of continous and discrete variables (q[i] <-> qi)

# Example:
```jldoctest
using QuantizedSystemSolver

symDict= Dict{Symbol, Expr}(:qi => :(q[i]), :q10 => :(q[10]), :q2 => :(q[2]), :qiminus1 => :(q[i - 1]), :q1 => :(q[1]))
coefExpr=:(1.5qiminus1) 
  newEx=QuantizedSystemSolver.restoreRef(coefExpr, symDict)


# output

:(1.5 * (q[i - 1])[0]) 
```
"""
function restoreRef(coefExpr,symDict)
  newEx=postwalk(coefExpr) do element# 
    if element isa Symbol && !(element in (:+,:-,:*,:/)) && haskey(symDict, element) 
        if String(element)[1] == 'q' 
          element=symDict[element]
          element=changeExprToFirstValue(element)# change q[1] to q[1][0]
        else  #element== :p or any other constant parameter
            element=symDict[element]  
        end
    end
    return element
  end#end postwalk
  newEx
end


function is_custom_function(expr::Expr)
    return expr.args[1] isa Symbol &&
           !(expr.args[1] in (:+, :-, :*, :/, :^, :%, :&, :|, :!, :(=), :(==), :!=, :<, :>, :<=, :>=)) &&
           !(isdefined(Base, expr.args[1]) && getfield(Base, expr.args[1]) isa Function)
end



function useNumber_inFuncalls(ex::Expr)
  newEx=postwalk(ex) do element#postwalk to change var names and parameters
      if element isa Expr && element.head == :call && element.args[1] isa Symbol 
        sym=element.args[1]
        if !(sym in (:+, :-, :*, :/, :^, :%, :&, :|, :!, :(=),:(==), :!=, :<, :>, :<=, :>=)) && !(isdefined(Base, sym) && getfield(Base, sym) isa Function)
          element=changeExprToFirstValue(element)
        end
      end
      return element
    end#end postwalk
  newEx
end



"""
    extractJacDepNormal(varNum::Int,rhs::Expr,jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},exactJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},jac_mode ::Symbol,symDict::Dict{Symbol,Expr},dD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 

Extracts the jacobian dependency (jac) as well as the exacte symbolic jacobian expression (exactJacExpr) and the dependency of state derivatives to discrete variables (dD), in the form of dictionaries, from the simple differential equations.\n

  For the continuous part, similar to [`extractJacDepNormal`](@ref)  function, this function sarts by looking for the 'i' in q[i] in the RHS and storing this 'i' in a jacSet for the varNum. Then, it changes q[i] to qi for symbolic differentiation. After finding ``\\frac{\\partial f_i}{\\partial q_i} `` as the exact jacobian entry, it changes back qi to q[i]. Also, any mute variable from the differential equations is changed to 'i'.\n
  For the discrete part, the function puts the index of the differential equation in a set, and stores this set in a dictionary dD with the key being the index of the discrete variable.
  # Example:
```jldoctest
using QuantizedSystemSolver
(varNum, rhs, jac, exactJacExpr, symDict, dD) = (1, :(p[2] - 2.0 * q[1] * p[2]), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(), Dict{Expr, Union{Float64, Int64, Expr, Symbol}}(), Dict{Symbol, Expr}(:q10 => :(q[10]), :p2 => :(p[2]), :qiminus1 => :(q[i - 1]), :p1 => :(p[1]), :q1 => :(q[1])), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}())
QuantizedSystemSolver.extractJacDepNormal(varNum, rhs, jac, exactJacExpr,:symbolic, symDict, dD )
(jac, exactJacExpr, dD) 

# output

(Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(1 => Set([1])), Dict{Expr, Union{Float64, Int64, Expr, Symbol}}(:((1, 1)) => :(-2.0 * p[2])), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(2 => Set([1])))
```

"""
function extractJacDepNormal(varNum::Int,rhs::Expr,jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},exactJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},jac_mode ::Symbol,symDict::Dict{Symbol,Expr},dD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 
  jacSet=Set{Union{Int,Symbol,Expr}}()
  #jacDiscrSet is not needed unlike  jacset. from jacset we construct jac and sd. for discrete variables there is not jac, there is only dD, so its construction will be direct.
  m=postwalk(rhs) do a   #
      if a isa Expr && a.head == :ref 
        if a.args[1]==:q#  e same as the continuous problem
          push!(jacSet,  (a.args[2]))  
          #a=eliminateRef(a)#q[i] -> qi
        elseif a.args[1]==:p# 
            dDset=Set{Union{Int,Symbol,Expr}}()
            if haskey(dD, (a.args[2]))    # dict dD already contains key a.args[2] (in this case var i) # optional check but to skip the get function
                dDset=get(dD,(a.args[2]),dDset) # if var di first time to influence some var, dDset is empty, otherwise get its set of influences 
            end
            push!(dDset,  varNum)      # p... also influences varNum, so update the set
            dD[(a.args[2])]=dDset       #update the dict
        end
        #a=symbolFromRef(a.args[1],a.args[2])
        symarg=symbolFromRef(a.args[1],a.args[2]) #after getting the index i in previous line, we can change the expression
        symDict[symarg]=a
        a=symarg
      end 
      return a 
  end
  if jac_mode==:symbolic
    # extract the jac (continuous part)
    basi = convert(Basic, m) # m ready: all refs are symbols
    for i in jacSet  # jacset contains vars in RHS
      symarg=symbolFromRef(:q,i) # specific to elements in jacSet: get q1 from 1 for exple
      coef = diff(basi, symarg) # symbolic differentiation: returns type Basic
      coefstr=string(coef);coefExpr=Meta.parse(coefstr)#convert from basic to expression
      jacEntry=restoreRef(coefExpr,symDict)# get back ref: qi->q[i]
      exactJacExpr[:(($varNum,$i))]=jacEntry # entry (varNum,i) is jacEntry
    end
  end
  if length(jacSet)>0 jac[varNum]=jacSet end
  #= if rhs isa Expr =# useNumber_inFuncalls(rhs) #= end =# # prevent using taylor0 in helper function calls because their body was not changed to zero-allocation ops.
end

# like above except (b,niter) instead of varNum
"""
    extractJacDepLoop(b::Int,niter::Int,rhs::Expr,jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},exactJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},jac_mode ::Symbol,symDict::Dict{Symbol,Expr},dD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 

This function is similar to the [`extractJacDepNormal`](@ref)  function by using the tuple (b,niter) instead of the integer varNum. Itextracts the jacobian dependency (jac) as well as the exacte symbolic jacobian expression (exactJacExpr) and the dependency of state derivatives to discrete variables (dD), in the form of dictionaries, from the differential equations that are written in a loop.\n

  For the continuous part, it sarts by looking for the 'i' in q[i] in the RHS and storing this 'i' in a jacSet. Then, it changes q[i] to qi for symbolic differentiation. After finding ``\\frac{\\partial f_i}{\\partial q_i} `` as the exact jacobian entry, it changes back qi to q[i]. Also, any mute variable from the differential equations is changed to 'i'.\n
  For the discrete part, the function puts the the tuple (b,niter) in a set, and stores this set in a dictionary dD with the key being the index of the discrete variable.
    # Example:
```jldoctest
using QuantizedSystemSolver
(b, niter, rhs, jac, exactJacExpr, symDict, dD) = (2, 9, :(p[1] * q[i - 1] * 1.5), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(1 => Set([1])), Dict{Expr, Union{Float64, Int64, Expr, Symbol}}(:((1, 1)) => :(-2.0 * p[2])), Dict{Symbol, Expr}(:q10 => :(q[10]), :p2 => :(p[2]), :qiminus1 => :(q[i - 1]), :p1 => :(p[1]), :q1 => :(q[1])), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(2 => Set([1])))
QuantizedSystemSolver.extractJacDepLoop(b, niter, rhs, jac, exactJacExpr,:symbolic, symDict, dD )
(jac, exactJacExpr, dD) 

# output

(Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(:((2, 9)) => Set([:(i - 1)]), 1 => Set([1])), Dict{Expr, Union{Float64, Int64, Expr, Symbol}}(:((1, 1)) => :(-2.0 * p[2]), :(((2, 9), i - 1)) => :(1.5 * p[1])), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(2 => Set([1]), 1 => Set([:((2, 9))])))
```
"""
function extractJacDepLoop(b::Int,niter::Int,rhs::Expr,jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},exactJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},jac_mode ::Symbol,symDict::Dict{Symbol,Expr},dD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 
  jacSet=Set{Union{Int,Symbol,Expr}}()
  m=postwalk(rhs) do a   
      if a isa Expr && a.head == :ref        
        if a.args[1]==:q# 
              push!(jacSet,  (a.args[2]))  #             
        elseif a.args[1]==:p
          if a.args[2] isa Int  
            dDset=Set{Union{Int,Symbol,Expr}}()
            if haskey(dD, (a.args[2]))
                dDset=get(dD,(a.args[2]),dDset)
            end
            push!(dDset,  :(($b,$niter)))
            dD[(a.args[2])]=dDset
          end
        end 
        #a=symbolFromRef(a.args[1],a.args[2])
        symarg=symbolFromRef(a.args[1],a.args[2]) #after getting the index i in previous line, we can change the expression
        symDict[symarg]=a
        a=symarg
      end
      return a
  end
  if jac_mode==:symbolic
    basi = convert(Basic, m)
    for i in jacSet
      symarg=symbolFromRef(:q,i);
      coef = diff(basi, symarg)
      coefstr=string(coef);
      coefExpr=Meta.parse(coefstr)
      jacEntry=restoreRef(coefExpr,symDict)
      exactJacExpr[:((($b,$niter),$i))]=jacEntry
    end
  end
  jac[:(($b,$niter))]=jacSet
  #= if rhs isa Expr =# useNumber_inFuncalls(rhs)#=  end  =#
end






"""
    createJacVect(jac:: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},::Val{T}) where {T}

  constructs the jacobian dependency as a vector from the existing dictionary jac resulted from [`extractJacDepNormal`](@ref) and [`extractJacDepLoop`](@ref) functions. \n

  This function just collects the data from the value of the dictionary if the key of the dictionary is an integer. (a dictionary contains(key=>value),...). In the case it is an expression :(b,niter), the function uses a 'for' loop to replace each 'b' by its corresponding integer. This approach depends on the size of the problem, but it runs one time.

# Example:   
```jldoctest
using QuantizedSystemSolver
jac=Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(1 => Set([2, 1]),:((2, 9)) => Set([:(i - 1), :i]),10 => Set([1, 10]))
jacVect=QuantizedSystemSolver.createJacVect(jac,Val(10) )
string(jacVect)

# output

"[[2, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [10, 1]]"

```

"""
function createJacVect(jac:: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},::Val{T}) where {T}
  jacVect = Vector{Vector{Int}}(undef, T)
  for i=1:T
      jacVect[i]=Vector{Int}()# define it so i can push elements as i find them below
  end
  for dictElement in jac
      if dictElement[1] isa Int  #jac[varNum]=jacSet
         jacVect[dictElement[1]]=[Int(eval(fa)) for fa in dictElement[2] ]# collect(dictElement[2]) # collect: set -> vector
      elseif dictElement[1] isa Expr  #jac[:(($b,$niter))]=jacSet
          for j_=(dictElement[1].args[1]):(dictElement[1].args[2])  
              temp=Vector{Int}()
              for element in dictElement[2]
                  if element isa Expr || element isa Symbol#can split symbol alone since no need to postwalk
                      fa= postwalk(a -> a isa Symbol && a==:i ? j_ : a, element) # change each symbol i to exact number 
                      push!(temp,Int(eval(fa))) # change i+1 for exmple to the exact number when i known
                  else #element is int
                      push!(temp,element)
                  end
              end
              jacVect[j_]=temp #fill the jac
          end
      end     
  end
  jacVect
end
"""
    createSDVect(jac:: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},::Val{T}) where {T}

 constructs the State to derivative dependency (opposite of jacobian dependency) as a vector from the existing dictionary jac resulted from the [`extractJacDepNormal`](@ref) and the [`extractJacDepLoop`](@ref) functions. It is the opposite in the sense that here we collect the keys into some vectors whereas in the jacobian dependency we collect the values of the dictionary in some vectors. \n

  If the key of the dictionary is an integer, then for all elements 'k' in the value of the dictionary (a set), the key is pushed into a new vector indexed at 'k'. In the case the key is an expression :(b,niter), the function uses a 'for' loop to replace each 'b' by its corresponding integer. This approach depends on the size of the problem, but it runs one time.
# Example:   
```jldoctest
using QuantizedSystemSolver
jac=Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(1 => Set([2, 1]),:((2, 9)) => Set([:(i - 1), :i]),10 => Set([1, 10]));
SD=QuantizedSystemSolver.createSDVect(jac,Val(10) );
string(SD)

# output

"[[10, 2, 1], [2, 3, 1], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9], [10]]"

```
"""
function createSDVect(jac:: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},::Val{T}) where {T}
  sdVect = Vector{Vector{Int}}(undef, T)
  for i=1:T
      sdVect[i]=Vector{Int}()# define it so I can push elements as I find them below
  end
  for dictElement in jac
      if dictElement[1] isa Int # key is an int
          for k in dictElement[2]   #elments values 
              push!(sdVect[Int(eval(k))],dictElement[1])
          end
      elseif dictElement[1] isa Expr # key is an expression
          for j_=(dictElement[1].args[1]):(dictElement[1].args[2])  
              for element in dictElement[2]
                  if element isa Expr || element isa Symbol#element=
                  fa= postwalk(a -> a isa Symbol && a==:i ? j_ : a, element)
                  push!(sdVect[Int(eval(fa))],j_)
                  else#element is int
                      push!(sdVect[element],j_)
                  end
              end
             
          end
      end     
  end
  sdVect
end

"""
    createExactJacFun(otherCode::Expr,Exactjac:: Dict{Expr,Union{Float64,Int,Symbol,Expr}},funName::Symbol,f) 

 constructs the exact jacobian entries as a function from the existing dictionary Exactjac resulted from extractJacDepNormal and extractJacDepLoop functions. \n

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
function createExactJacFun(otherCode::Expr,Exactjac:: Dict{Expr,Union{Float64,Int,Symbol,Expr}},funName::Symbol,f) 
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
  push!(otherCode.args,myex1)
  
  Base.remove_linenums!(otherCode)
  def1=Dict{Symbol,Any}() 
  def1[:head] = :function
  def1[:name] = Symbol(:exactJac,funName)  
  def1[:args] = [:(q::Vector{Taylor0}),:(p::Vector{Float64}),:(cache::AbstractVector{Float64}),:(i::Int),:(j::Int),:(t::Float64),:(f_)]
  def1[:body] = otherCode
  functioncode1=combinedef(def1)
end
 