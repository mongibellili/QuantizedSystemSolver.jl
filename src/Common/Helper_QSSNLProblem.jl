function changeExprToFirstValue(ex::Expr)##
  newEx=postwalk(ex) do a  # change u[1] to u[1][0]
      if a isa Expr && a.head == :ref && a.args[1]==:q # change q[1] to q[1][0]
           outerRef=Expr(:ref)
          push!(outerRef.args,a)
          push!(outerRef.args,:(0))
          a=outerRef
      end
      return a
  end
  newEx
end
function eliminateRef(a)#q[i] -> qi
  if a.args[2] isa Expr 
    if a.args[2].args[1]==:+
      a=Symbol((a.args[1]),(a.args[2].args[2]), "plus",(a.args[2].args[3]))
    elseif a.args[2].args[1]==:-
      a=Symbol((a.args[1]),(a.args[2].args[2]), "minus",(a.args[2].args[3]))
    elseif a.args[2].args[1]==:*
      a=Symbol((a.args[1]),(a.args[2].args[2]), "times",(a.args[2].args[3]))
    elseif a.args[2].args[1]==:/
      a=Symbol((a.args[1]),(a.args[2].args[2]), "over",(a.args[2].args[3]))
    end
  else
    a=Symbol((a.args[1]),(a.args[2]))
  end
  return a
end
function symbolFromRef(refEx)#refEx is i+1 in q[i+1] for example
  if refEx isa Expr #
    if refEx.args[1]==:+
      refEx=Symbol("q",(refEx.args[2]), "plus",(refEx.args[3]))
    elseif refEx.args[1]==:-
      refEx=Symbol("q",(refEx.args[2]), "minus",(refEx.args[3]))
    elseif refEx.args[1]==:*
      refEx=Symbol("q",(refEx.args[2]), "times",(refEx.args[3]))
    elseif refEx.args[1]==:/
      refEx=Symbol("q",(refEx.args[2]), "over",(refEx.args[3]))
    end
  else
    refEx=Symbol("q",(refEx))
  end
  return refEx
end
function symbolFromRefd(refEx)#refEx is i+1 in q[i+1] for example
  if refEx isa Expr #
    if refEx.args[1]==:+
      refEx=Symbol("d",(refEx.args[2]), "plus",(refEx.args[3]))
    elseif refEx.args[1]==:-
      refEx=Symbol("d",(refEx.args[2]), "minus",(refEx.args[3]))
    elseif refEx.args[1]==:*
      refEx=Symbol("d",(refEx.args[2]), "times",(refEx.args[3]))
    elseif refEx.args[1]==:/
      refEx=Symbol("d",(refEx.args[2]), "over",(refEx.args[3]))
    end
  else
    refEx=Symbol("d",(refEx))
  end
  return refEx
end
function restoreRef(coefExpr,symDict)
  newEx=postwalk(coefExpr) do element# 
    if element isa Symbol && !(element in (:+,:-,:*,:/)) && haskey(symDict, element) && element != :d 
      element=symDict[element]
      element=changeExprToFirstValue(element)# change u[1] to u[1][0]
    elseif element== :d 
      element=symDict[element]
    end
    return element
  end#end postwalk
  newEx
end
function changeVarNames_params(ex::Expr,stateVarName::Symbol,muteVar::Symbol,param::Dict{Symbol,Union{Float64,Expr}},symDict::Dict{Symbol,Expr})#
  newEx=postwalk(ex) do element#postwalk to change var names and parameters
      if element isa Symbol   
          if haskey(param, element)#symbol is a parameter
              element=copy(param[element]) # copy needed in the case symbol id=expression substitued in equations...do not want all eqs reference same expression...ie if 1 eq changes, other eqs change
          elseif element==stateVarName #symbol is a var
              element=:q 
          elseif element==:discrete #symbol is a discr var
              element=:d
          elseif element==muteVar #symbol is a mute var
              element=:i
          #= else  # + - * /
               =#
          end
      elseif element isa Expr && element.head == :ref && element.args[1]==:q# 
            symarg=symbolFromRef(element.args[2])  #q[i] -> qi
            symDict[symarg]=element #store this translation  q[i] <-> qi for later use
      elseif element isa Expr && element.head == :ref && element.args[1]==:d#   
        symarg=symbolFromRefd(element.args[2])  #d[i] -> di
        symDict[symarg]=element #store this translation  d[i] <-> di   
      end
      return element
    end#end postwalk
  newEx
end
function changeVarNames_params(ex::Expr,stateVarName::Symbol,muteVar::Symbol,param::Dict{Symbol,Union{Float64,Expr}})######special for if statements and events
  newEx=postwalk(ex) do element#postwalk to change var names and parameters
      if element isa Symbol   
          if haskey(param, element)#symbol is a parameter
              element=copy(param[element])
          elseif element==stateVarName #symbol is a var
              element=:q 
          elseif element==:discrete #symbol is a discr var
              element=:d
          elseif element==muteVar #symbol is a mute var
              element=:i
          end
      end
      return element
    end#end postwalk
  newEx
end

# these 2 function handle continuous problems only
function extractJacDepNormal(varNum::Int,rhs::Union{Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}, exacteJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},symDict::Dict{Symbol,Expr}) 
  jacSet=Set{Union{Int,Symbol,Expr}}()
  m=postwalk(rhs) do a   #
      if a isa Expr && a.head == :ref # 
              push!(jacSet,  (a.args[2]))  # du[varNum=1]=rhs=u[5]+u[2] : 2 and 5 are stored in jacset
              a=eliminateRef(a)#q[i] -> qi
      end
      return a 
  end
  basi = convert(Basic, m)
  for i in jacSet
    symarg=symbolFromRef(i) # specific to elements in jacSet: get q1 from 1 for exple
    coef = diff(basi, symarg) # symbolic differentiation: returns type Basic
    coefstr=string(coef);coefExpr=Meta.parse(coefstr)#convert from basic to expression
    jacEntry=restoreRef(coefExpr,symDict)# get back ref: qi->q[i][0]  ...0 because later in exactJac fun cache[1]::Float64=jacEntry
    exacteJacExpr[:(($varNum,$i))]=jacEntry # entry (varNum,i) is jacEntry
  end
  if length(jacSet)>0 jac[varNum]=jacSet end # jac={1->(2,5)}
  #@show jac
end

function extractJacDepLoop(b::Int,niter::Int,rhs::Union{Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}} ,exacteJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},symDict::Dict{Symbol,Expr}) 
  jacSet=Set{Union{Int,Symbol,Expr}}()
  m=postwalk(rhs) do a   
      if a isa Expr && a.head == :ref # 
              push!(jacSet,  (a.args[2]))  #
              a=eliminateRef(a)
      end
      return a
  end
  basi = convert(Basic, m)
  for i in jacSet
    symarg=symbolFromRef(i);
    coef = diff(basi, symarg)
    coefstr=string(coef);
    coefExpr=Meta.parse(coefstr)
    jacEntry=restoreRef(coefExpr,symDict)
    exacteJacExpr[:((($b,$niter),$i))]=jacEntry
  end
  if length(jacSet)>0 jac[:(($b,$niter))]=jacSet end
end