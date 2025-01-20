
"""
    symbolFromRefdiscrete(refEx)

 Similar to [`symbolFromRef`](@ref)  but for discrete variables. It gets a symbol diplusNumber, diminusNumber, or ditimesNumber from expressions like i+Number, i+Number, i+Number.
"""
function symbolFromRefdiscrete(refEx)#refEx is i+1 in p[i+1] for example
  if refEx isa Expr #
    if refEx.args[1]==:+
      refEx=Symbol("p",(refEx.args[2]), "plus",(refEx.args[3]))
    elseif refEx.args[1]==:-
      refEx=Symbol("p",(refEx.args[2]), "minus",(refEx.args[3]))
    elseif refEx.args[1]==:*
      refEx=Symbol("p",(refEx.args[2]), "times",(refEx.args[3]))
    elseif refEx.args[1]==:/
      Error("parameter_index division is not implemented Yet")
    end
  else
    refEx=Symbol("p",(refEx))
  end
  return refEx
end

"""
    EventDependencyStruct
 A struct that holds the event dependency information. It has the following fields:
  - `id::Int:` the id of the event
  - `evCont::Vector{Int}:` the index tracking used for HD & HZ. Also it is used to update q,quantum,recomputeNext when x is modified in an event
  - `evDisc::Vector{Int}:` the index tracking used for HD & HZ.
  - `evContRHS::Vector{Int}:` the index tracking used to update other Qs before executing the event

"""
struct EventDependencyStruct
  id::Int
  evCont::Vector{Int} #index tracking used for HD & HZ. Also it is used to update q,quantum,recomputeNext when x is modified in an event
  evDisc::Vector{Int} #index tracking used for HD & HZ.
  evContRHS::Vector{Int} #index tracking used to update other Qs before executing the event
end
# the following functions handle discrete problems
"""
    extractJacDepNormalDiscrete(varNum::Int,rhs::Union{Symbol,Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},exactJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},symDict::Dict{Symbol,Expr},dD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 

Extracts the jacobian dependency (jac) as well as the exacte symbolic jacobian expression (exactJacExpr) and the dependency of state derivatives to discrete variables (dD), in the form of dictionaries, from the simple differential equations.\n

  For the continuous part, similar to [`extractJacDepNormal`](@ref)  function, this function sarts by looking for the 'i' in q[i] in the RHS and storing this 'i' in a jacSet for the varNum. Then, it changes q[i] to qi for symbolic differentiation. After finding ``\\frac{\\partial f_i}{\\partial q_i} `` as the exact jacobian entry, it changes back qi to q[i]. Also, any mute variable from the differential equations is changed to 'i'.\n
  For the discrete part, the function puts the index of the differential equation in a set, and stores this set in a dictionary dD with the key being the index of the discrete variable.
  # Example:
```jldoctest
using QuantizedSystemSolver
(varNum, rhs, jac, exactJacExpr, symDict, dD) = (1, :(p[2] - 2.0 * q[1] * p[2]), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(), Dict{Expr, Union{Float64, Int64, Expr, Symbol}}(), Dict{Symbol, Expr}(:q10 => :(q[10]), :p2 => :(p[2]), :qiminus1 => :(q[i - 1]), :p1 => :(p[1]), :q1 => :(q[1])), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}())
QuantizedSystemSolver.extractJacDepNormalDiscrete(varNum, rhs, jac, exactJacExpr, symDict, dD )
(jac, exactJacExpr, dD) 

# output

(Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(1 => Set([1])), Dict{Expr, Union{Float64, Int64, Expr, Symbol}}(:((1, 1)) => :(-2.0 * p[2])), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(2 => Set([1])))
```

"""
function extractJacDepNormalDiscrete(varNum::Int,rhs::Union{Symbol,Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},exactJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},symDict::Dict{Symbol,Expr},dD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 
  jacSet=Set{Union{Int,Symbol,Expr}}()
  #jacDiscrSet is not needed unlike  jacset. from jacset we construct jac and sd. for discrete variables there is not jac, there is only dD, so its construction will be direct.
  m=postwalk(rhs) do a   #
      if a isa Expr && a.head == :ref && a.args[1]==:q#  these 3 lines are the same as the continuous problem
          push!(jacSet,  (a.args[2]))  
          a=eliminateRef(a)#q[i] -> qi
      elseif a isa Expr && a.head == :ref && a.args[1]==:p# 
          dDset=Set{Union{Int,Symbol,Expr}}()
          if haskey(dD, (a.args[2]))    # dict dD already contains key a.args[2] (in this case var i) # optional check but to skip the get function
              dDset=get(dD,(a.args[2]),dDset) # if var di first time to influence some var, dDset is empty, otherwise get its set of influences 
          end
          push!(dDset,  varNum)      # p... also influences varNum, so update the set
          dD[(a.args[2])]=dDset       #update the dict
          a=eliminateRef(a)#p[i] -> di
      end
      return a 
  end
  # extract the jac (continuous part)
  basi = convert(Basic, m) # m ready: all refs are symbols
  for i in jacSet  # jacset contains vars in RHS
    symarg=symbolFromRef(i) # specific to elements in jacSet: get q1 from 1 for exple
    coef = diff(basi, symarg) # symbolic differentiation: returns type Basic
    coefstr=string(coef);coefExpr=Meta.parse(coefstr)#convert from basic to expression
    jacEntry=restoreRef(coefExpr,symDict)# get back ref: qi->q[i]
    exactJacExpr[:(($varNum,$i))]=jacEntry # entry (varNum,i) is jacEntry
  end
  if length(jacSet)>0 jac[varNum]=jacSet end
end

# like above except (b,niter) instead of varNum
"""
    extractJacDepLoopDiscrete(b::Int,niter::Int,rhs::Union{Symbol,Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},exactJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},symDict::Dict{Symbol,Expr},dD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 

This function is similar to the [`extractJacDepNormalDiscrete`](@ref)  function by using the tuple (b,niter) instead of the integer varNum. Itextracts the jacobian dependency (jac) as well as the exacte symbolic jacobian expression (exactJacExpr) and the dependency of state derivatives to discrete variables (dD), in the form of dictionaries, from the differential equations that are written in a loop.\n

  For the continuous part, it sarts by looking for the 'i' in q[i] in the RHS and storing this 'i' in a jacSet. Then, it changes q[i] to qi for symbolic differentiation. After finding ``\\frac{\\partial f_i}{\\partial q_i} `` as the exact jacobian entry, it changes back qi to q[i]. Also, any mute variable from the differential equations is changed to 'i'.\n
  For the discrete part, the function puts the the tuple (b,niter) in a set, and stores this set in a dictionary dD with the key being the index of the discrete variable.
    # Example:
```jldoctest
using QuantizedSystemSolver
(b, niter, rhs, jac, exactJacExpr, symDict, dD) = (2, 9, :(p[1] * q[i - 1] * 1.5), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(1 => Set([1])), Dict{Expr, Union{Float64, Int64, Expr, Symbol}}(:((1, 1)) => :(-2.0 * p[2])), Dict{Symbol, Expr}(:q10 => :(q[10]), :p2 => :(p[2]), :qiminus1 => :(q[i - 1]), :p1 => :(p[1]), :q1 => :(q[1])), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(2 => Set([1])))
QuantizedSystemSolver.extractJacDepLoopDiscrete(b, niter, rhs, jac, exactJacExpr, symDict, dD )
(jac, exactJacExpr, dD) 

# output

(Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(:((2, 9)) => Set([:(i - 1)]), 1 => Set([1])), Dict{Expr, Union{Float64, Int64, Expr, Symbol}}(:((1, 1)) => :(-2.0 * p[2]), :(((2, 9), i - 1)) => :(1.5 * p[1])), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(2 => Set([1]), 1 => Set([:((2, 9))])))
```
"""
function extractJacDepLoopDiscrete(b::Int,niter::Int,rhs::Union{Symbol,Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},exactJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},symDict::Dict{Symbol,Expr},dD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 
  jacSet=Set{Union{Int,Symbol,Expr}}()
  m=postwalk(rhs) do a   
      if a isa Expr && a.head == :ref && a.args[1]==:q# 
              push!(jacSet,  (a.args[2]))  #
              a=eliminateRef(a)#q[i] -> qi
      elseif a isa Expr && a.head == :ref && a.args[1]==:p
        if a.args[2] isa Int  
          dDset=Set{Union{Int,Symbol,Expr}}()
          if haskey(dD, (a.args[2]))
              dDset=get(dD,(a.args[2]),dDset)
          end
          push!(dDset,  :(($b,$niter)))
          dD[(a.args[2])]=dDset
        end 
        a=eliminateRef(a)#q[i] -> qi
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
    exactJacExpr[:((($b,$niter),$i))]=jacEntry
  end
  jac[:(($b,$niter))]=jacSet
end
"""
    extractZCJacDepNormal(counter::Int,zcf::Expr,zcjac :: Vector{Vector{Int}},SZ ::Dict{Int,Set{Int}},dZ :: Dict{Int,Set{Int}}) 

Extracts the zero-crossing jacobian dependency as a vector (zcjac), the dependency of the zero-crossing functions to continuous (SZ) and discrete variables (dZ) in the form of dictionaries, from the 'if-statements' (zcf).\n

 The zcjac is a vector of vectors, where each vector contains the indices of the continuous variables that the zero-crossing function depends on. The SZ dictionary contains the indices of the zero-crossing functions as values and the indices of the continuous variables as keys. Similarly, the dZ dictionary contains the indices of the zero-crossing functions as values and the indices of the discrete variables as keys.
   # Example:
```jldoctest
using QuantizedSystemSolver
(counter, zcf, zcjac, SZ, dZ) = (2, :(q[2] - p[1]), [[1]], Dict{Int64, Set{Int64}}(1 => Set([1])), Dict{Int64, Set{Int64}}())
QuantizedSystemSolver.extractZCJacDepNormal(counter, zcf, zcjac, SZ, dZ)
(zcjac, SZ, dZ) 

# output

([[1], [2]], Dict{Int64, Set{Int64}}(2 => Set([2]), 1 => Set([1])), Dict{Int64, Set{Int64}}(1 => Set([2])))
```
"""
function extractZCJacDepNormal(counter::Int,zcf::Expr,zcjac :: Vector{Vector{Int}},SZ ::Dict{Int,Set{Int}},dZ :: Dict{Int,Set{Int}}) 
  zcjacSet=Set{Int}()
 #zcjacDiscrSet not needed unlike  zcjacSet. from zcjacSet we construct zcjac and SZ. for discrete variables there is not zcjac, there is only dZ, so its construction will be direct.
  postwalk(zcf) do a   #
      if a isa Expr && a.head == :ref && a.args[1]==:q# 
          push!(zcjacSet,  (a.args[2]))  #
          SZset=Set{Int}()
          if haskey(SZ, (a.args[2]))
              SZset=get(SZ,(a.args[2]),SZset)
          end
          push!(SZset,  counter)
          SZ[(a.args[2])]=SZset
      elseif a isa Expr && a.head == :ref && a.args[1]==:p# 
          dZset=Set{Int}()
          if haskey(dZ, (a.args[2]))
              dZset=get(dZ,(a.args[2]),dZset)
          end
          push!(dZset,  counter)
          dZ[(a.args[2])]=dZset
      end
      return a 
  end
  push!(zcjac,collect(zcjacSet))#convert set to vector
end


"""
    createSZvect(SZ :: Dict{Int64, Set{Int64}},::Val{T}) where{T} 

constructs the zero-crossing dependency to state variables as a vector from the existing dictionary SZ resulted from the [`extractZCJacDepNormal`](@ref) function. The continuous variables are the keys and the zero-crossing are the values.\n
   # Example:
```jldoctest
using QuantizedSystemSolver
(SZ, T) = (Dict{Int64, Set{Int64}}(2 => Set([2]), 1 => Set([1])), 10)
szVec=QuantizedSystemSolver.createSZvect(SZ, Val(T))
string(szVec)

# output

"[[1], [2], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[]]"
```
"""
function createSZvect(SZ :: Dict{Int64, Set{Int64}},::Val{T}) where{T}
  szVect = Vector{Vector{Int}}(undef, T)
  for ii=1:T
      szVect[ii]=Vector{Int}()# define it so i can push elements as i find them below
  end
  for ele in SZ
      szVect[ele[1]]=collect(ele[2])
  end
  szVect
end

"""
    createdDvect(dD::Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}},::Val{D}) where{D}

constructs the State to derivative dependency to discrete variables as a vector from the existing dictionary dD resulted from the [`extractJacDepNormalDiscrete`](@ref) and the [`extractJacDepLoopDiscrete`](@ref) functions. The discrete variables are the keys and the differential equations are the values. This dependency is needed only in [`createDependencyToEventsDiscr`](@ref)  .\n
   # Example:
```jldoctest
using QuantizedSystemSolver
(dD, D) = (Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(2 => Set([10, 1]), 1 => Set([:((2, 9))])), 2)
dDVect=QuantizedSystemSolver.createdDvect(dD, Val(D) )
string(dDVect)

# output

"[[2, 3, 4, 5, 6, 7, 8, 9], [10, 1]]"
```
"""
function createdDvect(dD::Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}},::Val{D}) where{D}
  dDVect = Vector{Vector{Int}}(undef, D)
  for ii=1:D
      dDVect[ii]=Vector{Int}()# define it so i can push elements as i find them below
  end
  for dictElement in dD
      temp=Vector{Int}()
      for k in dictElement[2] 
          if k isa Int
              push!(temp,k)
          elseif k isa Expr #tuple
              for j_=(k.args[1]):(k.args[2])  
                  push!(temp,j_)
              end
          end
      end
     dDVect[dictElement[1]]=temp
  end
  dDVect
end

"""
    createDependencyToEventsDiscr(dD::Vector{Vector{Int}},dZ::Dict{Int64, Set{Int64}},eventDep::Vector{EventDependencyStruct}) 

constructs the dependency of zero-crossing functions and state derivatives to events using only discrete variables.
# arguments
- `dD::Vector{Vector{Int}}:` the dependency of state derivatives to discrete variables as a vector
- `dZ::Dict{Int64, Set{Int64}}:` the dependency of zero-crossing functions to discrete variables as a dictionary
- `eventDep::Vector{EventDependencyStruct}:` the event dependency information as a vector of structs
# returns
- `HZ1::Vector{Vector{Int}}:` the dependency of zero-crossing functions to events using only discrete variables
- `HD1::Vector{Vector{Int}}:` the dependency of differential equations to events using only discrete variables

   # Example:
```jldoctest
using QuantizedSystemSolver
(dD, dZ, eventDep) = ([[2, 3, 4, 5, 6, 7, 8, 9], [10, 1]], Dict{Int64, Set{Int64}}(1 => Set([2])), QuantizedSystemSolver.EventDependencyStruct[QuantizedSystemSolver.EventDependencyStruct(1, Int64[], [1], Int64[]), QuantizedSystemSolver.EventDependencyStruct(2, Int64[], Int64[], Int64[]), QuantizedSystemSolver.EventDependencyStruct(3, [3], [2], [3, 1, 2]), QuantizedSystemSolver.EventDependencyStruct(4, Int64[], Int64[], Int64[])])
(HZ1, HD1) =QuantizedSystemSolver.createDependencyToEventsDiscr(dD, dZ, eventDep )
(HZ1, HD1) 

# output

([[2], Int64[], Int64[], Int64[]], [[5, 4, 6, 7, 2, 9, 8, 3], Int64[], [10, 1], Int64[]])
```
"""
function createDependencyToEventsDiscr(dD::Vector{Vector{Int}},dZ::Dict{Int64, Set{Int64}},eventDep::Vector{EventDependencyStruct}) 
    Y=length(eventDep) # number of events
    lendD=length(dD) # number of discrete variables
    HD1 = Vector{Vector{Int}}(undef, Y) # dependency of differential equations to events using only discrete variables
    HZ1 = Vector{Vector{Int}}(undef, Y) # dependency of zero-crossing functions to events using only discrete variables
      for ii=1:Y
        HD1[ii] =Vector{Int}()# define it so i can push elements as i find them below
        HZ1[ii] =Vector{Int}()# define it so i can push elements as i find them below
      end
    for j=1:Y #for each event j
      hdSet=Set{Int}()
      hzSet=Set{Int}()
      evdiscrete=eventDep[j].evDisc # vector of discrete variables affected by event j
      for i in evdiscrete # i is the index of the discrete variable
          if i<=lendD  
            for k in dD[i] # k is the index of the differential equation affected by discrete variable i
              push!(hdSet,k)#hdSet is the set of differential equations affected by all discrete variables in evdiscrete
            end
          end
          tempSet=Set{Int}()# index of zero-crossing functions affected discrete variable i
          if haskey(dZ, i)
            tempSet=get(dZ,i,tempSet)
          end
          for kk in tempSet
            push!(hzSet,kk)#hzset is the set of zero-crossing functions affected by all discrete variables in evdiscrete
          end
      end
      HD1[j] =collect(hdSet)
      HZ1[j] =collect(hzSet)
    end #end for j  (events)
    return (HZ1,HD1)
end 
"""
    createDependencyToEventsCont(SD::Vector{Vector{Int}},sZ::Dict{Int64, Set{Int64}},eventDep::Vector{EventDependencyStruct}) 

 
constructs the dependency of zero-crossing functions and state derivatives to events using only continuous variables.
# arguments
- `SD::Vector{Vector{Int}}:` the dependency of state derivatives to continuous variables as a vector
- `sZ::Dict{Int64, Set{Int64}}:` the dependency of zero-crossing functions to continuous variables as a dictionary
- `eventDep::Vector{EventDependencyStruct}:` the event dependency information as a vector of structs
# returns
- `HZ2::Vector{Vector{Int}}:` the dependency of zero-crossing functions to events using only continuous variables
- `HD2::Vector{Vector{Int}}:` the dependency of differential equations to events using only continuous variables

   # Example:
```jldoctest
using QuantizedSystemSolver
(SD, sZ, eventDep) = ([[10, 2, 1], [3], [4], [5], [6], [7], [8], [9], Int64[], [10]], Dict{Int64, Set{Int64}}(2 => Set([2]), 1 => Set([1])), QuantizedSystemSolver.EventDependencyStruct[QuantizedSystemSolver.EventDependencyStruct(1, Int64[], [1], Int64[]), QuantizedSystemSolver.EventDependencyStruct(2, Int64[], Int64[], Int64[]), QuantizedSystemSolver.EventDependencyStruct(3, [3], [2], [3, 1, 2]), QuantizedSystemSolver.EventDependencyStruct(4, Int64[], Int64[], Int64[])])
(HZ2, HD2) =QuantizedSystemSolver.createDependencyToEventsCont(SD, sZ, eventDep)
(HZ2, HD2) 

# output

([Int64[], Int64[], Int64[], Int64[]], [Int64[], Int64[], [4], Int64[]])
```
"""
function createDependencyToEventsCont(SD::Vector{Vector{Int}},sZ::Dict{Int64, Set{Int64}},eventDep::Vector{EventDependencyStruct}) 
  Y=length(eventDep)
  HD2 = Vector{Vector{Int}}(undef, Y)
  HZ2 = Vector{Vector{Int}}(undef, Y)
    for ii=1:Y
      HD2[ii] =Vector{Int}()# define it so i can push elements as i find them below
      HZ2[ii] =Vector{Int}()# define it so i can push elements as i find them below
    end
  for j=1:Y
    hdSet=Set{Int}()
    hzSet=Set{Int}()
      evContin=eventDep[j].evCont#vector of continuous variables affected by event j
      for i in evContin # i is the index of the continuous variable
            for k in SD[i] # k is the index of the differential equation affected by continuous variable i
              push!(hdSet,k)#hdSet is the set of differential equations affected by all continuous variables in evContin
            end
            tempSet=Set{Int}()# index of zero-crossing functions affected continuous variable i
            if haskey(sZ, i)
              tempSet=get(sZ,i,tempSet)
            end
            for kk in tempSet
              push!(hzSet,kk)#hzset is the set of zero-crossing functions affected by all continuous variables in evContin
            end
      end
      HD2[j] =collect(hdSet)# 
      HZ2[j] =collect(hzSet)
  end #end for j  (events)
  return (HZ2,HD2)
end 

"""
    unionDependency(HZD1::Vector{Vector{Int}},HZD2::Vector{Vector{Int}})
 merges the state derivatives and zero-crossing functions dependencies to events using both continuous and discrete variables.
 
   # Example:
```jldoctest
using QuantizedSystemSolver
(HD1, HD2) = ([[5, 4, 6, 7, 2, 9, 8, 3], Int64[], [10, 1], Int64[]], [Int64[], Int64[], [4], Int64[]])
HD=QuantizedSystemSolver.unionDependency(HD1, HD2)
string(HD)

# output

"[[5, 4, 6, 7, 2, 9, 8, 3], Int64[], [4, 10, 1], Int64[]]"
```
"""
function unionDependency(HZD1::Vector{Vector{Int}},HZD2::Vector{Vector{Int}})
  Y=length(HZD1)
  HZD = Vector{Vector{Int}}(undef, Y)
    for ii=1:Y
      HZD[ii] =Vector{Int}()# define it so i can push elements as i find them below
    end
  for j=1:Y
    hzSet=Set{Int}()
    for kk in HZD1[j]
      push!(hzSet,kk)
    end
    for kk in HZD2[j]
      push!(hzSet,kk)
    end
    HZD[j]=collect(hzSet)
  end
  HZD
end 


function createDiscEqFun(otherCode::Expr,equs::Dict{Union{Int,Expr},Union{Int,Symbol,Expr}},zcequs::Vector{Expr},eventequs::Vector{Expr},fname::Symbol,f::F) where{F}
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
  push!(otherCode.args,myex1)
   ##############ZCF###################
  if length(zcequs)>0
      s="if zc==1  $(zcequs[1]) ;return nothing"
      for i=2:length(zcequs)
          s*= " elseif zc==$i $(zcequs[i]) ;return nothing"
      end
      s*= " end "
      myex2=Meta.parse(s)
      push!(otherCode.args,myex2)
  end
   #############events#################
  if length(eventequs)>0
      s= "if ev==1  $(eventequs[1]) ;return nothing"
      for i=2:length(eventequs)
          s*= " elseif ev==$i $(eventequs[i]) ;return nothing"
      end
      s*= " end "
      myex3=Meta.parse(s)
      push!(otherCode.args,myex3)
  end
 
  Base.remove_linenums!(otherCode)
  def=Dict{Symbol,Any}()
  def[:head] = :function
  def[:name] = fname  
  def[:args] = [:(i::Int),:(zc::Int),:(ev::Int),:(q::Vector{Taylor0}),:(p::Vector{Float64}), :(t::Taylor0),:(cache::Vector{Taylor0}),:(f_::F)]
  def[:body] = otherCode 
  functioncode=combinedef(def)
 # @show functioncode

end