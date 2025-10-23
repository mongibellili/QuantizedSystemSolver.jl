

"""
    extractJacDep(b::Int,niter::Int,rhs::Expr,jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},exactJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},jac_mode ::Symbol,symDict::Dict{Symbol,Expr},dD :: Dict{Union{Int,Symbol,Expr},Set{Union{Int,Symbol,Expr}}}) 

Extracts the jacobian dependency (jac) as well as the exacte symbolic jacobian expression (exactJacExpr) and the dependency of state derivatives to discrete variables (dD), in the form of dictionaries, from the simple differential equations.\n

  For the continuous part, similar to [`extractJacDep`](@ref)  function, this function sarts by looking for the 'i' in q[i] in the RHS and storing this 'i' in a jacSet for the varNum. Then, it changes q[i] to qi for symbolic differentiation. After finding ``\\frac{\\partial f_i}{\\partial q_i} `` as the exact jacobian entry, it changes back qi to q[i]. Also, any mute variable from the differential equations is changed to 'i'.\n
  For the discrete part, the function puts the index of the differential equation in a set, and stores this set in a dictionary dD with the key being the index of the discrete variable.
  # Example:
```jldoctest
using QuantizedSystemSolver
(varNum, rhs, jac, exactJacExpr, symDict, dD) = (1, :(p[2] - 2.0 * q[1] * p[2]), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(), Dict{Expr, Union{Float64, Int64, Expr, Symbol}}(), Dict{Symbol, Expr}(:q10 => :(q[10]), :p2 => :(p[2]), :qiminus1 => :(q[i - 1]), :p1 => :(p[1]), :q1 => :(q[1])), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}())
QuantizedSystemSolver.extractJacDep(varNum, rhs, jac, exactJacExpr,:symbolic, symDict, dD )
(jac, exactJacExpr, dD) 

# output

(Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(1 => Set([1])), Dict{Expr, Union{Float64, Int64, Expr, Symbol}}(:((1, 1)) => :(-2.0 * p[2])), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(2 => Set([1])))
```

"""
function extractJacDep(b::Int,niter::Int,rhs::Expr,jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},dD :: Dict{Union{Int,Symbol,Expr},Set{Union{Int,Symbol,Expr}}}) 
  
  #jacDiscrSet is not needed unlike  jacset. from jacset we construct jac and sd. for discrete variables there is not jac, there is only dD, so its construction will be direct.
  #postwalk

  jacSet=Set{Union{Int,Symbol,Expr}}()
  postwalk(rhs) do a   #
      if a isa Expr && a.head == :ref 
        base, idx1 =  a.args[1], a.args[2]  #parse_ref_upto2(a)
        if base==:q#  
          push!(jacSet,  idx1)  
        elseif base==:p# 
            dDset=Set{Union{Int,Symbol,Expr}}() # dDset cannot contain a symbol either a normal equation (int) or loop (expr)-------> delete symbol later
            if haskey(dD, idx1)    # dict dD already contains key a.args[2] (in this case var i) # optional check but to skip the get function
                dDset=get(dD,idx1,dDset) # if var di first time to influence some var, dDset is empty, otherwise get its set of influences 
            end
            if b!=-1
              push!(dDset,  :(($b,$niter)))
            else
              push!(dDset,  niter)  
            end
            dD[idx1]=dDset       #update the dict
        end
      
      end 
      return a 
  end

  
  if length(jacSet)>0
    if b!=-1
      jac[:(($b,$niter))]=jacSet
    else
      jac[niter]=jacSet
    end
  end

  jacSet

   

   #useNumber_inFuncalls(rhs)  # prevent using taylor0 in helper function calls because their body was not changed to zero-allocation ops.
   #rhs
end

function extractJacExpression(b::Int,niter::Int,rhs::Expr,jacSet::Set{Union{Int,Symbol,Expr}},exactJacExpr :: Dict{Expr,ScopedEquation},helperAssignments::Vector{AbstractODEStatement},symDict::Dict{Symbol,Expr}) 
  
 
    empty!(symDict) # in each equation fk, we change the symbols, find the dfk/dxi where i ∈ jacset

    m=postwalk(rhs) do a   #
      if a isa Expr && a.head == :ref 
        base, idx1 =  a.args[1], a.args[2] 
        symarg=Symbol(string(expr_to_flat_name(base), _idxToString(idx1)))
        symDict[symarg]=a
        a=symarg
      elseif a isa Expr && a.head == :vect 
         symarg=Symbol(expr_to_flat_name(a)) 
        symDict[symarg]=a
        a=symarg
      end 
      return a 
  end



    # extract the jac (continuous part)
    basi = convert(Basic, m) # m ready: all refs are symbols
    for i in jacSet  # jacset contains vars in RHS
      symarg=Symbol(string(:q, _idxToString(i)))#symbolFromRef(:q,i) # specific to elements in jacSet: get q1 from 1 for exple
      coef = diff(basi, symarg) # symbolic differentiation: returns type Basic
      coefstr=string(coef);
      coefExpr=Meta.parse(coefstr)#convert from basic to expression
      jacEntry=restoreRef(coefExpr,symDict)# get back ref: 
      jacEntry = changeExprToFirstValue(jacEntry)  # qi->q[i][0]
      jacEntry=changeT(jacEntry) # t -> t[0]
      if b!=-1
        exactJacExpr[:((($b,$niter),$i))]= ScopedEquation(helperAssignments,jacEntry )
      else
        exactJacExpr[:(($niter,$i))]=ScopedEquation(helperAssignments,jacEntry )
      end
    end
 

   #useNumber_inFuncalls(rhs)  # prevent using taylor0 in helper function calls because their body was not changed to zero-allocation ops.
   #rhs
end

"""
    createJacVect(jac:: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},::Val{T}) where {T}

  constructs the jacobian dependency as a vector from the existing dictionary jac resulted from [`extractJacDep`](@ref) and [`extractJacDepLoop`](@ref) functions. \n

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

 constructs the State to derivative dependency (opposite of jacobian dependency) as a vector from the existing dictionary jac resulted from the [`extractJacDep`](@ref) and the [`extractJacDepLoop`](@ref) functions. It is the opposite in the sense that here we collect the keys into some vectors whereas in the jacobian dependency we collect the values of the dictionary in some vectors. \n

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
    extractZCJacDep(counter::Int,zcf::Expr,zcjac :: Vector{Vector{Int}},SZ ::Dict{Int,Set{Int}},dZ :: Dict{Int,Set{Int}}) 

Extracts the zero-crossing jacobian dependency as a vector (zcjac), the dependency of the zero-crossing functions to continuous (SZ) and discrete variables (dZ) in the form of dictionaries, from the 'if-statements' (zcf).\n

 The zcjac is a vector of vectors, where each vector contains the indices of the continuous variables that the zero-crossing function depends on. The SZ dictionary contains the indices of the zero-crossing functions as values and the indices of the continuous variables as keys. Similarly, the dZ dictionary contains the indices of the zero-crossing functions as values and the indices of the discrete variables as keys.
   # Example:
```jldoctest
using QuantizedSystemSolver
(counter, zcf, zcjac, SZ, dZ) = (2, :(q[2] - p[1]), [[1]], Dict{Int64, Set{Int64}}(1 => Set([1])), Dict{Int64, Set{Int64}}())
QuantizedSystemSolver.extractZCJacDep(counter, zcf, zcjac, SZ, dZ)
(zcjac, SZ, dZ) 

# output

([[1], [2]], Dict{Int64, Set{Int64}}(2 => Set([2]), 1 => Set([1])), Dict{Int64, Set{Int64}}(1 => Set([2])))
```
"""
function extractZCJacDep(counter::Int,zcf::Expr,zcjac :: Vector{Vector{Int}},SZ ::Dict{Int,Set{Int}},dZ :: Dict{Int,Set{Int}}) 
  zcjacSet=Set{Int}()
 # @show zcf
 #zcjacDiscrSet not needed unlike  zcjacSet. from zcjacSet we construct zcjac and SZ. for discrete variables there is not zcjac, there is only dZ, so its construction will be direct.
 #postwalk
 prewalk(zcf) do a   #
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

function extractZCJacDep(counter::Int,zcf::Symbol,zcjac :: Vector{Vector{Int}},SZ ::Dict{Int,Set{Int}},dZ :: Dict{Int,Set{Int}}) # case if t>0
end




"""
    createSZVect(SZ :: Dict{Int64, Set{Int64}},::Val{T}) where {T} 

constructs the zero-crossing dependency to state variables as a vector from the existing dictionary SZ resulted from the [`extractZCJacDep`](@ref) function. The continuous variables are the keys and the zero-crossing are the values.\n
   # Example:
```jldoctest
using QuantizedSystemSolver
(SZ, T) = (Dict{Int64, Set{Int64}}(2 => Set([2]), 1 => Set([1])), 10)
szVec=QuantizedSystemSolver.createSZVect(SZ, Val(T))
string(szVec)

# output

"[[1], [2], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[]]"
```
"""
function createSZVect(SZ :: Dict{Int64, Set{Int64}},::Val{T}) where {T}
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
    createdDVect(dD::Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}},::Val{D}) where {D}

constructs the State to derivative dependency to discrete variables as a vector from the existing dictionary dD resulted from the [`extractJacDepNormal`](@ref) and the [`extractJacDepLoop`](@ref) functions. The discrete variables are the keys and the differential equations are the values. This dependency is needed only in [`createDependencyToEventsDiscr`](@ref)  .\n
   # Example:
```jldoctest
using QuantizedSystemSolver
(dD, D) = (Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(2 => Set([10, 1]), 1 => Set([:((2, 9))])), 2)
dDVect=QuantizedSystemSolver.createdDVect(dD, Val(D) )
string(dDVect)

# output

"[[2, 3, 4, 5, 6, 7, 8, 9], [10, 1]]"
```
"""
function createdDVect(dD::Dict{Union{Int64,Symbol,Expr}, Set{Union{Int64, Expr, Symbol}}},::Val{D}) where {D}
  dDVect = Vector{Vector{Int}}(undef, D)
  @show dDVect, D
  for ii=1:D
      dDVect[ii]=Vector{Int}()# define it so i can push elements as i find them below
  end
  for dictElement in dD
    if dictElement[1] isa Int
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
     dDVect[dictElement[1]]=temp ###################################
    elseif dictElement[1] isa Symbol
       setElement=first(dictElement[2])
          for j_=(setElement.args[1]):(setElement.args[2])  
           push!(dDVect[j_],j_)          ########################
       end
      elseif dictElement[1] isa Expr
        setElement=first(dictElement[2])
          for j_=(setElement.args[1]):(setElement.args[2])  
            fa= postwalk(a -> a isa Symbol && a==:i ? j_ : a, dictElement[1]) # change each symbol i to exact number 
          push!(dDVect[Int(eval(fa))],j_)
       end
    end
  end
  dDVect
end

"""
  handleEvents(argI::Expr,eventequs::Vector{Expr},length_zcequs::Int64,evsArr::Vector{EventDependencyStruct})

Handles events in the quantized system solver.

# Arguments
- `argI::Expr`: An expression representing the 'if-statement'.
- `eventequs::Vector{Expr}`: A vector of expressions representing the event equations.
- `length_zcequs::Int64`: current number of treated zero-crossing equations. Usezd as index to store the events in order.
- `evsArr`: An array containing EventDependencyStruct  objects.


"""
function handleEvents(argI::Expr,eventequs::Vector{Expr},length_zcequs::Int64,evsArr::Vector{EventDependencyStruct})
                      # argI.args[2] is the condition, which we do not handle here (already handled as zcf in qssProblem)
            # each 'if-statmets' has 2 events (arg[2]=posEv and arg[3]=NegEv) each pos or neg event has a function...
              if length(argI.args)==2  #if user only wrote the positive evnt, here I added the negative event wich does nothing
                #nothingexpr = quote nothing end # neg dummy event 
                nothingexpr = Expr(:block,nothing) # does not require remove_linenums
                push!(argI.args, nothingexpr)
               # Base.remove_linenums!(argI.args[3])

                nothingexpr = Expr(:block,nothing)
            end 
            #pos event
            newPosEventExprToFunc=changeExprToFirstValue(argI.args[2])  #change u[1] to u[1][0]  # pos ev can't be a symbol ...later maybe add check anyway
           newPosEventExprToFunc=changeT(newPosEventExprToFunc) # t -> t[0]
            #@show newPosEventExprToFunc
            push!(eventequs,newPosEventExprToFunc) 
            #neg eve
            if argI.args[3].args[1] isa Expr # argI.args[3] isa expr and is either an equation or :block symbol end ...so lets check .args[1]
                newNegEventExprToFunc=changeExprToFirstValue(argI.args[3])
                newNegEventExprToFunc=changeT(newNegEventExprToFunc) # t -> t[0]
                push!(eventequs,newNegEventExprToFunc) 
            else
                push!(eventequs,argI.args[3]) #symbol nothing
            end
            #after constructing the equations we move to dependencies: we need to change A[n] to An so that they become symbols
            posEvExp =  argI.args[2]
            negEvExp =  argI.args[3]
           # @show posEvExp
            indexPosEv = 2 * length_zcequs  - 1 # store events in order
            indexNegEv =  indexPosEv + 1 
                                  #------------------pos Event--------------------#
            posEv_disArrLHS= Set{Int}()  
            posEv_conArrLHS= Set{Int}() 
            posEv_conArrRHS=Set{Int}() 
            #posEv_conArrRHS=Set{Int}()    #to be used inside intgrator to updateOtherQs (intgrateState) before executing the event there is no discArrRHS because p is not changing overtime to be updated      
            
            assigns = find_assignments(posEvExp)
            for j in assigns
            #for j = 1:length(posEvExp.args)  # j coressponds the number of statements under one posEvent
             # @show posEvExp.args[j] 
                #if (posEvExp.args[j]  isa Expr &&  posEvExp.args[j].head in [:(=),:+=, :-=, :*=, :/=]) 
                        poslhs=j.args[1];posrhs=j.args[2]
                       # @show poslhs.args[1]
                    if (poslhs  isa Expr &&  poslhs.head == :ref && (poslhs.args[1]==:q || poslhs.args[1]==:p))  
                     # @show poslhs    
                       if poslhs.args[1]==:q
                            push!(posEv_conArrLHS,poslhs.args[2])
                            if j.head in [:+=, :-=, :*=, :/=] # lhs+=rhs is lhs=lhs+rhs so lhs need to be accounted for
                              push!(posEv_conArrRHS,  (poslhs.args[2])) 
                            end
                        else # lhs is a disc var 
                            push!(posEv_disArrLHS,poslhs.args[2])
                        end
                        postwalk(posrhs) do a   #
                            if a isa Expr && a.head == :ref && a.args[1]==:q# 
                                push!(posEv_conArrRHS,  (a.args[2]))  #       ....             
                            end
                            return a 
                        end
                    end
                #end
            end
                                        #------------------neg Event--------------------#
            negEv_disArrLHS= Set{Int}()#
            negEv_conArrLHS= Set{Int}()# 
            negEv_conArrRHS=Set{Int}()#to be used inside intgrator to updateOtherQs (intgrateState) before executing the event there is no discArrRHS because p is not changing overtime to be updated      
            if negEvExp.args[1] != :nothing
                assigns = find_assignments(negEvExp)
                for j in assigns
                #for j = 1:length(negEvExp.args)  # j coressponds the number of statements under one negEvent
                  #if (negEvExp.args[j]  isa Expr &&  negEvExp.args[j].head in [:(=),:+=, :-=, :*=, :/=])  
                    neglhs=j.args[1];negrhs=j.args[2]                   
                    if (neglhs  isa Expr &&  neglhs.head == :ref && (neglhs.args[1]==:q || neglhs.args[1]==:p))  
                      
                        if neglhs.args[1]==:q
                            push!(negEv_conArrLHS,neglhs.args[2])
                            if j.head in [:+=, :-=, :*=, :/=] # lhs+=rhs is lhs=lhs+rhs so lhs need to be accounted for
                              push!(negEv_conArrRHS,  (neglhs.args[2])) 
                            end
                        else # lhs is a disc var 
                            push!(negEv_disArrLHS,neglhs.args[2])
                        end
                        postwalk(negrhs) do a   #
                            if a isa Expr && a.head == :ref && a.args[1]==:q# 
                                push!(negEv_conArrRHS,  (a.args[2]))  #                    
                            end
                            return a 
                        end#end postwalk
                    end#end if
                  #end
                end#end for
            end 
            structposEvent = EventDependencyStruct(indexPosEv, collect(posEv_conArrLHS), collect(posEv_disArrLHS),collect(posEv_conArrRHS)) # posEv_conArr is vect 
            push!(evsArr, structposEvent)
            structnegEvent = EventDependencyStruct(indexNegEv, collect(negEv_conArrLHS), collect(negEv_disArrLHS),collect(negEv_conArrRHS))
            push!(evsArr, structnegEvent)

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
