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
  evCont::Vector{Int} #LHS: index tracking used for HD & HZ. Also it is used to update q,quantum,recomputeNext when x is modified in an event
  evDisc::Vector{Int} #LHS: index tracking used for HD & HZ.
  evContRHS::Vector{Int} #Here we look at the RHS: index tracking used to update other Qs before executing the event
end

function changeEventExprToFirstValue(ex::Expr)
  newEx=postwalk(ex) do a  
      if (a==:t)|| (a isa Expr && a.head == :ref && a.args[1]==:q) 
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
  handleEvents(argI::Expr,eventequs::Vector{Expr},length_zcequs::Int64,evsArr::Vector{EventDependencyStruct})

Handles events in the quantized system solver.

# Arguments
- `argI::Expr`: An expression representing the 'if-statement'.
- `eventequs::Vector{Expr}`: A vector of expressions representing the event equations.
- `length_zcequs::Int64`: current number of treated zero-crossing equations. Usezd as index to store the events in order.
- `evsArr`: An array containing EventDependencyStruct  objects.


"""
function handleEvents(argI::Expr,eventequs::Vector{Expr},length_zcequs::Int64,evsArr::Vector{EventDependencyStruct})
                      
            # each 'if-statmets' has 2 events (arg[2]=posEv and arg[3]=NegEv) each pos or neg event has a function...later i can try one event for zc
              if length(argI.args)==2  #if user only wrote the positive evnt, here I added the negative event wich does nothing
                #nothingexpr = quote nothing end # neg dummy event 
                nothingexpr = Expr(:block,nothing) # does not require remove_linenums
                push!(argI.args, nothingexpr)
               # Base.remove_linenums!(argI.args[3])

                nothingexpr = Expr(:block,nothing)
            end 
            #pos event
            newPosEventExprToFunc=changeEventExprToFirstValue(argI.args[2])  #change u[1] to u[1][0]  # pos ev can't be a symbol ...later maybe add check anyway
           #@show newPosEventExprToFunc
            push!(eventequs,newPosEventExprToFunc) 
            #neg eve
            if argI.args[3].args[1] isa Expr # argI.args[3] isa expr and is either an equation or :block symbol end ...so lets check .args[1]
                newNegEventExprToFunc=changeEventExprToFirstValue(argI.args[3])
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
            posEv_disArrLHS= Vector{Int}()  
            posEv_conArrLHS= Vector{Int}() 
            posEv_conArrRHS=Set{Int}() 
            #posEv_conArrRHS=Set{Int}()    #to be used inside intgrator to updateOtherQs (intgrateState) before executing the event there is no discArrRHS because p is not changing overtime to be updated      
            for j = 1:length(posEvExp.args)  # j coressponds the number of statements under one posEvent
             # @show posEvExp.args[j] 
                if (posEvExp.args[j]  isa Expr &&  posEvExp.args[j].head in [:(=),:+=, :-=, :*=, :/=]) 
                        poslhs=posEvExp.args[j].args[1];posrhs=posEvExp.args[j].args[2]
                       # @show poslhs.args[1]
                    if (poslhs  isa Expr &&  poslhs.head == :ref && (poslhs.args[1]==:q || poslhs.args[1]==:p))  
                     # @show poslhs    
                       if poslhs.args[1]==:q
                            push!(posEv_conArrLHS,poslhs.args[2])
                            if posEvExp.args[j].head in [:+=, :-=, :*=, :/=] # lhs+=rhs is lhs=lhs+rhs so lhs need to be accounted for
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
                end
            end
            #------------------neg Event--------------------#
            negEv_disArrLHS= Vector{Int}()#
            negEv_conArrLHS= Vector{Int}()# 
            negEv_conArrRHS=Set{Int}()#to be used inside intgrator to updateOtherQs (intgrateState) before executing the event there is no discArrRHS because p is not changing overtime to be updated      
            if negEvExp.args[1] != :nothing
                for j = 1:length(negEvExp.args)  # j coressponds the number of statements under one negEvent
                  if (negEvExp.args[j]  isa Expr &&  negEvExp.args[j].head in [:(=),:+=, :-=, :*=, :/=])  
                    neglhs=negEvExp.args[j].args[1];negrhs=negEvExp.args[j].args[2]                   
                    if (neglhs  isa Expr &&  neglhs.head == :ref && (neglhs.args[1]==:q || neglhs.args[1]==:p))  
                      
                        if neglhs.args[1]==:q
                            push!(negEv_conArrLHS,neglhs.args[2])
                            if negEvExp.args[j].head in [:+=, :-=, :*=, :/=] # lhs+=rhs is lhs=lhs+rhs so lhs need to be accounted for
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
                  end
                end#end for
            end 
            structposEvent = EventDependencyStruct(indexPosEv, posEv_conArrLHS, posEv_disArrLHS,collect(posEv_conArrRHS)) # posEv_conArr is vect 
            push!(evsArr, structposEvent)
            structnegEvent = EventDependencyStruct(indexNegEv, negEv_conArrLHS, negEv_disArrLHS,collect(negEv_conArrRHS))
            push!(evsArr, structnegEvent)

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
function createdDVect(dD::Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}},::Val{D}) where {D}
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


function createDiscEqFun(otherCode::Expr,equs::Dict{Union{Int,Expr},Union{Int,Symbol,Expr}},zcequs::Vector{Expr},eventequs::Vector{Expr},fname::Symbol,f)# where{F}
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
  otherCodeCopy=copy(otherCode)# to avoid changing the original
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
  def[:args] = [:(i::Int),:(zc::Int),:(ev::Int),:(q::Vector{Taylor0}),:(p::Vector{Float64}), :(t::Taylor0),:(cache::Vector{Taylor0}),:(f_)]
  def[:body] = otherCodeCopy 
  functioncode=combinedef(def)
 # @show functioncode

end



