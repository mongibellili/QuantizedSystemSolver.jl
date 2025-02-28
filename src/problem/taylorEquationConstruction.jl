

"""
    transformFSimplecase(ex::Union{Float64,Int64,Expr,Symbol})

 transforms expressions of the right hand side of differential equations and zero-crossing functions to personalized ones that use caching to a form that can be used in the [`createT`](@ref) function. The right hand side of the equations should be a number or a q[i] term.
# Example:   
```jldoctest
using QuantizedSystemSolver
ex=:(q[2])
newEx=QuantizedSystemSolver.transformFSimplecase(ex);


# output

:(createT(q[2], cache[1])) 

```
"""  
function transformFSimplecase(ex::Union{Float64,Int64,Expr,Symbol})
    #  it s easier and healthier to leave the big case alone in one prewalk (below) when returning the size of the distributed cache
    ex=Expr(:call, :createT,ex)# case where rhs of eq is a number needed to change to taylor (in cache) to be used for qss ||  #case where rhs is q[i]...return cache that contains it
    cachexpr = Expr(:ref, :cache)   # add cache[1] to createT(ex,....)
    push!(cachexpr.args,1)
    push!( ex.args, cachexpr)    
    return ex
end 

"""
    transformF(ex::Expr)

 transforms expressions of the right hand side of differential equations and zero-crossing functions to personalized ones that use caching to a form that can be used in functions like  [`addT`](@ref), [`subT`](@ref), [`mulT`](@ref), [`muladdT`](@ref). The right hand side of the equations can be any form of expression.
# Example:   
```jldoctest
using QuantizedSystemSolver
ex=:(q[2] - 2.0*q[1]*q[2],1)
newEx=QuantizedSystemSolver.transformF(ex);


# output

:((subT(q[2], mulTT(2.0, q[1], q[2], cache[2], cache[3]), cache[1]), 3))

```
"""  
function transformF(ex::Expr)# name to be changed later....i call this funciton in the docs the function that does the transformation
   cachexpr_lengthtracker = Expr(:mongi)# an expr to track the number of distributed caches
    prewalk(ex) do x
      #############################minus sign#############################################
      if x isa Expr && x.head == :call && x.args[1] == :- && length(x.args) == 3 && !(x.args[2] isa Number && x.args[3] isa Number) && !(x.args[2] isa Int || x.args[3] isa Int)  #last 2 to avoid changing u[i-1]  #MacroTools.isexpr(x, :call) replaces the begining
              if x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :- && length(x.args[2].args) == 3
                  push!(x.args, x.args[3])#  args[4]=c
                  x.args[3] = x.args[2].args[3]# args[3]=b
                  x.args[2] = x.args[2].args[2]# args[2]=a
                  x.args[1] = :subsub
                  push!(cachexpr_lengthtracker.args,:b) #b is anything ...not needed, just to fill vector tracker
                  cachexpr = Expr(:ref, :cache) #prepare cache
                  push!(cachexpr.args,length(cachexpr_lengthtracker.args))#constrcut cache with index cache[1]
                  push!(x.args, cachexpr)
              elseif x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :+ && length(x.args[2].args) == 3 
                  push!(x.args, x.args[3])#  args[4]=c
                  x.args[3] = x.args[2].args[3]# args[3]=b
                  x.args[2] = x.args[2].args[2]# args[2]=a
                  x.args[1] = :addsub        
                  push!(cachexpr_lengthtracker.args,:b)
                  cachexpr = Expr(:ref, :cache)
                  push!(cachexpr.args,length(cachexpr_lengthtracker.args))
                  #cachexpr.args[2] = index[i]
                  push!(x.args, cachexpr)
              elseif x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :* && length(x.args[2].args) == 3
                  push!(x.args, x.args[3])#  args[4]=c
                  x.args[3] = x.args[2].args[3]# args[3]=b
                  x.args[2] = x.args[2].args[2]# args[2]=a
                  x.args[1] = :mulsub 
                  push!(cachexpr_lengthtracker.args,:b)
                  cachexpr1 = Expr(:ref, :cache)
                  push!(cachexpr1.args,length(cachexpr_lengthtracker.args))
                  push!(x.args, cachexpr1)
              else
                  x.args[1] = :subT  # symbol changed cuz avoid type taylor piracy    
                  push!(cachexpr_lengthtracker.args,:b)
                  cachexpr = Expr(:ref, :cache)
                  push!(cachexpr.args,length(cachexpr_lengthtracker.args))
                  push!(x.args, cachexpr)
              end
      elseif x isa Expr && x.head == :call && x.args[1] == :- && length(x.args) == 2  && !(x.args[2] isa Number)#negate
              x.args[1] = :negateT  # symbol changed cuz avoid type taylor piracy
              push!(cachexpr_lengthtracker.args,:b)
              cachexpr = Expr(:ref, :cache)
              push!(cachexpr.args,length(cachexpr_lengthtracker.args))
              push!(x.args, cachexpr)
      ############################### plus sign#######################################
      elseif x isa Expr && x.head == :call && x.args[1] == :+ && length(x.args) == 3 && !(x.args[2] isa Number && x.args[3] isa Number) && !(x.args[2] isa Int || x.args[3] isa Int)#&& typeof(x.args[2])!=Int && typeof(x.args[3])!=Int #last 2 to avoid changing u[i+1]
            if x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :- && length(x.args[2].args) == 3        
                  push!(x.args, x.args[3])#  args[4]=c
                  x.args[3] = x.args[2].args[3]# args[3]=b
                  x.args[2] = x.args[2].args[2]# args[2]=a
                  x.args[1] = :subadd#:µ  # £  § ....
                  push!(cachexpr_lengthtracker.args,:b)
                  cachexpr = Expr(:ref, :cache)
                  push!(cachexpr.args,length(cachexpr_lengthtracker.args))
                  push!(x.args, cachexpr)
            elseif x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :* && length(x.args[2].args) == 3
                  push!(x.args, x.args[3])#  args[4]=c
                  x.args[3] = x.args[2].args[3]# args[3]=b
                  x.args[2] = x.args[2].args[2]# args[2]=a
                  x.args[1] = :muladdT#:µ  # £  § ....
                  push!(cachexpr_lengthtracker.args,:b) 
                  cachexpr1 = Expr(:ref, :cache)
                  push!(cachexpr1.args,length(cachexpr_lengthtracker.args))
                  push!(x.args, cachexpr1)
            else
                  x.args[1] = :addT
                  push!(cachexpr_lengthtracker.args,:b)
                  cachexpr = Expr(:ref, :cache)
                  push!(cachexpr.args,length(cachexpr_lengthtracker.args))
                  push!(x.args, cachexpr)
            end
      elseif x isa Expr && x.head == :call && x.args[1] == :+ && (4 <= length(x.args) <= 9) # testing it was allocating anyway after 9
                x.args[1] = :addT
                push!(cachexpr_lengthtracker.args,:b)
                cachexpr = Expr(:ref, :cache) 
                push!(cachexpr.args,length(cachexpr_lengthtracker.args))
                push!(x.args, cachexpr)
     ############################### multiply sign#######################################
      elseif x isa Expr && x.head == :call && (x.args[1] == :*) && (3 == length(x.args))  && !(x.args[2] isa Number && x.args[3] isa Number) && !(x.args[2] isa Int || x.args[3] isa Int)#&& typeof(x.args[2])!=Int && typeof(x.args[3])!=Int #last 2 to avoid changing u[i*1]
             x.args[1] = :mulT
             push!(cachexpr_lengthtracker.args,:b)
             cachexpr1 = Expr(:ref, :cache)
             push!(cachexpr1.args,length(cachexpr_lengthtracker.args))
             push!(x.args, cachexpr1)
      elseif x isa Expr && x.head == :call && (x.args[1] == :*) && (4 <= length(x.args) <= 7)# testing it was allocating anyway after 7
              x.args[1] = :mulTT
              push!(cachexpr_lengthtracker.args,:b)
              cachexpr1 = Expr(:ref, :cache)
              push!(cachexpr1.args,length(cachexpr_lengthtracker.args))
              push!(x.args, cachexpr1)
              push!(cachexpr_lengthtracker.args,:b)
              cachexpr2 = Expr(:ref, :cache)   #multiply needs two caches
              push!(cachexpr2.args,length(cachexpr_lengthtracker.args))
              push!(x.args, cachexpr2)
      elseif x isa Expr && x.head == :call && (x.args[1] == :/) && !(x.args[2] isa Number && x.args[3] isa Number) && !(x.args[2] isa Int || x.args[3] isa Int)#&& typeof(x.args[2])!=Int && typeof(x.args[3])!=Int #last 2 to avoid changing u[i/1]
                x.args[1] = :divT   
                push!(cachexpr_lengthtracker.args,:b)
                cachexpr = Expr(:ref, :cache)
                push!(cachexpr.args,length(cachexpr_lengthtracker.args))
                push!(x.args, cachexpr)
      elseif x isa Expr && x.head == :call && (x.args[1] == :exp ||x.args[1] == :log ||x.args[1] == :sqrt ||x.args[1] == :abs )
                push!(cachexpr_lengthtracker.args,:b)
                cachexpr = Expr(:ref, :cache)
                push!(cachexpr.args,length(cachexpr_lengthtracker.args))
                push!(x.args, cachexpr)                
      elseif x isa Expr && x.head == :call && (x.args[1] == :^ )
                x.args[1] = :powerT  # should be like above but for lets keep it now...in test qss.^ caused a problem...fix later
                push!(cachexpr_lengthtracker.args,:b)
                cachexpr = Expr(:ref, :cache)
                push!(cachexpr.args,length(cachexpr_lengthtracker.args))
                push!(x.args, cachexpr)
      elseif x isa Expr && x.head == :call && (x.args[1] == :cos ||x.args[1] == :sin ||x.args[1] == :tan||x.args[1] == :atan )
                push!(cachexpr_lengthtracker.args,:b)
                cachexpr = Expr(:ref, :cache)
                push!(cachexpr.args,length(cachexpr_lengthtracker.args))
                push!(x.args, cachexpr)
                #second cache
                push!(cachexpr_lengthtracker.args,:b)#index2 
                cachexpr2 = Expr(:ref, :cache)   
                push!(cachexpr2.args,length(cachexpr_lengthtracker.args))#construct cache[2]
                push!(x.args, cachexpr2)
     
      elseif x isa Expr && x.head == :call && (x.args[1] == :acos ||x.args[1] == :asin)
                #did not change symbol here
                push!(cachexpr_lengthtracker.args,:b)
                cachexpr = Expr(:ref, :cache)
                push!(cachexpr.args,length(cachexpr_lengthtracker.args))
                push!(x.args, cachexpr)
                #second cache
                push!(cachexpr_lengthtracker.args,:b)#index2
                cachexpr2 = Expr(:ref, :cache)   
                push!(cachexpr2.args,length(cachexpr_lengthtracker.args))#construct cache[2]
                push!(x.args, cachexpr2)
                #third cache
                push!(cachexpr_lengthtracker.args,:b)#index3
                cachexpr2 = Expr(:ref, :cache)   
                push!(cachexpr2.args,length(cachexpr_lengthtracker.args))#construct cache[2]
                push!(x.args, cachexpr2)        
      elseif false #holder for next symbol for other functions..
         
      end
      return x # this is the line that actually enables modification of the original expression (prewalk only)
    end#end prewalk
    ex.args[2]=length(cachexpr_lengthtracker.args)# return the number of caches distributed...to be used in the cache creation 
    return ex
end#end function
  