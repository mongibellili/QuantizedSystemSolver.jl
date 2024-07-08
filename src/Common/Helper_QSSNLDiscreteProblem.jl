struct EventDependencyStruct
  id::Int
  evCont::Vector{Int} #index tracking used for HD & HZ. Also it is used to update q,quantum,recomputeNext when x is modified in an event
  evDisc::Vector{Int} #index tracking used for HD & HZ.
  evContRHS::Vector{Int} #index tracking used to update other Qs before executing the event
end
# the following functions handle discrete problems
 # to extract jac & dD....SD can be extracted from jac later
function extractJacDepNormal(varNum::Int,rhs::Union{Symbol,Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},exacteJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},symDict::Dict{Symbol,Expr},dD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 
  jacSet=Set{Union{Int,Symbol,Expr}}()
  #jacDiscrSet=Set{Union{Int,Symbol,Expr}}()
  m=postwalk(rhs) do a   #
      if a isa Expr && a.head == :ref && a.args[1]==:q#  q[2]
          push!(jacSet,  (a.args[2]))  # du[varNum=1]=rhs=u[5]+u[2] : 2 and 5 are stored in jacset one at a time
          a=eliminateRef(a)#q[i] -> qi
      elseif a isa Expr && a.head == :ref && a.args[1]==:d# d[3]
         # push!(jacDiscrSet,  (a.args[2]))  #
          dDset=Set{Union{Int,Symbol,Expr}}()
          if haskey(dD, (a.args[2]))    # dict dD already contains key a.args[2] (in this case var 3)
              dDset=get(dD,(a.args[2]),dDset) # if var d3 first time to influence some var, dDset is empty, otherwise get its set of influences 
          end
          push!(dDset,  varNum)      # d3 also influences varNum
          dD[(a.args[2])]=dDset       #update the dict
          a=eliminateRef(a)#d[i] -> di
      end
      return a 
  end
  # extract the jac
  basi = convert(Basic, m) # m ready: all refs are symbols
  for i in jacSet  # jacset contains vars in RHS
    symarg=symbolFromRef(i) # specific to elements in jacSet: get q1 from 1 for exple
    coef = diff(basi, symarg) # symbolic differentiation: returns type Basic
    coefstr=string(coef);coefExpr=Meta.parse(coefstr)#convert from basic to expression
    jacEntry=restoreRef(coefExpr,symDict)# get back ref: qi->q[i][0]  ...0 because later in exactJac fun cache[1]::Float64=jacEntry
    exacteJacExpr[:(($varNum,$i))]=jacEntry # entry (varNum,i) is jacEntry
  end
  if length(jacSet)>0 jac[varNum]=jacSet end
  #if length(jacDiscrSet)>0 jacDiscr[varNum]=jacDiscrSet end
end

# like above except (b,niter) instead of varNum
function extractJacDepLoop(b::Int,niter::Int,rhs::Union{Symbol,Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},exacteJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},symDict::Dict{Symbol,Expr},dD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 
  jacSet=Set{Union{Int,Symbol,Expr}}()
  m=postwalk(rhs) do a   
      if a isa Expr && a.head == :ref && a.args[1]==:q# 
              push!(jacSet,  (a.args[2]))  #
              a=eliminateRef(a)#q[i] -> qi
      elseif a isa Expr && a.head == :ref && a.args[1]==:d
        if a.args[2] isa Int  #for now allow only d[integer]
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
    exacteJacExpr[:((($b,$niter),$i))]=jacEntry
  end
  jac[:(($b,$niter))]=jacSet
end

function extractZCJacDepNormal(counter::Int,zcf::Expr,zcjac :: Vector{Vector{Int}},SZ ::Dict{Int,Set{Int}},dZ :: Dict{Int,Set{Int}}) 
  zcjacSet=Set{Int}()
 #zcjacDiscrSet=Set{Int}()
  postwalk(zcf) do a   #
      if a isa Expr && a.head == :ref && a.args[1]==:q# 
          push!(zcjacSet,  (a.args[2]))  #
          SZset=Set{Int}()
          if haskey(SZ, (a.args[2]))
              SZset=get(SZ,(a.args[2]),SZset)
          end
          push!(SZset,  counter)
          SZ[(a.args[2])]=SZset
      elseif a isa Expr && a.head == :ref && a.args[1]==:d# 
         # push!(zcjacDiscrSet,  (a.args[2]))  #
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
  #push!(zcjacDiscr,collect(zcjacDiscrSet))
end

function createDependencyToEventsDiscr(dD::Vector{Vector{Int}},dZ::Dict{Int64, Set{Int64}},eventDep::Vector{EventDependencyStruct}) 
    Y=length(eventDep)
    lendD=length(dD)
    HD2 = Vector{Vector{Int}}(undef, Y)
    HZ2 = Vector{Vector{Int}}(undef, Y)
      for ii=1:Y
        HD2[ii] =Vector{Int}()# define it so i can push elements as i find them below
        HZ2[ii] =Vector{Int}()# define it so i can push elements as i find them below
      end
    for j=1:Y
      hdSet=Set{Int}()
      hzSet=Set{Int}()
        evdiscrete=eventDep[j].evDisc
        for i in evdiscrete
             if i<=lendD  
              for k in dD[i]
                push!(hdSet,k)
              end
            end
              tempSet=Set{Int}()
              if haskey(dZ, i)
                tempSet=get(dZ,i,tempSet)
              end
              for kk in tempSet
                push!(hzSet,kk)
              end
        end
        HD2[j] =collect(hdSet)# define it so i can push elements as i find them below
        HZ2[j] =collect(hzSet)
    end #end for j  (events)
    return (HZ2,HD2)
end 

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
      evContin=eventDep[j].evCont
      for i in evContin
            for k in SD[i]
              push!(hdSet,k)
            end
            tempSet=Set{Int}()
            if haskey(sZ, i)
              tempSet=get(sZ,i,tempSet)
            end
            for kk in tempSet
              push!(hzSet,kk)
            end
      end
      HD2[j] =collect(hdSet)# define it so i can push elements as i find them below
      HZ2[j] =collect(hzSet)
  end #end for j  (events)
  return (HZ2,HD2)
end 

#function unionDependency(HZ1::SVector{Y,SVector{Z,Int}},HZ2::SVector{Y,SVector{Z,Int}})where{Z,Y}
function unionDependency(HZ1::Vector{Vector{Int}},HZ2::Vector{Vector{Int}})
  Y=length(HZ1)
  HZ = Vector{Vector{Int}}(undef, Y)
    for ii=1:Y
      HZ[ii] =Vector{Int}()# define it so i can push elements as i find them below
    end
  for j=1:Y
    hzSet=Set{Int}()
    for kk in HZ1[j]
      push!(hzSet,kk)
    end
    for kk in HZ2[j]
      push!(hzSet,kk)
    end
    HZ[j]=collect(hzSet)
  end
  HZ 
end
