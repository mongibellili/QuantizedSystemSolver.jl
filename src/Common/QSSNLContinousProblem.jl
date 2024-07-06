"""NLODEContProblem{PRTYPE,T,Z,Y,CS}
A struct that holds the Problem of a system of ODEs. It has the following fields:\n
    - prname: The name of the problem\n
    - prtype: The type of the problem\n
    - a: The size of the problem\n
    - b: The number of zero crossing functions\n
    - c: The number of discrete events\n
    - cacheSize: The size of the cache\n
    - initConditions: The initial conditions of the problem\n
    - eqs: The function that holds all the ODEs\n
    - jac: The Jacobian dependency\n
    - SD: The state derivative dependency\n
    - exactJac: The exact Jacobian function\n
    - jacDim: The Jacobian dimension (if sparsity to be exploited)
"""
struct NLODEContProblem{PRTYPE,T,Z,Y,CS}<: NLODEProblem{PRTYPE,T,Z,Y,CS} 
    prname::Symbol # problem name used to distinguish printed results
    prtype::Val{PRTYPE} # problem type: not used but created in case in the future we want to handle problems differently
    a::Val{T} #problem size based on number of vars: T is used not a: 'a' is a mute var
    b::Val{Z} #number of Zero crossing functions (ZCF) based on number of 'if statements': Z is used not b: 'b' is a mute var
    c::Val{Y} #number of discrete events=2*ZCF: Z is used not c: 'c' is a mute var
    cacheSize::Val{CS}# CS= cache size is used  : 'cacheSize' is a mute var
    initConditions::Vector{Float64}  # 
    eqs::Function#function that holds all ODEs
    jac::Vector{Vector{Int}}#Jacobian dependency..I have a der and I want to know which vars affect it...opposite of SD...is a vect for direct method (later @resumable..closure..for saved method)
    SD::Vector{Vector{Int}}#  I have a var and I want the der that are affected by it
    exactJac::Function  # used only in the implicit intgration
    #map::Function  # if sparsity to be exploited: it maps a[i][j] to a[i][γ] where γ usually=1,2,3...(a small int)
    jacDim::Function # if sparsity to be exploited: gives length of each row
end


# to create NLODEContProblem above
function NLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{0},::Val{0}, initConditions::Vector{Float64} ,du::Symbol,symDict::Dict{Symbol,Expr})where {T}
    if VERBOSE println("nlodeprobfun  T= $T") end
    equs=Dict{Union{Int,Expr},Expr}() #du[int] or du[i]==du[a<i<b]==du[(a:b)]
    jac = Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}()# datastrucutre 'Set' used because we do not want to insert an existing varNum
    exacteJacExpr = Dict{Expr,Union{Float64,Int,Symbol,Expr}}()
    # NO need to constrcut SD...SDVect will be extracted from Jac
    num_cache_equs=1#initial cachesize::will hold the number of caches (vectors) of the longest equation
  
    for argI in odeExprs.args
        #only diff eqs: du[]= number || ref || call 
        if argI isa Expr &&  argI.head == :(=)  && argI.args[1] isa Expr && argI.args[1].head == :ref && argI.args[1].args[1]==du#expr LHS=RHS and LHS is du
            y=argI.args[1];rhs=argI.args[2]
            varNum=y.args[2] # order/index of variable
            if rhs isa Number # rhs of equ =number  
                equs[varNum]=:($((transformFSimplecase(:($(rhs)))))) #change rhs from N to cache=taylor0=[[N,0,0],2] for order 2 for exple
            elseif rhs isa Symbol # case du=t   
                #equs[varNum ]=quote $rhs end
                equs[varNum ]=:($((transformFSimplecase(:($(rhs))))))#
            elseif rhs.head==:ref #rhs is only one var
                extractJacDepNormal(varNum,rhs,jac,exacteJacExpr ,symDict ) #extract jacobian and exactJac approx form from normal equation
                equs[varNum ]=:($((transformFSimplecase(:($(rhs))))))#change rhs from q[i] to cache=q[i] ...just put taylor var in cache
            else #rhs head==call...to be tested later for  math functions and other possible scenarios or user erros                 
                extractJacDepNormal(varNum,rhs,jac,exacteJacExpr,symDict ) 
                temp=(transformF(:($(rhs),1))).args[2]  #temp=number of caches distibuted...rhs already changed inside
                if num_cache_equs<temp 
                        num_cache_equs=temp
                end 
                equs[varNum]=rhs
            end 
        elseif @capture(argI, for counter_ in b_:niter_ loopbody__ end) #case where diff equation is an expr in for loop
             specRHS=loopbody[1].args[2]
             extractJacDepLoop(b,niter,specRHS,jac,exacteJacExpr,symDict  )  #extract jacobian and SD dependencies from loop 
             temp=(transformF(:($(specRHS),1))).args[2]
                if num_cache_equs<temp 
                    num_cache_equs=temp
                end 
               equs[:(($b,$niter))]=specRHS            
        else#end of equations and user enter something weird...handle later
               # error("expression $x: top level contains only expressions 'A=B' or 'if a b' or for loop... ")#
        end #end for x in   
    end #end for args #########################################################################################################

    fname= :f #default problem name
    #path="./temp.jl" #default path 
    if odeExprs.args[1] isa Expr && odeExprs.args[1].args[2] isa Expr && odeExprs.args[1].args[2].head == :tuple#user has to enter problem info in a tuple
        fname= odeExprs.args[1].args[2].args[1]
        #path=odeExprs.args[1].args[2].args[2]
    end
 
    exacteJacfunction=createExactJacFun(exacteJacExpr,fname)
    
   #=  open("./temp.jl", "a") do io    
        println(io,string(exacteJacfunction)) 
    end =#

    exactJacfunctionF=@RuntimeGeneratedFunction(exacteJacfunction)
    
   
    diffEqfunction=createContEqFun(equs,fname)# diff equations before this are stored in a dict:: now we have a giant function that holds all diff equations
    jacVect=createJacVect(jac,Val(T)) #jacobian dependency
    SDVect=createSDVect(jac,Val(T))  # state derivative dependency
   # mapFun=createMapFun(jac,fname)  # for sparsity
    jacDimFunction=createJacDimensionFun(jac,fname) #for sparsity
    diffEqfunctionF=@RuntimeGeneratedFunction(diffEqfunction) # @RuntimeGeneratedFunction changes a fun expression to actual fun without running into world age problems
   # mapFunF=@RuntimeGeneratedFunction(mapFun)
    jacDimFunctionF=@RuntimeGeneratedFunction(jacDimFunction)
    prob=NLODEContProblem(fname,Val(1),Val(T),Val(0),Val(0),Val(num_cache_equs),initConditions,diffEqfunctionF,jacVect,SDVect,exactJacfunctionF,jacDimFunctionF)# prtype type 1...prob not saved and struct contains vects

end



function createContEqFun(equs::Dict{Union{Int,Expr},Expr},funName::Symbol)
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
    def=Dict{Symbol,Any}()
    def[:head] = :function
    def[:name] = funName   
    def[:args] = [:(i::Int),:(q::Vector{Taylor0}),:(t::Taylor0),:(cache::Vector{Taylor0})]
    def[:body] = myex1  
    functioncode=combinedef(def)
   # @show functioncode;functioncode
end

function createJacDimensionFun(jac:: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},funName::Symbol)
    ss="if i==0 return 0\n"
    for dictElement in jac
    #=  Base.remove_linenums!(dictElement[1])
        Base.remove_linenums!(dictElement[2]) =#
       # counterJac=1
        if dictElement[1] isa Int
            ss*="elseif i==$(dictElement[1])  \n"
            ss*="return $(length(dictElement[2]))  \n"
           # ss*=" return nothing \n"
        elseif dictElement[1] isa Expr
            ss*="elseif $(dictElement[1].args[1])<=i<=$(dictElement[1].args[2])  \n"
            
                ss*="return $(length(dictElement[2]))  \n"
              
        end     
    end
    ss*=" end \n"         
    myex1=Meta.parse(ss)
    Base.remove_linenums!(myex1)
    def1=Dict{Symbol,Any}() #any changeto Union{expr,Symbol}  ????
    def1[:head] = :function
    def1[:name] = Symbol(:jacDimension,funName)  
    def1[:args] = [:(i::Int)]
    def1[:body] = myex1
    functioncode1=combinedef(def1)
end

function createMapFun(jac:: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},funName::Symbol)
    ss="if i==0 return nothing\n"
    for dictElement in jac
    #=  Base.remove_linenums!(dictElement[1])
        Base.remove_linenums!(dictElement[2]) =#
       # counterJac=1
        if dictElement[1] isa Int
            ss*="elseif i==$(dictElement[1]) \n"
            ns=sort!(collect(dictElement[2]))
            ss*="if j==0 return nothing\n"
            iter=1
            for k in ns
                ss*="elseif j==$k \n"
                ss*="cache[1]=$iter \n"
                ss*=" return nothing \n"
                iter+=1
            end


            ss*=" end \n"
        elseif dictElement[1] isa Expr
            ss*="elseif $(dictElement[1].args[1])<=i<=$(dictElement[1].args[2])  \n"
            @show dictElement[2]
            ns=sort!(collect(dictElement[2]))
            ss*="if j==0 return nothing\n"
            iter=1
            for k in ns
                ss*="elseif j==$k \n"
                ss*="cache[1]=$iter \n"
                ss*=" return nothing \n"
                iter+=1
            end


            ss*=" end \n"
        end     
    end
    ss*=" end \n"         
    myex1=Meta.parse(ss)
    Base.remove_linenums!(myex1)
    def1=Dict{Symbol,Any}() #any changeto Union{expr,Symbol}  ????
    def1[:head] = :function
    def1[:name] = Symbol(:map,funName)  
    def1[:args] = [:(cache::MVector{1,Int}),:(i::Int),:(j::Int)]
    def1[:body] = myex1
    functioncode1=combinedef(def1)
end

function createExactJacFun(jac:: Dict{Expr,Union{Float64,Int,Symbol,Expr}},funName::Symbol)
    ss="if i==0 return nothing\n"
    for dictElement in jac
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
    def1=Dict{Symbol,Any}() #any changeto Union{expr,Symbol}  ????
    def1[:head] = :function
    def1[:name] = Symbol(:exactJac,funName)  
    def1[:args] = [:(q::Vector{Taylor0}),:(d::Vector{Float64}),:(cache::MVector{1,Float64}),:(i::Int),:(j::Int),:(t::Float64)]
    def1[:body] = myex1
    functioncode1=combinedef(def1)
end

function createExactJacDiscreteFun(jac:: Dict{Expr,Union{Float64,Int,Symbol,Expr}},funName::Symbol)
    ss="if i==0 return nothing\n"
    for dictElement in jac
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
    def1=Dict{Symbol,Any}() #any changeto Union{expr,Symbol}  ????
    def1[:head] = :function
    def1[:name] = Symbol(:exactJac,funName)  
    def1[:args] = [:(q::Vector{Taylor0}),:(d::Vector{Float64}),:(cache::MVector{1,Float64}),:(i::Int),:(j::Int),:(t::Float64)] # in jac t does not need to be a taylor
    def1[:body] = myex1
    functioncode1=combinedef(def1)
end

function createJacVect(jac:: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},::Val{T}) where{T}# 
    jacVect = Vector{Vector{Int}}(undef, T)
    for i=1:T
        jacVect[i]=Vector{Int}()# define it so i can push elements as i find them below
    end
    for dictElement in jac
        if dictElement[1] isa Int  #jac[varNum]=jacSet
           jacVect[dictElement[1]]=collect(dictElement[2])
        elseif dictElement[1] isa Expr  #jac[:(($b,$niter))]=jacSet
            for j_=(dictElement[1].args[1]):(dictElement[1].args[2])  
                temp=Vector{Int}()
                for element in dictElement[2]
                    if element isa Expr || element isa Symbol#can split symbol alone since no need to postwalk
                        fa= postwalk(a -> a isa Symbol && a==:i ? j_ : a, element) # change each symbol i to exact number 
                        push!(temp,eval(fa))
                    else #element is int
                        push!(temp,element)
                    end
                end
                jacVect[j_]=temp
            end
        end     
    end
    jacVect
end
function createSDVect(jac:: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},::Val{T}) where{T}
    sdVect = Vector{Vector{Int}}(undef, T)
    for ii=1:T
        sdVect[ii]=Vector{Int}()# define it so i can push elements as i find them below
    end
    #SD = Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(2 => Set([1]), 10 => Set([10]), :((2, 9)) => Set([:k, :(k - 1), :(k + 1)]), 9 => Set([10]), 1 => Set([1]))
    for dictElement in jac
        if dictElement[1] isa Int # key is an int
            for k in dictElement[2]   #elments values 
                push!(sdVect[k],dictElement[1])
            end
        elseif dictElement[1] isa Expr # key is an expression
            for j_=(dictElement[1].args[1]):(dictElement[1].args[2])  # j_=b:N this can be expensive when N is large::this is why it is recommended to use a function createsd (save) for large prob
                for element in dictElement[2]
                    if element isa Expr || element isa Symbol#element=
                    fa= postwalk(a -> a isa Symbol && a==:i ? j_ : a, element)
                    push!(sdVect[eval(fa)],j_)
                    else#element is int
                        push!(sdVect[element],j_)
                    end
                end
               
            end
        end     
    end
    sdVect
end


#helper funs used in lines 136 and 150
function Base.isless(ex1::Expr, ex2::Expr)#:i is the mute var that prob is now using
    fa= postwalk(a -> a isa Symbol && a==:i ? 1 : a, ex1)# check isa symbol not needed
    fb=postwalk(a -> a isa Symbol && a==:i ? 1 : a, ex2)
    eval(fa)<eval(fb)
  end
  function Base.isless(ex1::Expr, ex2::Symbol)
    fa= postwalk(a -> a isa Symbol && a==:i ? 1 : a, ex1)
    eval(fa)<1
  end
  function Base.isless(ex1::Symbol, ex2::Expr)
    fa= postwalk(a -> a isa Symbol && a==:i ? 1 : a, ex2)
    1<eval(fa)
  end