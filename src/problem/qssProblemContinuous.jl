
"""
    NLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{D},::Val{0},initConditions::Vector{Float64},du::Symbol,symDict::Dict{Symbol,Expr},tspan::Tuple{Float64, Float64},discrVars::Vector{EM},prbName::Symbol) where {T,D,EM}

This function continues building a continuous problem. it receives an expression and useful info from the main interface. it calls the transform function from the taylorEquationConstruction.jl file to change the AST of all operations to personlized ones and update the needed cache size. It also construct via helper functions the Exact jacobian function, the jacobian dependecy and the state-derivative dependency (opposite of jacobian) as vectors. Finally, it groups all differential equations in one function, and constructs a continous problem from the qssProblemDefinition.jl file.
# Arguments
- `odeExprs::Expr`: The expression of the whole user code in the function defining the problem with names modified and parameters plugged in.
- `Val{T}`: the dimensions of the system of differential equations.    
- `Val{0}`: No zero-crossing functions. pure continous problem.   
- `Val{0}`: No events functions. pure continous problem.  
- `initConditions::Vector{Float64}`: No zero-crossing functions. pure continous problem.  
- `du::Symbol`: to distinguish the start of a differential equations.  
- `symDict::Dict{Symbol,Expr}`: maps a reference expression to a symbol (qi->q[i]).  
- `tspan::Tuple{Float64, Float64}`: stores the initial time and final time of the simulation.  
- `prbName::Symbol`: The problem name as chosen by the user to be carried to the solution for displaying purposes.  
"""
function NLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{D},::Val{0},initConditions::Vector{Float64},du::Symbol,symDict::Dict{Symbol,Expr},tspan::Tuple{Float64, Float64},discrVars::Vector{EM},prbName::Symbol) where {T,D,EM} # 
    if VERBOSE println("nlodeprobfun  T= $T") end
    equs=Dict{Union{Int,Expr},Expr}() #du[int] or du[i]==du[a<i<b]==du[(a:b)]
    jac = Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}()# datastrucutre 'Set' used because we do not want to insert an existing varNum
    exacteJacExpr = Dict{Expr,Union{Float64,Int,Symbol,Expr}}()
    # NO need to constrcut SD...SDVect will be extracted from Jac
    num_cache_equs=1#initial cachesize::will hold the number of caches (vectors) of the longest equation
    for argI in odeExprs.args
        if argI isa Expr &&  argI.head == :(=) && (argI.args[1]== :discrete || argI.args[1]==:p)# old approach use the keyword discrete to define discrete variables
            discrVars = Vector{Float64}(argI.args[2].args)
        #only diff eqs: du[]= number || ref || call 
        elseif argI isa Expr &&  argI.head == :(=)  && argI.args[1] isa Expr && argI.args[1].head == :ref && argI.args[1].args[1]==du#expr LHS=RHS and LHS is du
            y=argI.args[1];rhs=argI.args[2]
            varNum=y.args[2] # order/index of variable
            if rhs isa Number # rhs of equ =number  
                equs[varNum]=transformFSimplecase(rhs) #change rhs from N to cache=taylor0=[[N,0,0],2] for order 2 for exple
            elseif rhs isa Symbol # case du=t   
                equs[varNum ]=transformFSimplecase(rhs)#
            elseif rhs.head==:ref #rhs is only one var
                extractJacDepNormal(varNum,rhs,jac,exacteJacExpr ,symDict ) #extract jacobian and exactJac approx form from normal equation
                equs[varNum ]=transformFSimplecase(rhs)#change rhs from q[i] to cache=q[i] ...just put taylor var in cache
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
        end #end for x in   
    end #end for args #########################################################################################################

    fname= prbName # problem name received as outside function name as written by user
    if odeExprs.args[1] isa Expr && odeExprs.args[1].args[2] isa Expr && odeExprs.args[1].args[2].head == :tuple#user has to enter problem info in a tuple (old maco)
        fname= odeExprs.args[1].args[2].args[1]
    end
    exacteJacfunction=createExactJacFun(exacteJacExpr,fname)
    exactJacfunctionF=@RuntimeGeneratedFunction(exacteJacfunction)
    diffEqfunction=createContEqFun(equs,fname)# diff equations before this are stored in a dict:: now we have a giant function that holds all diff equations
    jacVect=createJacVect(jac,Val(T)) #jacobian dependency
    SDVect=createSDVect(jac,Val(T))  # state derivative dependency
    diffEqfunctionF=@RuntimeGeneratedFunction(diffEqfunction) # @RuntimeGeneratedFunction changes a fun expression to actual fun without running into world age problems
    prob=NLODEContProblemSpan(fname,Val(1),Val(T),Val(D),Val(0),Val(num_cache_equs),initConditions,discrVars,diffEqfunctionF,jacVect,SDVect,exactJacfunctionF,tspan)# prtype type 1...prob not saved and struct contains vects
end
 

#old interface without tspan
function NLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{D},::Val{0},initCond::Vector{Float64},du::Symbol,symDict::Dict{Symbol,Expr})where {T,D}
    discrVars=Vector{Float64}()
    tspan = (0.0,1.0)
    prbName=:_
    probspan=NLodeProblemFunc(odeExprs,Val(T),Val(D),Val(0),initCond,du,symDict,tspan,discrVars,prbName)  
    myodeProblem =NLODEContProblem(probspan.prname,Val(1),Val(T),Val(D),Val(0),probspan.cacheSize,probspan.initConditions,probspan.discreteVars,probspan.eqs,probspan.jac,probspan.SD,probspan.exactJac)
end
