
"""
    NLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{D},::Val{0},initConditions::Vector{Float64},du::Symbol,tspan::Tuple{Float64, Float64},discrVars::Union{Vector{EM}, Tuple{Vararg{EM}}},prbName::Symbol) where {T,D,EM} 

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
function NLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{D},::Val{0},initConditions::Vector{Float64},discrVars::Union{Vector{EM}, Tuple{Vararg{EM}}},preProcessData::PreProcessData) where {T,D,EM} 
   #@show odeExprs
   du=preProcessData.du
   tspan= preProcessData.tspan
   fname= preProcessData.prbName
   mod= preProcessData.mod
   is_top_level= preProcessData.is_top_level
   numHelperFunCalls=preProcessData.numHelperFunCalls
    symDict=Dict{Symbol,Expr}()
    equs=Dict{Union{Int,Expr},Union{Int,Symbol,Expr}}() #du[int] or du[i]==du[a<i<b]==du[(a:b)]
    jac = Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}()# datastrucutre 'Set' used because we do not want to insert an existing varNum
    exacteJacExpr = Dict{Expr,Union{Float64,Int,Symbol,Expr}}()
    # NO need to constrcut SD...SDVect will be extracted from Jac
    num_cache_equs=1#initial cachesize::will hold the number of caches (vectors) of the longest equation
#=     functionEx=quote end #function to be inserted inside the main function
    otherCode=quote end #other code to be inserted inside the main function =#
    functionEx=Expr(:block) #function to be inserted inside the main function
    otherCode=Expr(:block) #other code to be inserted inside the main function
    for argI in odeExprs.args
        if argI isa Expr &&  argI.head == :(=) && (argI.args[1]== :discrete || argI.args[1]==:p)# old approach use the keyword discrete to define discrete variables
            discrVars = Vector{Float64}(argI.args[2].args)
        #only diff eqs: du[]= number || ref || call 
        elseif argI isa Expr &&  argI.head == :(=)  && argI.args[1] isa Expr && argI.args[1].head == :ref && argI.args[1].args[1]==du#expr LHS=RHS and LHS is du
            rhs=argI.args[2];varNum=argI.args[1].args[2] # varnum is the order/index of variable
            
           #=  if !(varNum isa Int) #case du[i]=number
                @show "! ",varNum
                varNum=eval(varNum) #convert to int
                
            else
                @show varNum
            end =#
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
        elseif argI isa Expr && argI.head==:function
            push!(functionEx.args,argI)
           # push!(otherCode.args,argI)
          #=   @show argI
            f = (closure_gen(@__MODULE__))(argI)
            push!(otherCode.args,f) =#
            
        elseif @capture(argI, for counter_ in b_:niter_ loopbody__ end) #case where diff equation is an expr in for loop
           #the 3 commented lines below are for allowing the user to enter many du in a loop
            #= for lpbody in loopbody[1].args
             specRHS=lpbody.args[2]
            end =#
             specRHS=loopbody[1].args[2]
             extractJacDepLoop(Int(eval(b)),Int(eval(niter)),specRHS,jac,exacteJacExpr,symDict  )   #extract jacobian and SD dependencies from loop 
             temp=(transformF(:($(specRHS),1))).args[2]
                if num_cache_equs<temp 
                    num_cache_equs=temp
                end 
               equs[:(($b,$niter))]=specRHS    
        elseif argI.args[1] isa Symbol && argI.args[2] isa Number  # already changed in main and plugged in diff eqs
        elseif  argI.args[1] isa Symbol && argI.args[2] isa Expr && (argI.args[2].head==:call || argI.args[2].head==:ref) #already changed in main and plugged in diff eqs
        elseif argI.args[2] isa Symbol && argI.args[1] isa Expr && argI.args[1].head==:tuple #already changed in main and plugged in diff eqs
        else# keep any other code written by user
            push!(otherCode.args,argI)
        end #end for x in   
    end #end for args #########################################################################################################


   
    jacVect=createJacVect(jac,Val(T)) #jacobian dependency
    SDVect=createSDVect(jac,Val(T))  # state derivative dependency
    #diffEqfunctionF() = (@eval _f()=$diffEqfunction; Base.invokelatest(_f))  
  #=   function diffEqfunctionF(i, q,t,d,cache)
        f = eval(:((i, q,t,d,cache) -> $diffEqfunction))
        return Base.invokelatest(f,i, q,t,d,cache)
    end  =#
    #diffEqfunctionF=RuntimeGeneratedFunction(@__MODULE__,diffEqfunction,opaque_closures=false)
    if length(functionEx.args)==0
        closurefunc=0
        diffEqfunction=createContEqFun(otherCode,equs,fname,closurefunc)# diff equations before this are stored in a dict:: now we have a giant function that holds all diff equations
        exactJacfunction=createExactJacFun(otherCode,exacteJacExpr,fname,closurefunc) 
    elseif length(functionEx.args)==1
        #@warn("Defining a helper function inside the main model may hurt performance. Consider moving it to the top level.")
        closurefunc=@RuntimeGeneratedFunction(functionEx.args[1]) 
        diffEqfunction=createContEqFun(otherCode,equs,fname,closurefunc)# diff equations before this are stored in a dict:: now we have a giant function that holds all diff equations
        exactJacfunction=createExactJacFun(otherCode,exacteJacExpr,fname,closurefunc)
    else
        error("Error: Currently only one helper function is allowed inside the model function. Consider moving them to the top level.")
    end
  
    if numHelperFunCalls>length(functionEx.args)# the case when the user has defined a helper function inside the model function and one helperF outside and only called the outside helperF is considered a user mistake (either user made a typo or forgot to remove).
            if is_top_level                    
                RuntimeGeneratedFunctions.init(Main)
                exactJacfunctionF=RuntimeGeneratedFunction(mod, mod,exactJacfunction)
                diffEqfunctionF=RuntimeGeneratedFunction(mod, mod,diffEqfunction)   
            else
              @warn("Simulation is not running on same level with helper functions. This may hurt performance. Consider placing the solve function on top level, or placing the helper function inside the model.")
             exactJacfunctionF1 = Base.eval(mod, exactJacfunction)
             exactJacfunctionF = (args...) -> Base.invokelatest(exactJacfunctionF1, args...)
             diffEqfunctionF1 = Base.eval(mod, diffEqfunction)
             diffEqfunctionF = (args...) -> Base.invokelatest(diffEqfunctionF1, args...)
            end
    else
        exactJacfunctionF=@RuntimeGeneratedFunction(exactJacfunction)
        diffEqfunctionF=@RuntimeGeneratedFunction(diffEqfunction) 
    end

    prob=NLODEContProblemSpan(fname,Val(1),Val(T),Val(D),Val(0),Val(num_cache_equs),initConditions,collect(discrVars),diffEqfunctionF,jacVect,SDVect,exactJacfunctionF,tspan,[closurefunc])# prtype type 1...prob not saved and struct contains vects
end
 
