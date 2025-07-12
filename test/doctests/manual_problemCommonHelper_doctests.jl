using QuantizedSystemSolver

ex=:(q[i - 1]) 
newEx=QuantizedSystemSolver.changeExprToFirstValue(ex)
@show newEx  #  :((q[i - 1])[0]) 

ex=:(i - 1) 
newEx=QuantizedSystemSolver.symbolFromRef(:q,ex)
@show (ex,newEx)  #(ex, newEx) = (:(i - 1), :qiminus1)

symDict= Dict{Symbol, Expr}(:qi => :(q[i]), :q10 => :(q[10]), :q2 => :(q[2]), :qiminus1 => :(q[i - 1]), :q1 => :(q[1]))
coefExpr=:(1.5qiminus1) 
newEx=QuantizedSystemSolver.restoreRef(coefExpr, symDict)
@show newEx  # newEx = :(1.5 * (q[i - 1])[0])

(varNum, rhs, jac, exactJacExpr, symDict, dD) = (1, :(p[2] - 2.0 * q[1] * p[2]), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(), Dict{Expr, Union{Float64, Int64, Expr, Symbol}}(), Dict{Symbol, Expr}(:q10 => :(q[10]), :p2 => :(p[2]), :qiminus1 => :(q[i - 1]), :p1 => :(p[1]), :q1 => :(q[1])), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}())
QuantizedSystemSolver.extractJacDepNormal(varNum, rhs, jac, exactJacExpr,:symbolic, symDict, dD )
(jac, exactJacExpr, dD) 
@show (jac, exactJacExpr, dD)#(jac, exactJacExpr, dD) = (Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(1 => Set([1])), Dict{Expr, Union{Float64, Int64, Expr, Symbol}}(:((1, 1)) => :(-2.0 * p[2])), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(2 => Set([1])))

(b, niter, rhs, jac, exactJacExpr, symDict, dD) = (2, 9, :(p[1] * q[i - 1] * 1.5), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(1 => Set([1])), Dict{Expr, Union{Float64, Int64, Expr, Symbol}}(:((1, 1)) => :(-2.0 * p[2])), Dict{Symbol, Expr}(:q10 => :(q[10]), :p2 => :(p[2]), :qiminus1 => :(q[i - 1]), :p1 => :(p[1]), :q1 => :(q[1])), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(2 => Set([1])))
QuantizedSystemSolver.extractJacDepLoop(b, niter, rhs, jac, exactJacExpr,:symbolic, symDict, dD )
(jac, exactJacExpr, dD) 
@show (jac, exactJacExpr, dD)#(jac, exactJacExpr, dD) = (Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(:((2, 9)) => Set([:(i - 1)]), 1 => Set([1])), Dict{Expr, Union{Float64, Int64, Expr, Symbol}}(:((1, 1)) => :(-2.0 * p[2]), :(((2, 9), i - 1)) => :(1.5 * p[1])), Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(2 => Set([1]), 1 => Set([:((2, 9))])))

jac=Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(1 => Set([2, 1]),:((2, 9)) => Set([:(i - 1), :i]),10 => Set([1, 10]))
jacVect=QuantizedSystemSolver.createJacVect(jac,Val(10) )
@show string(jacVect)#string(jacVect) = "[[2, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [10, 1]]"

jac=Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(1 => Set([2, 1]),:((2, 9)) => Set([:(i - 1), :i]),10 => Set([1, 10]));
SD=QuantizedSystemSolver.createSDVect(jac,Val(10) );
@show string(SD)#string(SD) = "[[10, 2, 1], [2, 3, 1], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9], [10]]"

exacteJacExpr=Dict{Expr,Union{Float64,Int,Symbol,Expr}}(:((1, 1)) => :(-2.0 * (q[2])[0]), :((1, 2)) => :(1 - 2.0 * (q[1])[0]),:(((2, 9), i - 1)) => :((q[i])[0]), :(((2, 9), i)) => :((q[i - 1])[0]),:((10, 10)) => -1, :((10, 1)) => 1);
exactJac=QuantizedSystemSolver.createExactJacFun(:(),exacteJacExpr,:f,0);
@show exactJac

#= exactJac = :(function exactJacf(q::Vector{Taylor0}, p::Vector{Float64}, cache::AbstractVector{Float64}, i::Int, j::Int, t::Float64, f_)
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
  end) =#


  equs = Dict{Union{Int,Expr},Union{Int,Symbol,Expr}}(10 => :(subT(q[1], q[10], cache[1])), :((2, 9)) => :(mulT(q[i], q[i - 1], cache[1])), 1 => :(subT(q[2], mulTT(2.0, q[1], q[2], cache[2], cache[3]), cache[1])));
diffEqfun=QuantizedSystemSolver.createContEqFun(:(),equs,:f,0); 
@show diffEqfun

#= diffEqfun = :(function f(i::Int, q::Vector{Taylor0}, p::Vector{Float64}, t::Taylor0, cache::Vector{Taylor0}, f_)
      (if i == 0
              return nothing
          elseif i == 10
              subT(q[1], q[10], cache[1])
              return nothing
          elseif 2 <= i <= 9
              mulT(q[i], q[i - 1], cache[1])
              return nothing
          elseif i == 1
              subT(q[2], mulTT(2.0, q[1], q[2], cache[2], cache[3]), cache[1])
              return nothing
          end,)
  end) =#

ex=:(q[2])
newEx=QuantizedSystemSolver.transformFSimplecase(ex);
@show newEx # :(createT(q[2], cache[1]))

ex=:(q[2] - 2.0*q[1]*q[2],1)
newEx=QuantizedSystemSolver.transformF(ex);
@show newEx # :((subT(q[2], mulTT(2.0, q[1], q[2], cache[2], cache[3]), cache[1]), 3))