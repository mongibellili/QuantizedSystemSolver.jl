using QuantizedSystemSolver
ex=:(u < 1 || u > 10)
flattened = Expr[]
QuantizedSystemSolver.recurse(ex, flattened)
@show flattened #Expr[:(u < 1), :(u > 10)] 

ex=:(u < 1 || u > 10)
(kind,flattened)=QuantizedSystemSolver.decompose_condition(ex)
@show (kind,flattened) #(:or, Expr[:(u < 1), :(u > 10)])

ex1 = :(u > 10)
ex2 = :(u < 1)
zcf1 = QuantizedSystemSolver.to_zcf(ex1)
zcf2 = QuantizedSystemSolver.to_zcf(ex2)
(zcf1, zcf2)
@show (zcf1, zcf2) #(:(u - 10), :(1 - u))

(ex, stateVarName, discrParamName,muteVar, param) = (:(du[k] = u[k] * u[k - 1] * coef2), :u,:p, :k, Dict{Symbol, Union{Float64, Int64,Expr,Symbol}}(:coef1 => 2.0, :coef2 => 1.5))

  newEx=QuantizedSystemSolver.changeVarNames_params(ex, stateVarName,discrParamName, muteVar, param,Set([:f]))
@show (newEx, stateVarName, muteVar, param)#(:(du[i] = q[i] * q[i - 1] * 1.5), :u, :k, Dict{Symbol, Union{Float64, Int64, Expr, Symbol}}(:coef1 => 2.0, :coef2 => 1.5))