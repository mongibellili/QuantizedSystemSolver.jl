module SimpleModelIR
const SKIP_SYMBOLS = (:+, :-, :*, :/, :^, :%, :&, :|, :!, :(=), :(==), :!=, :<, :>, :<=, :>=)
@enum InlineMode MANUAL AUTO FULL
using MacroTools: postwalk 
import Base: string
export  AbstractODEStatement, AssignStatement, IfStatement, ForStatement, WhileStatement, ExprStatement, ODEFunctionIR, problem_to_normalized_ir, build_ir, normalize_ir
export changeExprToFirstValue # handle_events and jac_expression use this

# these expoprts are for testing
export recurse, decompose_condition, changeVarNames_params, to_zcf, process_if_condition, process_if_block, probInfo, process_if_expr

include("ir_types.jl")
include("symbolTableFunctions.jl")
include("build_ir.jl")
include("normalize_ir_ifstatement.jl")
include("normalize_ir.jl")
end