module SimpleModelIR
using MacroTools: postwalk 
include("types.jl")
include("build_ir.jl")
include("normalize_ir.jl")
export
    AbstractODEStatement,
    AssignStatement,
    IfStatement,
    ForStatement,
    WhileStatement,
    ExprStatement,
    ODEFunctionIR,
    problem_to_normalized_ir,
    build_ir,
    normalize_ir,
    recurse,
    decompose_condition,
    changeVarNames_params,
    to_zcf,
    process_if_condition,
    process_if_block,
    probInfo,
    process_if_expr
end
