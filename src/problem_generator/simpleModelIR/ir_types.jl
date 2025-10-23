# ================================
# Structs for Intermediate IR
# ================================

"""
    AbstractODEStatement

An abstract type representing a single statement in the user code.
This serves as a base type for various specific types of statements that can be part of an ODE function's intermediate representation (IR).
It is used to encapsulate different kinds of operations, such as assignments, conditional statements, loops, and expressions, allowing for a structured representation of the ODE function's logic.
"""
abstract type AbstractODEStatement end

"""
    AssignStatement 

Represents an assignment statement within an ODE statement in the intermediate representation (IR) of a simple model.
This mutable struct is used to store information about lhs and rhs parts of an assignment.


"""
mutable struct AssignStatement <: AbstractODEStatement
    lhs::Union{Symbol, Expr}   # Expr(:ref, ...) allowed
    rhs::Union{Expr, Symbol, Number}
    keep_assignment::Bool # If true, the assignment is kept in the final code; if false, it is optimized out.
end
# Default constructor with keep_assignment = true
function AssignStatement(lhs, rhs) 
    if !(lhs isa Symbol || lhs isa Expr)
        throw(TypeError(:AssignStatement,"Invalid lhs type: Currently lhs must be Symbol or Expr",lhs))
    end
    if !(rhs isa Expr || rhs isa Symbol || rhs isa Number)
        throw(TypeError(:AssignStatement,"Invalid rhs type: Currently rhs must be Expr, Symbol, or Number",rhs))
    end
    AssignStatement(lhs, rhs, true)
end

"""
    IfStatement 

Represents a conditional (if) statement within the intermediate representation (IR) of a simple model ODE system.

# Fields
- `condition`: The condition expression to evaluate. it represents the zero-crossing function for an event.
- `body`: The statements to execute if the condition is true. it contains an expression of the whole if-statment. It will be used as the actual execution of the event.

# Usage
Used to model control flow in the IR for ODE problem generation.
"""
mutable struct IfStatement <: AbstractODEStatement
    condition::Union{Symbol, Expr}
    body::Expr  # Previously Vector{AbstractODEStatement}: inside if (1) no taylor vars (2) minimal code removal: done manually in process_if_block function.
end

"""
    ForStatement

Represents a `for` loop statement within the intermediate representation (IR) of an ODE problem.



# Description
This mutable struct is a subtype of `AbstractODEStatement` and is used to represent a `for` loop that is used to define differential equations.

# Example
"""
mutable struct ForStatement <: AbstractODEStatement
    var::Symbol
    start::Union{Int, Float64, Symbol, Expr}
    stop::Union{Int, Float64, Symbol, Expr}
    body::Vector{AbstractODEStatement}
    statement::Nothing                      # maybe future use
    loop_iter::Union{Nothing, Expr}     # full `k:N` or similar
    loop_type::Symbol                   # :range, :array, etc.
end

mutable  struct WhileStatement <: AbstractODEStatement
    condition::Expr
    body::Vector{AbstractODEStatement}
end

"""
    ExprStatement 

A mutable struct representing an expression statement within the ODE problem intermediate representation (IR).


# Description
`ExprStatement` is used to encapsulate a single expression that was not handled by other specific statement types like `AssignStatement`, `IfStatement`, `ForStatement`, or `WhileStatement`.
It allows for the inclusion of arbitrary expressions in the IR, which can be useful for representing complex operations or computations that do not fit neatly into the other categories.
"""
mutable  struct ExprStatement <: AbstractODEStatement
    expr::Expr
end


_code(x) = sprint(show, x)   # Works for Expr, Symbol, Number, etc.

Base.string(a::AssignStatement) = string(_code(a.lhs), " = ", _code(a.rhs), a.keep_assignment ? "" : " (optimized out)")
Base.string(a::IfStatement)     = string("if ", _code(a.condition), " ", _code(a.body))
Base.string(a::ForStatement)    = string("for ", _code(a.var), " in ", _code(a.start), ":", _code(a.stop), " ", _code(a.body))
Base.string(a::WhileStatement)  = string("while ", _code(a.condition), " ", _code(a.body))
Base.string(a::ExprStatement)   = _code(a.expr)




# A mutable struct representing the intermediate representation (IR) of an ordinary differential equation (ODE) function. This type contains a vector that stores all previous IR statements. 
mutable struct ODEFunctionIR
    statements::Vector{AbstractODEStatement}
end
Base.show(io::IO, ir::ODEFunctionIR) = show(io, ir.statements)
function Base.show(io::IO, irVect::Vector{AbstractODEStatement})
    println(io,"")
    for (i, stmt) in enumerate(irVect)
        print(io, string(stmt))
        if i < length(irVect)
            print(io, "\n")  # newline between statements
        end
    end
end
#A helper struct that holds the final IR. in addition, it used to return the number of `if_statements` and the number of helper functions written by the user. These two numbers are extracted during the normalization
struct probInfo
    ir::ODEFunctionIR
    numZC::Int
    helperFunSymSet::Int64
end


# stack of tables of symbols
mutable struct SymbolEntry
    value::Union{Float64, Int64, Expr, Symbol, Nothing}  # Nothing if not safe to inline
    safe_to_inline::Bool
end

mutable struct SymbolTable
    entries::Dict{Symbol, SymbolEntry}  #     #if in future: allow @inline with lhs is expr, then entries::Dict{Union{Expr, Symbol}, SymbolEntry} 
end
SymbolTable() = SymbolTable(Dict{Symbol, SymbolEntry}())


struct SymbolTableStack  # for different scopes: A stack of SymbolTables to manage symbol scopes
    tables::Vector{SymbolTable}
end
SymbolTableStack() = SymbolTableStack([SymbolTable()])  # Initialize with a single empty SymbolTable


struct AssignAction
    register::Bool   # call add_symbol!
    swap::Bool       # allow rhs to be substituted in later expressions
    keep::Bool       # keep assignment in normalized code
    warn::Bool       # whether to emit a warning (e.g., custom function inlined)
end

RegisterNoSwapKeep()     = AssignAction(true,  false, true, false)
RegisterSwapKeep()       = AssignAction(true,  true,  true, false)
RegisterSwapRemove()     = AssignAction(true,  true,  false, false)
RegisterSwapRemoveWarn() = AssignAction(true,  true,  false, true)
NoRegisterKeep()         = AssignAction(false, false, true, false)