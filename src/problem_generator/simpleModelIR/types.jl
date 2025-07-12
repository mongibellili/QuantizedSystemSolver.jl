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
mutable  struct AssignStatement <: AbstractODEStatement
    lhs::Any
    rhs::Any
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
    condition::Expr
    body::Expr  # Previously Vector{AbstractODEStatement}
end

"""
    ForStatement

Represents a `for` loop statement within the intermediate representation (IR) of an ODE problem.



# Description
This mutable struct is a subtype of `AbstractODEStatement` and is used to represent a `for` loop that is used to define differential equations.

# Example
"""
mutable  struct ForStatement <: AbstractODEStatement
    var::Symbol
    start::Any
    stop::Any
    body::Vector{AbstractODEStatement}
    statement::Any
    loop_iter::Any  # raw expression (optional)
    loop_type::Symbol     # :range, :range_with_statement, :array
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



# A mutable struct representing the intermediate representation (IR) of an ordinary differential equation (ODE) function. This type contains a vector that stores all previous IR statements. 
mutable struct ODEFunctionIR
    statements::Vector{AbstractODEStatement}
end