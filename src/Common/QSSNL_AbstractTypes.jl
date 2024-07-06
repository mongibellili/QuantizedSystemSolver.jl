
"""NLODEProblem{PRTYPE,T,Z,Y,CS}
This is superclass for all NLODE problems. It is parametric on these:\n
    - The problem type PRTYPE.\n
    - The number of continuous variables T\n
    - The number of events (zero crossing functions) Z\n
    - The actual number of events (an if-else statment has one zero crossing functions and two events) Y\n
    - The cache size CS.
"""
abstract type NLODEProblem{PRTYPE,T,Z,Y,CS} end

"""ALGORITHM{N,O}
This is superclass for all QSS algorithms. It is parametric on these:\n
    - The name of the algorithm N\n
    - The order of the algorithm O
"""
abstract type ALGORITHM{N,O} end

"""Sol{T,O}
This is superclass for all QSS solutions. It is parametric on these:\n
    - The number of continuous variables T\n
    - The order of the algorithm O
"""
abstract type Sol{T,O} end







