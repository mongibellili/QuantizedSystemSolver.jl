
"""
    NLODEProblem{F,PRTYPE,T,D,Z,CS}
This is a superclass for all ODE problems. It is parametric on:\n
  - The problem type PRTYPE.
  - The number of continuous variables T
  - The number of discrete events D
  - The number of events (zero crossing functions) Z
  - The cache size CS.
"""
abstract type NLODEProblem{F,PRTYPE,T,D,Z,CS} end

"""
    Algorithm{N,O}
This is a superclass for all QSS algorithms. It is parametric on:\n
  - The name of the algorithm N
  - The order of the algorithm O
"""
abstract type Algorithm{N,O} end

"""
    Sol{T,O}
This is a superclass for all QSS solutions. It is parametric on:\n
  - The number of continuous variables T
  - The order of the algorithm O
"""
abstract type Sol{T,O} end







