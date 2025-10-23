
"""
    ODEProblemData{JACMODE,T,D,Z,CS,F,JAC,CLS}
This is a superclass for all ODE problems. It is parametric on:\n
  - The problem type JACMODE.
  - The number of continuous variables T
  - The number of discrete events D
  - The number of events (zero crossing functions) Z
  - The cache size CS.
"""
abstract type ODEProblemData{JACMODE,T,D,Z,CS,F,JAC,CLS} end

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


"""
  LiQSS_Data{O,M}

An abstract type representing the data structure for the LiQSS (Linearly Quantized State System) solver.

# Type Parameters
- `O`: The type parameter representing the order of the QSS method.
- `M`: The cycle detection mechanism number, which can be used to specify different cycle detection strategies.
"""
abstract type LiQSS_Data{O,M} end




