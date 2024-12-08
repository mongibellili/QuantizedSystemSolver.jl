
"""
    QSSAlgorithm{N,O}
This is the first subclass for QSS algorithms. It is parametric on:\n
  - The name of the algorithm N
  - The order of the algorithm O
"""
struct QSSAlgorithm{N,O}<:Algorithm{N,O}
    name::Val{N}
    order::Val{O}
end


function show(io::IO, a::Algorithm{N,O}) where {N,O}
 print(io, N,O)
end

"""
    qss1()
calls the explicit quantized state system solver with order 1
```julia
qss1()=QSSAlgorithm(Val(:qss),Val(1))
```
"""
qss1()=QSSAlgorithm(Val(:qss),Val(1))

"""
    qss2()
calls the explicit quantized state system solver with order 2
"""
qss2()=QSSAlgorithm(Val(:qss),Val(2))

"""
    nmliqss1()
calls the modified imlicit quantized state system solver with order 1.
It is efficient when the system contains large entries outside the main diagonal of the Jacobian .
"""
nmliqss1()=QSSAlgorithm(Val(:nmliqss),Val(1))
"""
    nmliqss2()
calls the modified imlicit quantized state system solver with order 2.
It is efficient when the system contains large entries outside the main diagonal of the Jacobian .
"""
nmliqss2()=QSSAlgorithm(Val(:nmliqss),Val(2))


"""
    liqss1()
calls the  imlicit quantized state system solver with order 1.
"""
liqss1()=QSSAlgorithm(Val(:liqss),Val(1))
"""
    liqss2()
calls the  imlicit quantized state system solver with order 2.
"""
liqss2()=QSSAlgorithm(Val(:liqss),Val(2))


