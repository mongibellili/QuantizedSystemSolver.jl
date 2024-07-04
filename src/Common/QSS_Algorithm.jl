#

struct QSSAlgorithm{N,O}#<:QSSAlgorithm{N,O}
    name::Val{N}
    order::Val{O}
end

"""qss1()
calls the explicit quantized state system solver with order 1
"""
qss1()=QSSAlgorithm(Val(:qss),Val(1))

"""qss2()
calls the explicit quantized state system solver with order 2
"""
qss2()=QSSAlgorithm(Val(:qss),Val(2))
"""qss3()
calls the explicit quantized state system solver with order 3
"""
qss3()=QSSAlgorithm(Val(:qss),Val(3))
"""nmliqss1()
calls the modified imlicit quantized state system solver with order 1
"""
nmliqss1()=QSSAlgorithm(Val(:nmliqss),Val(1))
"""nmliqss2()
calls the modified imlicit quantized state system solver with order 2
"""
nmliqss2()=QSSAlgorithm(Val(:nmliqss),Val(2))

"""nmliqss3()
calls the modified imlicit quantized state system solver with order 3
"""
nmliqss3()=QSSAlgorithm(Val(:nmliqss),Val(3))

#= nliqss1()=QSSAlgorithm(Val(:nliqss),Val(1))
nliqss2()=QSSAlgorithm(Val(:nliqss),Val(2))
nliqss3()=QSSAlgorithm(Val(:nliqss),Val(3))

mliqss1()=QSSAlgorithm(Val(:mliqss),Val(1))
mliqss2()=QSSAlgorithm(Val(:mliqss),Val(2))
mliqss3()=QSSAlgorithm(Val(:mliqss),Val(3)) =#
"""liqss1()
calls the  imlicit quantized state system solver with order 1
"""
liqss1()=QSSAlgorithm(Val(:liqss),Val(1))
"""liqss2()
calls the  imlicit quantized state system solver with order 2
"""
liqss2()=QSSAlgorithm(Val(:liqss),Val(2))
"""liqss3()
calls the  imlicit quantized state system solver with order 3
"""
liqss3()=QSSAlgorithm(Val(:liqss),Val(3))



#= sparse()=Val(true)
dense()=Val(false) =#

