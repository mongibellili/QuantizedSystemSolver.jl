using QuantizedSystemSolver

(counter, zcf, zcjac, SZ, dZ) = (2, :(q[2] - p[1]), [[1]], Dict{Int64, Set{Int64}}(1 => Set([1])), Dict{Int64, Set{Int64}}())
QuantizedSystemSolver.extractZCJacDep(counter, zcf, zcjac, SZ, dZ)
(zcjac, SZ, dZ) 
@show (zcjac, SZ, dZ) # (zcjac, SZ, dZ) = ([[1], [2]], Dict{Int64, Set{Int64}}(2 => Set([2]), 1 => Set([1])), Dict{Int64, Set{Int64}}(1 => Set([2])))


(SZ, T) = (Dict{Int64, Set{Int64}}(2 => Set([2]), 1 => Set([1])), 10)
szVec=QuantizedSystemSolver.createSZVect(SZ, Val(T))
@show string(szVec)  #"[[1], [2], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[]]"

(dD, D) = (Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(2 => Set([10, 1]), 1 => Set([:((2, 9))])), 2)
dDVect=QuantizedSystemSolver.createdDVect(dD, Val(D) )
@show string(dDVect) # "[[2, 3, 4, 5, 6, 7, 8, 9], [10, 1]]"


(dD, dZ, eventDep) = ([[2, 3, 4, 5, 6, 7, 8, 9], [10, 1]], Dict{Int64, Set{Int64}}(1 => Set([2])), QuantizedSystemSolver.EventDependencyStruct[QuantizedSystemSolver.EventDependencyStruct(1, Int64[], [1], Int64[]), QuantizedSystemSolver.EventDependencyStruct(2, Int64[], Int64[], Int64[]), QuantizedSystemSolver.EventDependencyStruct(3, [3], [2], [3, 1, 2]), QuantizedSystemSolver.EventDependencyStruct(4, Int64[], Int64[], Int64[])])
(HZ1, HD1) =QuantizedSystemSolver.createDependencyToEventsDiscr(dD, dZ, eventDep )
(HZ1, HD1) 
@show (HZ1, HD1) # ([[2], Int64[], Int64[], Int64[]], [[5, 4, 6, 7, 2, 9, 8, 3], Int64[], [10, 1], Int64[]])


(SD, sZ, eventDep) = ([[10, 2, 1], [3], [4], [5], [6], [7], [8], [9], Int64[], [10]], Dict{Int64, Set{Int64}}(2 => Set([2]), 1 => Set([1])), QuantizedSystemSolver.EventDependencyStruct[QuantizedSystemSolver.EventDependencyStruct(1, Int64[], [1], Int64[]), QuantizedSystemSolver.EventDependencyStruct(2, Int64[], Int64[], Int64[]), QuantizedSystemSolver.EventDependencyStruct(3, [3], [2], [3, 1, 2]), QuantizedSystemSolver.EventDependencyStruct(4, Int64[], Int64[], Int64[])])
(HZ2, HD2) =QuantizedSystemSolver.createDependencyToEventsCont(SD, sZ, eventDep)
(HZ2, HD2) 
@show (HZ2, HD2) # ([Int64[], Int64[], Int64[], Int64[]], [Int64[], Int64[], [4], Int64[]])

(HD1, HD2) = ([[5, 4, 6, 7, 2, 9, 8, 3], Int64[], [10, 1], Int64[]], [Int64[], Int64[], [4], Int64[]])
HD=QuantizedSystemSolver.unionDependency(HD1, HD2)
@show string(HD) #"[[5, 4, 6, 7, 2, 9, 8, 3], Int64[], [4, 10, 1], Int64[]]"