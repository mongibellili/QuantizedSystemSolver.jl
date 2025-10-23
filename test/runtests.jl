using QuantizedSystemSolver
using Test

#using RuntimeGeneratedFunctions
@testset "QuantizedSystemSolver.jl" begin
    include("./unit/all_unit_tests.jl")
    include("./integration/all_integration_tests.jl")

end
