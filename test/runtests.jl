using QuantizedSystemSolver
using Test
using Plots
using BSON
@testset "QuantizedSystemSolver.jl" begin
    include("./unit/taylor0Functions.jl")
    include("./unit/taylor0Arithmetic.jl")
    include("./unit/qssUnitTests.jl")
    include("./integration/exampleTest.jl")
    include("./integration/baseFunctions.jl")
    include("./integration/realSystemsTests.jl")
    include("./integration/helperFunctionUsage.jl")
    include("./integration/compositeConditionTest.jl")
end
