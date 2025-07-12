module QuantizedSystemSolver
  using RuntimeGeneratedFunctions
  using StaticArrays
  using SymEngine
  using ExprTools  #combineddef
  using MacroTools: postwalk,prewalk, @capture
  using CodeTracking
  using RecipesBase
  using Dates: now,year,month,day,hour,minute,second #fortimestamp
  RuntimeGeneratedFunctions.init(@__MODULE__)
 
  import Plots: plot!,plot
          ##### for taylorseries subcomponent   #######
  import Base: ==, +, -, *, /, ^             
  import Base: iterate, size, eachindex, firstindex, lastindex,
     length, getindex, setindex!
  import Base:  sqrt, exp, log, sin, cos, sincos, tan,
    asin, acos, atan, sinh, cosh, tanh, atanh, asinh, acosh,
    zero, one, zeros, ones, isinf, isnan, iszero,sign,
    convert,  show,abs,mod,rem,min,max
                ##### list of public (API) 
  export ODEProblemTest,ODEProblem,solve ,ODEProblemData,Detection# 
  export qss1,qss2,qss3,liqss1,liqss2,liqss3,saveat,nmliqss1,nmliqss2,nmliqss3
  export save_SolSum,solInterpolated,plot_SolSum,plot
  export getErrorByRefs,getAverageErrorByRefs,getError,getAverageError
  # public functions and structs used in documentation
  export Taylor0,mulT,mulTT,createT,addsub,negateT,subsub,subadd,subT,addT,muladdT,mulsub,divT,powerT,testTaylor,Derivative # for CI testing
  export QSSAlgorithm,Sol,EventDependencyStruct,Stats,CommonQSS_Data,LiQSS_Data ,ODEFunctionIR,PreProcessData

  include("qssAbstractTypes.jl")
                                          ############ problem generator ############
  ##### SimpleModelIR                                        
  include("problem_generator/simpleModelIR/SimpleModelIR.jl")
  using .SimpleModelIR
  #####  Taylor series  #########
  include("problem_generator/taylor/constructors.jl") 
  include("problem_generator/taylor/arithmetic.jl")
  include("problem_generator/taylor/arithmeticT.jl")
  include("problem_generator/taylor/functions.jl")
  include("problem_generator/taylor/functionsT.jl")
   #####  problem #########
  include("problem_generator/problem/taylorEquationConstruction.jl")
  include("problem_generator/problem/qssProblemCommonHelper.jl")
  include("problem_generator/problem/qssProblemContinuousHelper.jl")
  include("problem_generator/problem/qssProblemDiscreteHelper.jl")
  include("problem_generator/problem/qssProblemDefinition.jl")
  include("problem_generator/problem/qssProblemContinuous.jl")
  include("problem_generator/problem/qssProblemDiscrete.jl")
  
                                        ############## solution ################
  include("solution/solution.jl")
  include("solution/solutionPlot.jl")
  include("solution/solutionError.jl")

                                      ############## simulator ################
  ##### commonQSS #########           
  include("simulator/commonQSS/rootfinders.jl") 
  include("simulator/commonQSS/qssAlgorithm.jl")
  include("simulator/commonQSS/qssData.jl")
  include("simulator/commonQSS/scheduler.jl")
  ##### integrators  ###########
  include("simulator/integrators/qssIntegrator.jl")
  include("simulator/integrators/qssDiscreteIntegrator.jl")
  include("simulator/integrators/liqssIntegrator.jl")  # implicit integrator when large entries on the main diagonal of the jacobian
  include("simulator/integrators/liqssDiscreteIntegrator.jl")
  include("simulator/integrators/nmliqssIntegrator.jl")  # implicit integrator when large entries NOT on the main diagonal of the jacobian
  include("simulator/integrators/nmliqssDiscreteIntegrator.jl")
   ##### Quantizers ######
  include("simulator/quantizers/quantizerQss.jl")    # common and explicit updateQ
  include("simulator/quantizers/quantizerLiqss1.jl")  # implicit & single updateQ order 1
  include("simulator/quantizers/quantizerLiqss2.jl")  # implicit & single updateQ order 2
  include("simulator/quantizers/cycleDetectionConditions.jl")  # cycle detection
  include("simulator/quantizers/quantizerMliqss1.jl")   #  simultaneous update order 1
  include("simulator/quantizers/quantizerMliqss2.jl")     #  simultaneous update order 2

                                       ##### main entrance/ Interface #######
  include("interface/main.jl")   # initiate problem generator
  include("interface/solve.jl")  # initiate simulator
end # module



