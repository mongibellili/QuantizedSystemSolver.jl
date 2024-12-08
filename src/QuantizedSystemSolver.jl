module QuantizedSystemSolver
  const global VERBOSE=false   # later move to solve to allow user to use it.
  using RuntimeGeneratedFunctions
  using StaticArrays
  using SymEngine
  using ExprTools  #combineddef
  using MacroTools: postwalk,prewalk, @capture
  using CodeTracking
  using Plots: savefig
  using Dates: now,year,month,day,hour,minute,second #fortimestamp
  RuntimeGeneratedFunctions.init(@__MODULE__)
  import Plots: plot!,plot
          ##### for taylorseries subcomponent   #######
  import Base: ==, +, -, *, /, ^             
  import Base: iterate, size, eachindex, firstindex, lastindex,
     length, getindex, setindex!
  import Base:  sqrt, exp, log, sin, cos, sincos, tan,
    asin, acos, atan, sinh, cosh, tanh, atanh, asinh, acosh,
    zero, one, zeros, ones, isinf, isnan, iszero,
    convert,  show,abs
                ##### list of public (API) 
  export ODEProblem,@NLodeProblem,NLodeProblem,solve ,NLODEProblem# 
  export qss1,qss2,qss3,liqss1,liqss2,liqss3,saveat,nmliqss1,nmliqss2,nmliqss3
  export save_Sol,plot_Sol,getPlot,getPlot!,save_SolSum,solInterpolated,plot_SolSum,plot
  export getErrorByRefs,getAverageErrorByRefs,getError,getAverageError
  # public functions and structs used in documentation
  export Taylor0,mulT,mulTT,createT,addsub,negateT,subsub,subadd,subT,addT,muladdT,mulsub,divT,powerT,testTaylor # for CI testing
  export QSSAlgorithm,Sol,EventDependencyStruct,Stats,CommonQSS_Data,LiQSS_Data   
  #####  Taylor series  #########
  include("taylor/constructors.jl") 
  include("taylor/arithmetic.jl")
  include("taylor/arithmeticT.jl")
  include("taylor/functions.jl")
  include("taylor/functionsT.jl")
  include("taylor/power.jl")
  include("taylor/powerT.jl")
  
  ##### commonQSS #########           
  include("commonQSS/qssAbstractTypes.jl")
  include("commonQSS/rootFinders.jl") 
  include("commonQSS/qssAlgorithm.jl")
  include("commonQSS/qssData.jl")
  include("commonQSS/scheduler.jl")
  ##### solution #########
  include("solution/solution.jl")
  include("solution/solutionPlot.jl")
  include("solution/solutionError.jl")
  #####  problem #########
  include("problem/taylorEquationConstruction.jl")
  include("problem/qssProblemContinuousHelper.jl")
  include("problem/qssProblemDiscreteHelper.jl")
  include("problem/qssProblemDefinition.jl")
  include("problem/qssProblemContinuous.jl")
  include("problem/qssProblemDiscrete.jl")
                ##### integrators  ###########
  include("integrators/qssIntegrator.jl")
  include("integrators/qssDiscreteIntegrator.jl")
  # implicit integrator when large entries on the main diagonal of the jacobian
  include("integrators/liqssIntegrator.jl")
  include("integrators/liqssDiscreteIntegrator.jl")
  # implicit integrator when large entries NOT on the main diagonal of the jacobian
  include("integrators/nmliqssIntegrator.jl")
  include("integrators/nmliqssDiscreteIntegrator.jl")
                ##### Quantizers ######
  # commonQSS & explicit
  include("quantizers/quantizerQss.jl")
  # implicit & single updateQ
  include("quantizers/quantizerLiqss1.jl")
  include("quantizers/quantizerLiqss2.jl")  
  # cycle detection and simultaneous update
  include("quantizers/quantizerMliqss1.jl")
  include("quantizers/quantizerMliqss2.jl")

               ##### main entrance/ Interface #######
  include("interface/main.jl")
  include("interface/macroInterface.jl")
  include("interface/solve.jl")
end # module



