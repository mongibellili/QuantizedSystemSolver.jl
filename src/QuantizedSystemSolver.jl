module QuantizedSystemSolver
  const global VERBOSE=false   # later move to solve to allow user to use it.
  const global DEBUG=false
  using RuntimeGeneratedFunctions
  using StaticArrays
  using SymEngine
  using PolynomialRoots
  using ExprTools  #combineddef
  using MacroTools: postwalk,prewalk, @capture#, isexpr,
  using Plots: savefig
  using Dates: now,year,month,day,hour,minute,second #fortimestamp
  RuntimeGeneratedFunctions.init(@__MODULE__)
  import Plots: plot!,plot
   ##### this section belongs to taylorseries subcomponent
  import Base: ==, +, -, *, /, ^                
  #import Base: Base.gc_enable
  import Base: iterate, size, eachindex, firstindex, lastindex,
     length, getindex, setindex!
  import Base:  sqrt, exp, log, sin, cos, sincos, tan,
    asin, acos, atan, sinh, cosh, tanh, atanh, asinh, acosh,
    zero, one, zeros, ones, isinf, isnan, iszero,
    convert, promote_rule, promote, show,abs,show
   ##### list of public (API) 
  export NLodeProblem,solve ,NLODEProblem,QSSAlgorithm,Sol# docs used NLODEProblem, QSSAlgorithm
  export qss1,qss2,qss3,liqss1,liqss2,liqss3,saveat,nmliqss1,nmliqss2,nmliqss3
  export save_Sol,plot_Sol,getPlot,getPlot!,save_SolSum,solInterpolated,plot_SolSum
  export getError,getErrorByRodas,getAllErrorsByRefs,getAverageErrorByRefs,getErrorByRefs
  export Taylor0,mulT,mulTT,createT,addsub,negateT,subsub,subadd,subT,addT,muladdT,mulsub,divT,powerT,constructIntrval,getQfromAsymptote,iterationH,testTaylor # for testing
   ##### include section of Taylor series subcomponent
  include("ownTaylor/constructors.jl") 
  include("ownTaylor/arithmetic.jl")
  include("ownTaylor/arithmeticT.jl")
  include("ownTaylor/functions.jl")
  include("ownTaylor/functionsT.jl")
  include("ownTaylor/power.jl")
  include("ownTaylor/powerT.jl")
   ##### Utils
  include("Utils/rootfinders/SimUtils.jl") 
   ##### Common
  include("Common/TaylorEquationConstruction.jl")
  include("Common/QSSNL_AbstractTypes.jl")
  include("Common/Solution.jl")
  include("Common/SolutionPlot.jl")
  include("Common/SolutionError.jl")
  include("Common/Helper_QSSNLProblem.jl")
  include("Common/Helper_QSSNLDiscreteProblem.jl")
  include("Common/QSSNLContinousProblem.jl")
  include("Common/QSSNLdiscrProblem.jl")
  include("Common/QSS_Algorithm.jl")
  include("Common/QSS_data.jl")
  include("Common/Scheduler.jl")
  ##### integrators
  include("dense/NL_integrators/NL_QSS_Integrator.jl")
  include("dense/NL_integrators/NL_QSS_discreteIntegrator.jl")
  # implicit integrator when large entries on the main diagonal of the jacobian
  include("dense/NL_integrators/NL_LiQSS_Integrator.jl")
    include("dense/NL_integrators/NL_LiQSS_discreteIntegrator.jl")
  # implicit integrator when large entries NOT on the main diagonal of the jacobian
  include("dense/NL_integrators/NL_nmLiQSS_Integrator.jl")
  include("dense/NL_integrators/NL_nmLiQSS_discreteIntegrator.jl")
  ##### Quantizers
  include("dense/Quantizers/Quantizer_Common.jl")
  include("dense/Quantizers/QSS_quantizer.jl")
  include("dense/Quantizers/LiQSS_quantizer1.jl")
  include("dense/Quantizers/LiQSS_quantizer2.jl")
  #include("dense/Quantizers/LiQSS_quantizer3.jl")
  include("dense/Quantizers/mLiQSS_quantizer1.jl")
  include("dense/Quantizers/mLiQSS_quantizer2.jl")
  #include("dense/Quantizers/mLiQSS_quantizer3.jl")
  ##### main entrance/ Interface
  include("Interface/indexMacro.jl")
  include("Interface/QSS_Solve.jl")
end # module



