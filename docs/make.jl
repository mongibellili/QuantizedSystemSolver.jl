using QuantizedSystemSolver
using StaticArrays
using Documenter

# Function to deploy documentation

    makedocs(
        sitename = "Quantized System Solver",
        modules = [QuantizedSystemSolver],
        pages = [       
        
        "Home" => "index.md",
        "Getting started" => [
            "guide/userTutorial.md",
            "guide/userAPI.md",],
        
        "Background" => [
            # Theoretical background
            "background/introductoryResources.md",
            "background/qss.md",
            "background/liqss.md",
            "background/mliqss.md",
        ],
        "Examples" => [
            "examples/linearTimeInvariantSystem.md",
            "examples/vanderpol.md",
            "examples/tysonModel.md",
            "examples/adr.md",
            "examples/dosing.md",
        ],
        "Developer resources" => [
            "developer/devIntro.md",
            "developer/problem.md",
            "developer/algorithm.md",
            "developer/solve.md",
            "developer/integrator.md",
            "developer/quantizer.md",
            "developer/solution.md",
            "developer/taylor.md",
            "developer/utils.md",
        ],
      
          
        ]
    )

    deploydocs(
      repo   = "github.com/mongibellili/QuantizedSystemSolver.jl.git",
      branch = "gh-pages",
    
    )
  