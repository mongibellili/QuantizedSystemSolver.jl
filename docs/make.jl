#import Pkg; Pkg.add("Documenter");Pkg.develop(url="https://github.com/mongibellili/QuantizedSystemSolver.git")
using Documenter
using QuantizedSystemSolver

# Set up deployment environment variables
#const DOCUMENTER_KEY = ENV["DOCUMENTER_KEY"]  # Load the private key from environment variables

# Function to deploy documentation

    makedocs(
        sitename = "Quantized System Solver",
        modules = [QuantizedSystemSolver],
        pages = [
            "Home" => "index.md",
            "Tutorial" => "userGuide.md",
            "developper Guide" => "developperGuide.md",
            "Examples" => "examples.md",
        ]
    )

    deploydocs(
      repo   = "github.com/mongibellili/QuantizedSystemSolver.git",
      branch = "gh-pages",
      #  deploy_key = ENV["DOCUMENTER_KEY"]  # Provide the deployment key explicitly
    )
      # Perform deployment using SSH key authentication
  #  run(`git config --global user.email "bellilimongi@gmail.com"`)
  #  run(`git config --global user.name "mongibellili"`)

  
   # cd("C:/Users/belli/Documents/mongibellili.github.io")
   # run(`pwd`)
   # cd("C:/Users/belli/.julia/dev/QuantizedSystemSolver/docs/build/")
   # run(`pwd`)
    # Copy the generated documentation files to the repository
  #  run(`cp -r "C:/Users/belli/.julia/dev/QuantizedSystemSolver/docs/build/*"`)  # Adjust the path based on your actual build directory

    # Add, commit, and push the changes
   # run(`git add .`)
   # run(`git commit -m "Deploy documentation"`)
   # run(`git push origin gh-pages`)  # Push to the gh-pages branch (or your designated branch)





# Main execution block
