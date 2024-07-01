using Documenter
using QuantizedSystemSolver
makedocs(sitename="Quantized State System Solver",
    modules=[QuantizedSystemSolver],
    pages=[
        "Home" => "index.md"
    ])

    deploydocs(
        repo = "github.com/mongibellili/mongibellili.github.io.git",
        target = "gh-pages",  # Ensure this branch exists or is created
        branch = "gh-pages"
    )
