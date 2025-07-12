#get .cov files
#= using Pkg
Pkg.test("QuantizedSystemSolver"; coverage=true) # run this then comment and run next code =#
#= using Pkg
Pkg.test("QuantizedSystemSolver"; coverage=false) =#

 # get lcov from .cov files (summary)
 using Coverage
coverage = process_folder()
open("lcov.info", "w") do io
    LCOV.write(io, coverage)
end; 

 #clean up the folder when not needed
#=  using Coverage
 Coverage.clean_folder(".") =#


# generate docs
# PS C:\Users\belli\.julia\dev\QuantizedSystemSolver\docs> julia --project=docs  make.jl