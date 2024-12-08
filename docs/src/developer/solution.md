# Solution 
This is where the results of the simulation are stored, including:
Settings, Time points (times\[\]), Variables (vars\[\]), Statistics
(stats).

 
## Problem extension
The Solution entity can be extended by subclassing:
```@docs
Sol{T,O}
```

### What is needed with a new Solution:
Extending a solution is easier than extending the Problem and the Algorithm. Some cases are presented as follows:
  - Creating a new feature: no extra work is needed except the feature itself.
  - Creating a new solution similar to the LightSol: no extra work is needed except the definition of the new type subclassing Sol.
  - Creating a solution that differs from the LightSol: extend all functions that refer to fields that do not exist in the new type.


### Example
```julia
struct HeavySol{T,O}<:Sol{T,O}
  size::Val{T}
  order::Val{O}
  savedTimes::Vector{Vector{Float64}}
  savedVars::Vector{Vector{Taylor0}}
end
```
Here, the heavy solution contains all variables and their derivatives. The `evaluateSol` function can be reimplemented to use the derivatives for interpolation.

## Define and create a solution
An implemented concrete solution type is: 
```@docs
QuantizedSystemSolver.LightSol{T,O}
```
```@docs
QuantizedSystemSolver.createSol(::Val{T}, ::Val{O}, savedTimes::Vector{Vector{Float64}}, savedVars::Vector{Vector{Float64}}, solver::String, nameof_F::String, absQ::Float64, stats::Stats, ft::Float64) where {T,O}
```


## Further reading

**The solution Struct** The LightSol struct is specifically designed to
hold the solution of a system of ODEs. It contains several fields,
including size and order, which denote the number of continuous
variables and the order of the algorithm used, respectively. The
savedTimes and savedVars fields are vectors that store the time steps
and corresponding values of the continuous variables at which they were
recorded during the simulation. Other fields include algName and
sysName, which specify the names of the algorithm and system being
solved, alongside absQ, the absolute tolerance used during the
simulation. The stats field contains an instance of the Stats struct,
providing access to the performance metrics, while ft holds the final
time of the simulation. The Stats struct includes fields for totalSteps,
which counts the total number of steps taken throughout the simulation,
simulStepCount, representing the number of simultaneous updates during
the simulation, and evCount, which tallies the number of events that
occurred. Additionally, the numSteps vector keeps track of the number of
steps taken for each continuous variable involved in the simulation,
allowing for detailed analysis of the performance for individual
components.

The functions defined in conjunction with these structs enable various
operations related to the simulation's solutions. The createSol function
initializes a LightSol instance with the specified parameters, while the
getindex function provides a helper to access either the saved times or
variables easily. The evaluateSol function allows for the evaluation of
the solution at a specified time, supporting linear interpolation
between saved points. Furthermore, the solInterpolated functions enable
the generation of interpolated solutions across specified intervals,
creating new solution objects based on the interpolated values and
times. The show function is tailored to present the statistics in a
user-friendly format, offering insights into the performance of the
simulation. Overall, these constructs provide a robust framework for
managing and analyzing the outcomes of simulations involving ODEs.

**Plotting the Solution:** The plotSol function is designed to visually
represent the solution of a system of ODEs stored in the solution
object. It accepts a variety of optional parameters, allowing users to
customize their plots significantly. The function takes in a solution
object sol along with indices for the variables to be plotted (xvars).
Users can also specify a title, axis limits (xlims and ylims), and
visual styles, such as marker types and whether to display a legend. If
a title is not provided, it defaults to a composite of various
parameters, including the system name, algorithm name, absolute
tolerance, total steps taken, and more.Within the function, each
selected variable is plotted against its corresponding saved time,
utilizing different line styles for better distinction. If no specific
variables are selected, all variables in the solution are plotted by
default. After constructing the plot, the function checks for specified
axis limits and applies them accordingly. Finally, the function returns
the plot object, which can be further modified or saved. In addition to
basic plotting, the saveSol function allows users to save their plots as
PNG files, automatically generating filenames based on relevant
parameters and timestamps. This functionality is crucial for
documentation and sharing results in a reproducible manner. This
comprehensive plotting capability ensures that users can effectively
visualize the behavior of their systems, making the analysis of results
intuitive and informative.

**Finding the Error:** The getError function computes the relative error
of a solution compared to a specified reference function, enabling users
to assess the accuracy of their numerical results. The function iterates
over the saved time points, and it calculates the sum of the squares of
the differences between the numerical solution and the true values
provided by the reference function, as well as the sum of the squares of
the true values. The relative error is then derived as the square root
of the ratio of these sums, providing a quantitative measure of the
solution's deviation from the expected results. This comprehensive error
analysis is essential for validating the performance of numerical
methods and ensuring that the solutions produced are both reliable and
precise. By comparing the numerical results to known reference values,
users can gain insights into the accuracy and effectiveness of the
solver, enabling them to make informed decisions about their simulations
and analyses.