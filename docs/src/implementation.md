# Implementation
### Package Structure {#ch5:section:Package-description}

While the package is optimized to be fast, extensibility is not
compromised. It is divided into 3 entities that can be extended
separately: Problem, Algorithm, and Solution. The rest of the code is to
create these entities and glue them together as shown in the Figure.

::: center
![The QSS Solver Structure.](./assets/img/diagram.png){#fig:qssdiag width="10cm"
height="6cm"}
:::

The API was designed to provide an easier way to handle events than the
approach provided by the classic integration solvers. Inside an
NLodeProblem function, the user may introduce any parameters, variables,
equations, and events:

        odeprob = NLodeProblem(quote 
        name
        parameters
        discrete and continous variables
        helper expressions
        differential equations
        if-statments for events)

The output of this function is an object of type Problem, and it is
passed to the solve function along any other configuration arguments
such as the algorithm type, the time span and the tolerance. The solve
function dispatches on the given algorithm and start the numerical
integration.

        tspan = (start time, end time)
        sol= solve(odeprob,algorithm,tspan,abstol=...,reltol=...)    

A the end, a solution object is produced that can be queried and
plotted.

        sol(time,idxs=variable index) 
        sol.stats
        plot(sol)    


### Internal Code Explained {#ch5:section:code-explained}

The solver structure is designed for modularity, allowing for different
QSS algorithms and quantizer orders to be used in the simulation. The
core components of the solver include the NLodeProblem function, the
solve function, the scheduler, the quantizer, and the solution object.
These components work together to parse user-provided code, solve the
system of ODEs using QSS algorithms, and store the results for analysis.
The following sections provide an in-depth explanation of the key
components and their interactions within the solver framework.

#### creating the problem

The NLodeProblem function is the entry point for defining a new problem
to be solved by the QSS solver. It takes user-provided code, which
includes system parameters, variables, equations, and event logic, and
constructs a Problem object that encapsulates all the necessary
information for the solver to simulate the system such as problem
dimensions, dependencies, and equations. The function works by parsing
the user code and extracting relevant data to populate the Problem
object.

*NLODEDiscProblem{PRTYPE,T,Z,Y,CS}:* This is the struct that holds all
the necessary data for a nonlinear ordinary differential equation (ODE)
problem with discrete events. The structure includes various fields such
as initial conditions, discrete variables, Jacobians, event
dependencies, and other data related to how the problem is formulated.
This structure serves as the core data holder for the problem and will
be used in the solver. It is a parametric abstract type that has the
following parameters:

PRTYPE: The type of the problem (to distinguish between various types,
and allow future extension of the solver to handle new types).

T: The number of continuous variables (state variables).

Z: The number of zero-crossing functions, which are used to detect
events.

Y: The actual number of events.

CS: Cache size, which is used to store intermediate operations.

The use of abstract types in this context allows for flexibility and
extensibility in the solver. By defining these abstract types, the code
can be easily adapted to handle different types of problems, algorithms,
and solutions without needing to modify the core solver logic. This
design choice enhances the maintainability and scalability of the
solver, making it easier to add new features or support additional
problem types in the future.

**NLodeProblemFunc:** After an initial preparation performed by the The
NLodeProblem function, The function NLodeProblemFunc takes the resulting
expressions to continue constructing an instance of the NLODEDiscProblem
structure. It works in several key stages:

*Initialization:* The function begins by initializing vectors and
dictionaries that will hold equations (equs), Jacobian dependencies
(jac), zero-crossing functions (ZCjac), and event dependencies. These
serve to store the different types of equations and their relationships.

*Processing ODEs:* It loops through each of the ODE expressions provided
by the user. Depending on the type of expression (discrete variables,
differential equations, or loop constructs), it processes the right-hand
side (RHS) of the equation. For differential equations, it extracts
dependencies to build the Jacobian and transform the equations into a
more appropriate form for further use. Special cases are handled, such
as if the RHS is a number or a symbol.

*Handling Events:* The function also processes event-related constructs
(if conditions) that correspond to different points where the system
might undergo discrete changes. It process the RHS of the event
equations, transforms them into a suitable form, and builds the
necessary dependency structures. Specifically, it constructs how
discrete and continuous variables influence one another through the
events.

*Constructing the Function Code:* After processing all ODEs and events,
the function dynamically generates a Julia function code needed to store
the system of ODEs and events. This code is built into a function that
handles different cases (i.e., which equation to evaluate based on an
index of a state change or an event).

**Building Dependencies:** Several helper functions that build the
dependencies between variables, events. They build dependency vectors
that track how discrete and continous variables influence the system.
This is used to know what variables to update and determine when
specific events should be checked. By tracking the relationships between
variables and events, the solver can determine the appropriate actions
to take at each time step. The dependencies are stored in the following
vectors:

\-$jac$: It determines which variables affect a derivative.

\-$ZCjac$: It determines which variables affect a zero-crossing
function.

\-$SD$: It determines which derivatives that are affected by a given
variable.

\-$SZ$: It determines which zero-crossing functions that are affected by
a given variable.

\-$HZ$: It tells which Zero-crossing functions influenced by a given
event.

\-$HD$: It tells which derivatives influenced by a given event.

Here's a quick summary and what each helper function is doing:

*extractJacDepNormal:* It Extracts the dependencies for normal
(non-loop) expressions. It updates the Jacobian matrix $jac$ and a
dictionary $dD$ for tracking dependencies of derivatives to discrete
variables.

*extractJacDepLoop:* Similar to extractJacDepNormal, but specifically
for loop expressions. It tracks dependencies across loop iterations.

*extractZCJacDepNormal:* It Extracts zero-crossing Jacobian dependencies
for discrete variables ($dZ$), and it updates $zcjac$, $SZ$.

*createDependencyToEventsDiscr:* It maps discrete dependencies (dD, dZ)
to specific events, it and constructs dependency matrices HZ and HD from
the discrete variables only.

*createDependencyToEventsCont:* Similar to
createDependencyToEventsDiscr, but for continuous dependencies (SD, sZ),
and it updates the matrices HZ and HD from the continuous variables
only.

*unionDependency:* Merges the two previous sets of dependencies
(continuous and discrete) into the final matrices HZ and HD.

#### The solve function

The solve function is the primary interface for solving non-linear ODE
problems using various QSS (Quantized State Systems) algorithms. It does
this by dispatching on the problem and algorithm types to select the
right solver.

**QSS Algorithm:** The diagram shows three specific algorithms: (QSS
Algorithm, LiQSS Algorithm, mLiQSS Algorithm). The Algorithm type is
parametric on the solver name (N) and order (O). When a user doesn't
provide a solver explicitly, the solve function defaults to using the
modified second-order implicit algorith (mLiQSS2). This approach allows
flexibility: the user can select a different algorithm by providing an
alternative QSS Algorithm without needing to modify the solve function.
Based on the QSS Algorithm provided, it either selects a basic QSS
integration method or, in the case of LiQSS, constructs additional data
structures needed for implicit integration. The method defaults to
sparse handling as false (Val(false)), tolerances (abstol=1e-4 and
reltol=1e-3), and a maximum number of iterations (maxiters=1e7). These
parameters can be adjusted based on the problem's complexity and desired
accuracy.

**Helper Functions:**

*createCommonData:* This function sets up the common data required by
the QSS solver. It initializes all necessary vectors (x, q, tx, tq,
nextStateTime, etc.) used in the integration process, and it
pre-allocates a cache of Taylor series to avoid repeated memory
allocation during integration (taylorOpsCache).

*getClosure:* This function creates a closure for the Jacobian and the
dependency matrices, allowing in a flexible in a way that enables easy
extension. The closure is used to access the Jacobian matrix during the
integration process.

*createLiqssData:* For LiQSS solvers, this function sets up additional
data structures and auxiliary variables used for storing linear
approximation coefficients.

#### The scheduler function

The scheduler component retrieves the next index (a state change or an
event) for processing. It determines the next action by comparing three
types of timing events: state updates, events, and input updates. It
initializes minimum time variables for each action type to infinity and
iterates through the provided vectors to find the minimum times and
their corresponding indices. Based on these comparisons, it decides
whether the next scheduled action is an input update, an event, or a
state update. If no valid actions are found (indicated by a zero index),
the function assigns a default state indicating the next action is a
state update at time infinity. Finally, it returns a tuple containing
the index of the next action, its time, and a symbol representing the
type of action to take next (:STINPUT, :STEVENT, or :STSTATE). This
structure ensures that the simulation progresses accurately and
efficiently.

#### The Quantizer functions

The system handles different quantizer orders (Order 1 and Order 2). It
defines methods for state integration, derivative computation, event
time computation, updating quantized values, and cycle detection
updates.

**Utils:** This section contains utility functions like root finding
that assist in the solution process.

**computeNextTime:** for first-order, the state of the system changes at
a constant rate. The core calculation takes place when the first
derivative of the state, represented by $x[i][1]$, is non-zero. In this
case, the function determines the time to the next event by dividing a
quantum threshold by this derivative. Additionally, to prevent numerical
issues, it ensures that this calculated time-step does not fall below a
predefined minimum value, $absDeltaT$. If the first derivative is
extremely small or essentially zero, the function adjusts it to avoid
potential numerical instabilities that could arise from very small time
increments.

For second-order, the system state evolves with both a rate of change
(first derivative) and acceleration (second derivative). The core
calculation occurs when the second derivative of the state, represented
by $x[i][2]$, is non-zero. In this case, the function computes the time
to the next event by using the square root of the ratio between a
quantum threshold and the second derivative. This ensures that the
time-step reflects the influence of the system's acceleration.
Additionally, to prevent numerical issues (such as division by zero or
overly small time-steps), a minimum delta time (absDeltaT) is enforced.
If the computed time-step is smaller than this threshold, the function
adjusts the second derivative to maintain stability. If the second
derivative is zero but the first derivative is non-zero, the time to the
next event is calculated based on the first derivative. The function
ensures that the time-step does not drop below absDeltaT, adjusting the
first derivative if necessary. If both derivatives are zero, the system
is assumed to have no change, and the next event time is set to infinity
(Inf).

**recomputeNextTime:** The reComputeNextTime functions are used in the
explicit QSS algorithms to enable the recalculation of the next time
after interactions between different variables, such as variable i and
variable j. These functions determine the time until the system crosses
the quantum threshold by solving polynomial equations derived from the
difference between the quantized state and the actual state. In
situations where the quantum threshold has already been surpassed, they
promptly return a very small time increment (e.g., simt + 1e-12) to
trigger an immediate update of the system's state.

**LiqssrecomputeNextTime:** In the LIQSS methods of first-order, the
recomputeNextTime function calculates the next event time based on the
current state (x), its first derivative (x1), and the quantized state
(q). First, if the difference between the current state and the
quantized state exceeds twice the quantum size, the next event is
scheduled almost immediately with a small time-step. Otherwise, if the
derivative is non-zero, the function computes the time-step by dividing
the state difference by the derivative. If the result is positive, this
time is added to the current simulation time (simt). If negative, it
adjusts the time-step by either adding or subtracting twice the quantum
size based on the direction of change. If the derivative is zero,
indicating no change, the event time is set to infinity. Lastly, if the
computed time is in the past, the function resets it to a far future
time to prevent any premature events. This ensures that the system
evolves smoothly and accurately.

In the second-order LIQSS method, the recomputeNextTime function
calculates the next event time by considering both the state (x) and its
first (x1) and second (x2) derivatives, along with the quantized state
(q). The function constructs a polynomial using the current state,
derivative values, and the second derivative, then calculates the next
event time by finding the smallest positive root of this polynomial.

**computeNextInputTime:** The computeNextInputTime functions focus on
computing the next action when derivatives depend only on time. They
assess changes in the derivatives over a specified elapsed time to
compute the time increment until the next input action. If the
derivatives are null, the function reverts to handling the system in a
lower-order thus simplifying the calculation.

**computeNextEventTime:** The computeNextEventTime function calculates
the next event time based on the zero-crossing function ($ZCFun$) of a
system. It first checks if a sign change has occurred in the
zero-crossing function, indicating that the system is leaving zero and
should be considered an event, provided the previous value is
significantly different from zero (to prevent duplicate events). If a
sign change is detected, the function updates the event time to the
current simulation time. If both the old and new values of the
zero-crossing function are zero, it sets the next event time to
infinity, indicating no event should occur. For cases where the old and
new values have the same sign, it calculates the minimum positive root
of the zero-crossing function, representing the time of the next event,
ensuring that this time is not too close to zero to avoid spurious
events. The function then updates the old sign values for future
comparisons.

#### The solution

This is where the results of the simulation are stored, including:
Settings, Time points (times\[\]), Variables (vars\[\]), Statistics
(stats).

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