

# Tutorial

## Solving a Nonlinear ODE Problem with NLodeProblem in Julia
In this tutorial, we will go through the process of setting up, solving, querying, and plotting a nonlinear ordinary differential equation (ODE) problem using the NLodeProblem function. We will use a buck converter circuit model as an example.

## Step 1: Define the ODE Problem
The first step is to define the ODE problem using the NLodeProblem function. We provide the function with user code that specifies the parameters, variables, and differential equations.
Here is the example code for a buck converter circuit model with detailed explanations:

```julia
odeprob = NLodeProblem(quote
    name = (buck,) # Define the name of the problem for identification
    
    # Parameters
    C = 1e-4          # Capacitance in farads
    L = 1e-4          # Inductance in henrys
    R = 10.0          # Resistance in ohms
    U = 24.0          # Input voltage in volts
    T = 1e-4          # Switching period in seconds
    DC = 0.5          # Duty cycle
    ROn = 1e-5        # On-state resistance of the switch in ohms
    ROff = 1e5        # Off-state resistance of the switch in ohms
    
    # Discrete and continuous variables
    discrete = [1e5, 1e-5, 1e-4, 0.0, 0.0] # Initial discrete states
    u = [0.0, 0.0]                        # Initial continuous states
    
    # Rename for convenience
    rd = discrete[1]       # Diode resistance
    rs = discrete[2]       # Switch resistance
    nextT = discrete[3]    # Next switching time
    lastT = discrete[4]    # Last switching time
    diodeon = discrete[5]  # Diode state (on/off)
    il = u[1]              # Inductor current
    uc = u[2]              # Capacitor voltage
    
    # Helper equations
    id = (il * rs - U) / (rd + rs)  # Diode current calculation
    
    # Differential equations
    du[1] = (-id * rd - uc) / L     # Derivative of inductor current
    du[2] = (il - uc / R) / C       # Derivative of capacitor voltage
    
    # Events
    if t - nextT > 0.0 
        lastT = nextT             # Update last switching time
        nextT = nextT + T         # Schedule next switching time
        rs = ROn                  # Set switch to on-state
    end
    
    if t - lastT - DC * T > 0.0 
        rs = ROff                 # Set switch to off-state
    end
    
    if diodeon * id + (1.0 - diodeon) * (id * rd - 0.6) > 0
        rd = ROn                  # Diode is on
        diodeon = 1.0
    else
        rd = ROff                 # Diode is off
        diodeon = 0.0
    end
end)
```
### Explanation:
Parameters: These are constants used in the model.
C, L, R, U, T, DC, ROn, ROff are physical constants related to the buck converter circuit.\
Discrete and continuous variables: These represent the initial states of discrete and continuous variables in the model.\
Renaming variables: For convenience, the elements of discrete and u are renamed to more descriptive variable names (rd, rs, nextT, lastT, diodeon, il, uc).\
Helper equations: These are intermediate expressions needed for the differential equations.
id is the current through the diode.\
Differential equations: These represent the system's dynamics.
du[1] is the rate of change of the inductor current.
du[2] is the rate of change of the capacitor voltage.\
Events: These are conditions that modify the continous and discrete variables.
They handle the switching behavior of the buck converter and the diode's state transitions.
## Step 2: Solve the ODE Problem
Next, we solve the ODE problem using the solve function. We need to specify the problem, the algorithm, and the simulation settings.


```julia
# Define the time span for the simulation
tspan = (0.0, 0.001)  # Start at 0 seconds, end at 0.001 seconds
# Solve the ODE problem with the chosen algorithm and settings
sol = solve(odeprob, nmliqss2(), tspan, abstol=1e-4, reltol=1e-3)
```
### Explanation:
Time span: tspan defines the interval over which the solution is computed, from 0.0 to 0.001 seconds.\
Solver function: solve is used to compute the solution.\
odeprob is the ODE problem we defined.\
nmliqss2() specifies the algorithm used to solve the problem (e.g., qss2,nmliqss1... might be other algorithms).\
abstol and reltol are the absolute and relative tolerances for the solver, controlling the accuracy of the solution.
## Step 3: Query the Solution
After solving the problem, we can query the solution to extract useful information such as variable values at specific times, the number of steps, events, and more.

```julia
# Get the value of variable 2 at time 0.0005
value_at_time = sol(2, 0.0005)

# Get the total number of steps to end the simulation
total_steps = sol.totalSteps

# Get the number of simultaneous steps during the simulation
simul_step_count = sol.simulStepCount

# Get the total number of events during the simulation
event_count = sol.evCount

# Get the saved times and variables
saved_times = sol.savedTimes
saved_vars = sol.savedVars

# Print the results
println("Value of variable 2 at time 0.0005: ", value_at_time)
println("Total number of steps: ", total_steps)
println("Number of simultaneous steps: ", simul_step_count)
println("Total number of events: ", event_count)
```
### Explanation:
Value of variable 2 at a specific time: sol(2, 0.0005) returns the value of the second variable (capacitor voltage uc) at time 0.0005 seconds.
Total steps: sol.totalSteps gives the total number of steps taken by the solver to reach the end of the simulation.\
Simultaneous steps: sol.simulStepCount provides the number of steps where simultaneous updates occurred.\
Total events: sol.evCount gives the total number of events (e.g., switch state changes) during the simulation.\
Saved times and variables: sol.savedTimes and sol.savedVars store the time points and corresponding variable values computed during the simulation.
## Step 4: Plot the Solution
Finally, we can plot the solution to visualize the results. We have two options: generate a plot object for further customization or save the plot directly to a file.


``` julia
# Generate a plot object for all variables of the solution
plot_obj = plot_Sol(sol)
# Generate a plot object for variable 1 of the solution
plot_Sol(sol,1)
# Display the plot
display(plot_obj)


#plot  variables 1 and 2 of the solution
plot_Sol(sol,1,2,note=" ",xlims=(0.0,1.0),ylims=(-0.5,0.5),legend=false) 
#plot  the sum of variables 1 and 2 of the solution
plot_SolSum(sol,1,2)

# Save the plot to a file
save_Sol(sol, note=" ",xlims=(0.0,1.0),ylims=(-0.5,0.5),legend=false)

```
### Explanation:
Generate plot object: plot_Sol(sol) creates a plot object from the solution data.
Display the plot: display(plot_obj) shows the plot in the current environment.
Save the plot: save_Sol(sol, ...) saves the plot to a file *.png.



## User Documentation
More about the user documentation can be found in:

### [API](./Interface.md#application-programming-interface)

### [Available algorithms](./algorithm.md#available-algorithms)
