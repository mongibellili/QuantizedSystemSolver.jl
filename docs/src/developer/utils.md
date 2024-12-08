# Utils references

## Root finding

These methods find the minimum positive root of a linear and quadratic equations for order 1 and 2 respectively. For order 3, similar methods should be added to solve the cubic equation.

```@docs
QuantizedSystemSolver.minPosRoot(c::Float64,b::Float64, ::Val{1})
```

```@docs
QuantizedSystemSolver.minPosRoot(c::Float64,b::Float64,a::Float64,::Val{2})
```

```@docs
QuantizedSystemSolver.minPosRoot(coeff::Taylor0, ::Val{1})
```
```@docs
QuantizedSystemSolver.minPosRoot(coeff::Taylor0, ::Val{2})
```


## Scheduler


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


```@docs
QuantizedSystemSolver.updateScheduler(::Val{T}, nextStateTime::Vector{Float64}, nextEventTime::MVector{Z,Float64}, nextInputTime::Vector{Float64}) where {T, Z}
```