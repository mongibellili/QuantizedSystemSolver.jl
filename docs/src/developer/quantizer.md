# Quantizer 

## The Quantizer functions

The system handles different quantizer orders (Order 1 and Order 2...). It
defines methods for state integration, derivative computation, event
time computation, updating quantized values, and cycle detection
updates.


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

## Quantizer references

```@docs
QuantizedSystemSolver.integrateState(::Val{0}, x::Taylor0, elapsed::Float64)
```
```@docs
QuantizedSystemSolver.integrateState(::Val{1}, x::Taylor0, elapsed::Float64)
```
```@docs
QuantizedSystemSolver.integrateState(::Val{2}, x::Taylor0, elapsed::Float64)
```


```@docs
QuantizedSystemSolver.computeDerivative(::Val{1}, x::Taylor0, f::Taylor0)
```
```@docs
QuantizedSystemSolver.computeDerivative(::Val{2}, x::Taylor0, f::Taylor0)
```


```@docs
QuantizedSystemSolver.computeNextTime(::Val{1}, i::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0}, quantum::Vector{Float64})
```
```@docs
QuantizedSystemSolver.computeNextTime(::Val{2}, i::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0}, quantum::Vector{Float64})
```

```@docs
QuantizedSystemSolver.reComputeNextTime(::Val{1}, index::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0}, q::Vector{Taylor0}, quantum::Vector{Float64})
```

```@docs
QuantizedSystemSolver.reComputeNextTime(::Val{2}, index::Int, simt::Float64, nextTime::Vector{Float64}, x::Vector{Taylor0}, q::Vector{Taylor0}, quantum::Vector{Float64})
```

```@docs
QuantizedSystemSolver.computeNextInputTime(::Val{1}, i::Int, simt::Float64, elapsed::Float64, tt::Taylor0, nextInputTime::Vector{Float64}, x::Vector{Taylor0}, quantum::Vector{Float64})
```

```@docs
QuantizedSystemSolver.computeNextInputTime(::Val{2}, i::Int, simt::Float64, elapsed::Float64, tt::Taylor0, nextInputTime::Vector{Float64}, x::Vector{Taylor0}, quantum::Vector{Float64})
```


```@docs
QuantizedSystemSolver.computeNextEventTime(::Val{O},j::Int,ZCFun::Taylor0,oldsignValue::MMatrix{Z,2} ,simt::Float64,  nextEventTime :: MVector{Z,Float64}, quantum::Vector{Float64},absQ::Float64) where {O, Z}
```

```@docs
QuantizedSystemSolver.updateQ(::Val{1}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}, exactA::Function, d::Vector{Float64}, cacheA::MVector{1,Float64}, dxaux::Vector{MVector{1,Float64}}, qaux::Vector{MVector{1,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64})
```
```@docs
QuantizedSystemSolver.updateQInit(::Val{1}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}, exactA::Function, d::Vector{Float64}, cacheA::MVector{1,Float64}, dxaux::Vector{MVector{1,Float64}}, qaux::Vector{MVector{1,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64})
```
```@docs
QuantizedSystemSolver.Liqss_reComputeNextTime(::Val{1}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64})
```


```@docs
QuantizedSystemSolver.nmisCycle_and_simulUpdate(aij::Float64, aji::Float64, trackSimul, ::Val{1}, index::Int, j::Int, dirI::Float64, x::Vector{Taylor0}, q::Vector{Taylor0}, quantum::Vector{Float64}, exactA::Function, d::Vector{Float64}, cacheA::MVector{1,Float64}, dxaux::Vector{MVector{1,Float64}}, qaux::Vector{MVector{1,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64)
```



```@docs
QuantizedSystemSolver.updateQ(::Val{2}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}, exactA::Function, d::Vector{Float64}, cacheA::MVector{1,Float64}, dxaux::Vector{MVector{2,Float64}}, qaux::Vector{MVector{2,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64})
```
```@docs
QuantizedSystemSolver.updateQInit(::Val{2}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}, exactA::Function, d::Vector{Float64}, cacheA::MVector{1,Float64}, dxaux::Vector{MVector{2,Float64}}, qaux::Vector{MVector{2,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64})
```

```@docs
QuantizedSystemSolver.Liqss_reComputeNextTime(::Val{2}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64})
```

```@docs
QuantizedSystemSolver.nmisCycle_and_simulUpdate(aij::Float64, aji::Float64, trackSimul, ::Val{2}, index::Int, j::Int, dirI::Float64, x::Vector{Taylor0}, q::Vector{Taylor0}, quantum::Vector{Float64}, exactA::Function, d::Vector{Float64}, cacheA::MVector{1,Float64}, dxaux::Vector{MVector{2,Float64}}, qaux::Vector{MVector{2,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64)
``` 

## Index

```@index
Pages = ["quantizer.md"]
Order = [:type, :function]
```