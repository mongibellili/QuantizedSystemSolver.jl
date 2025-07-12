# Introduction to Quantized System Methods


QSS methods, including LIQSS, are alternatives to traditional integration methods for solving ordinary differential equations (ODE). These methods are specially designed to solve challenges posed by the increasing complexity of modern systems. These methods offer advantages such as reduced computational cost and improved performance compared to traditional numerical integration methods in large sparse and stiff systems with frequent discontinuities. 

The main idea behind QSS methods is to divide the system state space into quantized regions and represent the system state in terms of these quantized values. The quantization of the states leads to an asynchronous discrete-event simulation model instead of a discrete time difference equation model. The QSS methods update the system status only when certain thresholds or conditions are met, instead of continuously simulating the system's behavior. All state derivatives remain constant if the states or inputs do not cross the next threshold, and the time of the crossing can be computed explicitly. QSS algorithms are asynchronous because The time, at which a state variable reaches its next threshold is different for separate states. For instance, state variables with large gradients will get to their thresholds more frequently than states with small gradients. If a state or input gets to its next threshold, a discrete event is triggered, and information is passed on to other integrators that depend on it [[1]](@ref refs1).

In solving the following system of $n$ ODEs in Eq.(1), classic methods look for the values of all variables after a small elapsed time, while QSS methods search for the time it takes to change the value of one variable by a small quantity. 
```math
\begin{equation}
   \dot X=f(X,t) , 
\end{equation}
```

where $X=[x_1,x_2...,x_n]^T$ is the state vector, $f:\mathbb{R}^n \times \mathbb{R}^+  \rightarrow \mathbb{R}^n$ is the derivative function, and $t$ is the independent variable.
In classic methods, the difference between $t_k$ (the current time) and $t_{k+1}$ (the next time) is called the step size. 
In QSS, besides the step size, the difference between $x_i(t_k)$ (the current value) and $x_i(t_{k+1})$ (the next value) is called the quantum $\Delta_i$.
Depending on the type of the QSS method (explicit or implicit), a new variable $q_i$ is set to equal $x_i(t_k)$  or $x_i(t_{k+1})$ respectively. 
$q_i$ is called the quantized state of $x_i$.
The use of the quantized state $q_{i}$ as the future value $x_{i}$ to calculate the derivatives is what makes the method implicit. 
Instead of solving the system shown in Eq.(1), QSS methods solve the following system in Eq.(2):
```math
\begin{align}
  & \dot X=f(Q,t) , 
\end{align}
```
where $Q=[q_1,q_2...,q_n]^T$ is the Quantized state vector.

The general form of a problem composed of a set of ODEs and a set of events that QSS is able to solve is described in the following: 

**System of $n$ ODEs:**

```math
\begin{align*}
  & \dot X=f(X,D,t) , 
\end{align*}
```

**System of $v$ events:**
```math
\begin{align*}
& if \; zc_v(x_i...,d_p...,t) \; i \in [1,n]  \;  ; \; p  \in [1,m] \\
& \qquad x_i=H(x_i...,d_p...,t) \\
& \qquad \qquad...\\
& \qquad d_p=L(x_i...,d_p...,t)  \\
& \qquad \qquad...\\
\end{align*}
```

where $n$ and $m$ are the number of state variables and discrete variables of the system respectively. $D=[d_1,d_2...,d_m]^T$ is the vector of the system discrete variables. $v$ is the number of events and $zc$ is an event condition, $H$ and $L$ are functions used in the effects of the event $zc$.

QSS methods were shown to have nice stability and error bound properties and they outperformed some classic solvers [[2]](@ref refs1).




## [References](@id refs1)

[1] Cellier, F. and Kofman, E. (2006). Continuous system. Springer, simulation.New York.

[2]  Kofman, E. (2009). Relative error control in quantization based integration. Latin American
Applied Research, 39(no.3):pp.231–238.



