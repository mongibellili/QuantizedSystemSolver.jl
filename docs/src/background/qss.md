
# Quantized State System Methods

## QSS1
The QSS1 method is the simplest of the QSS methods. It is based on the assumption that the state variables are piecewise constant between the thresholds. The state variables are updated only when a threshold is crossed. The time of the next threshold crossing is computed explicitly. The QSS1 method is a first-order method, and the state variables are updated using a taylor expansion, similar to the Euler method, given by the following equations:
```math
    \begin{align*}
      & x_{i}=x_{i}+\dot x_i.e_i  
    \end{align*}
```
where $x_{i}$ is the $i^{th}$ state variable, $\dot x_i$ is the derivative of this state variable, and $e$ is the elapsed time since its last update. The next time of change is computed as $$t_i^n=t_{current}+\frac{\Delta_i}{|\dot x_i|}$$ where $\Delta_i$ is the quantum of the variable $x_i$. The QSS1 method is defined in the foloowing Algorithm.

**QSS1 algorithm**

 1. If a variable $i$ needs to change
    - Compute the elapsed time $e_i$ since the last update of variable $i$
    - Update the variable $x_{i}=x_{i}+\dot x_i.e_i$
    - Update the quantum $\Delta_i$
    - Update the Quantized variable $q_i$ = $x_i$
    - Compute the next time of change $t_i^n$
    - For any variable $j$ depends on $i$
      - Update the variable $x_{j}=x_{j}+\dot x_j.e_j$
      - Update the derivative $\dot x_j=f_j(Q,t)$ 
      - Compute the next time when $q_j=x_j \; or |q_j-x_j|=2\Delta_j$ 
    - For any zero crossing function $zc$ depends on $i$
        - Update $zc$
        - Compute the next event time of $zc$ 
 2. If an event needs to occur
    - recheck validity of the event
    - execute the event and update the related quantized variables
    - For any variable $j$ depends on the event
        - Update the variable $x_{j}=x_{j}+\dot x_j.e_j$
        - Update the derivative $\dot x_j=f_j(Q,t)$ 
        - Compute the next time when $q_j=x_j \; or |q_j-x_j|=2\Delta_j$ 
    - For any zero crossing function $zc$ depends on the event
        - Update $zc$
        - Compute the next event time of $zc$ 


## QSS2
The QSS2 method is a second-order method. The state variables are updated using a Taylor expansion given by the following equations:
```math
  \begin{align*}
    & x_{i}=x_{i}+\dot x_i.e_i+\frac{1}{2}\ddot x_i.e_i^2  \\
    & \dot x_{i}=\dot x_{i}+\ddot x_i.e_i 
  \end{align*}
```
where $\ddot x_i$ is the second derivative of this state variable. The next time of change is computed as $$t_i^n=t_{current}+sqrt(\frac{\Delta_i}{|\ddot x_i|})$$. The QSS2 method is defined in the following Algorithm.

**QSS2 algorithm**


1. If a variable $i$ needs to change
    - Compute the elapsed time $e_i$ since the last update of variable $i$
    - Update the variable $x_{i}=x_{i}+\dot x_i.e_i+\frac{1}{2}\ddot x_i.e_i^2$
    - Update the quantum $\Delta_i$
    - Update the variable $\dot x_{i}=\dot x_{i}+\ddot x_i.e$
    - Update the Quantized variable $q_i$ = $x_i$
    - Update the Quantized derivative $\dot q_i$ = $\dot x_i$
    - Compute the next time of change $t_i^n$
    - For any variable $j$ depends on $i$
        - Update the variable $x_{j}=x_{j}+\dot x_j.e_j+\frac{1}{2}\ddot x_j.e_j^2$
        - For any variable $k$ that $f_j$ depends upon
            - Update the Quantized variable $q_{k}=q_{k}+\dot q_k.e_k$
        - Update the derivative $\dot x_j=f_j(Q,t)$ 
        - Update the second derivative $\ddot x_j=\dot f_j(Q,t)$ 
        - Compute the next time when $q_j=x_j \; or |q_j-x_j|=2\Delta_j$ 
    - For any zero crossing function $zc$ depends on $i$
        - For any variable $k$ that $zc$ depends upon
            - Update the Quantized variable $q_{k}=q_{k}+\dot q_k.e_k$
        - Update $zc$
        - Compute the next event time of $zc$ 
2. If an event needs to occur
    - recheck validity of the event
    - execute the event and update the related quantized variables
    - For any variable $j$ depends on the event
        - Update the variable $x_{j}=x_{j}+\dot x_j.e_j+\frac{1}{2}\ddot x_j.e_j^2$
        - For any variable $k$ that $f_j$ depends upon
            - Update the Quantized variable $q_{k}=q_{k}+\dot q_k.e_k$
        - Update the derivative $\dot x_j=f_j(Q,t)$ 
        - Update the second derivative $\ddot x_j=\dot f_j(Q,t)$ 
        - Compute the next time when $q_j=x_j \; or |q_j-x_j|=2\Delta_j$ 
    - For any zero crossing function $zc$ depends on the event
        - For any variable $k$ that $zc$ depends upon
            - Update the Quantized variable $q_{k}=q_{k}+\dot q_k.e_k$
        - Update $zc$
        - Compute the next event time of $zc$ 







