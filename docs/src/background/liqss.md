# Linearly Implicit Quantized State System Methods
Despite their advantages, explicit QSS methods are inefficient in solving stiff systems. In order to deal with stiffness, LiQSS evaluates the state derivative at future instants using a future quantized state value, as in classic implicit algorithms. The quantized state becomes a future value of the state and it is calculated smarter [[2]](@ref refs2). 

## LIQSS1

In first order, we either add or subtract the quantum $\Delta$ to the quantized state, and when when we sense a change of the derivative, we set the derivative equals to zero. 
```math
\begin{equation*}
 q_i=\begin{cases}
  x_i\mp\Delta_i & if  \dot x_i . \tilde{f}_i(x_i\mp\Delta_i) >0  \\
  \frac{-u_{ii}}{a_{ii}} & otherwise
\end{cases} 
\end{equation*}
```
where $\tilde{f}_i=a_{ii}.q_i+u_{ii}$  is a linear approximation of the derivative function $f_i$ at the future instant. 


The LIQSS order 1 method is defined in the following Algorithm:
**LIQSS1 algorithm**
1. If a variable $i$ needs to change
    - Compute the elapsed time $e$ since the last update of variable $i$
    - Update its value using Taylor expansion: $x_{i}=x_{i}+\dot x_i.e$
    - Update the quantum $\Delta_i$
    - If $\dot x . (a_{ii}.(x_i+sign(\dot x).\Delta_i)+u_{ii}) >0$
        - Update the Quantized variable $q_i=x_i+sign(\dot x).\Delta_i$
    - Else
        - Update the Quantized variable $q_i=\frac{-u_{ii}}{a_{ii}}$   
    - Compute the next time when $q_i=x_i$ 
    - For {any variable $j$ depends on $i$}
        - Update the variable $x_{j}=x_{j}+\dot x_j.e_j$
        - Update the derivative $\dot x_{j}=f_j(q,t)$ 
        - Compute the next time when $q_j=x_j \; or |q_j-x_j|=2\Delta_i$ 
    - For {any zero crossing function $zc$ depends on $i$}
        - Update $zc$
        - Compute the next event time of $zc$ 
2. If an event needs to occur
    - Recheck validity of the event
    - Execute the event and update the related quantized variables
    - For any variable $j$ depends on the event
        - Update the variable $x_{j}=x_{j}+\dot x_j.e_j$
        - Update the derivative $\dot x_j=f_j(Q,t)$ 
        - Compute the next time when $q_j=x_j \; or |q_j-x_j|=2\Delta_j$ 
    - For any zero crossing function $zc$ depends on the event
        - Update $zc$
        - Compute the next event time of $zc$ 



## LIQSS2
In Order 2, the updates involve not only the first-order variables but also their second-order derivatives. Similar to Order 1, these updates depend on the future behavior of the state variables, but with a focus on both their positions and rates of change as reflected in the second-order derivative terms.
To define LIQSS2, we need to consider the second derivative of the state variables. The linear approximation of the future second derivative is given as follows:
$\ddot x_i=a_{ii}.\dot q_i+\dot u_{ii}$.

The LIQSS order 2 method is defined in the following Algorithm.
**LIQSS2 algorithm** 
1. If a variable $i$ needs to change
    - Compute the elapsed time $e$ since the last update of variable $i$
    - Update the variable $x_{i}=x_{i}+\dot x_i.e_i+\frac{1}{2}\ddot x_i.e_i^2$
    - Update the quantum $\Delta_i$
    - Update the derivative $\dot x_{i}=\dot x_{i}+\ddot x_i.e$
    - Update the affine coefficient $u_{ii} = \dot x_i - a_{ii} . (q_{i}+\dot q_i.e_i)$
    - Update the derivative of the affine coefficient $\dot u_{ii} = \ddot x_i - a_{ii} . \dot q_i$
    - The Quantized variable $q_i(h)=\frac{(x_i-h.a_{ii}.x_i-h^2.(a.u_{ii}+\dot u_{ii})/2)}{(1 - h . a_{ii} + h^2.a_{ii}^2 / 2)}$
    - Update the step size $h=Final\;Simulation\;Time-current\;Sim\;Time$;  $q_i=q_i(h)$
    - If {$|q_i-x_i|>\Delta_i$}
        - Update the Quantized variable $h=\sqrt{\frac{2\Delta_i}{|ddx_i|}}$;  $q_i=q_i(h)$
    - While $|q_i-x_i|>\Delta_i$
        - Update the Quantized variable $h=h\sqrt{\frac{\Delta_i}{|q_i-x_i|}}$;  $q_i=q_i(h)$
    - Update the Quantized derivative $\dot q_i=\frac{a_{ii}.q_i+u_{ii}+h.\dot u_{ii}}{1-h.a_{ii}}$
    - Compute the next time when $q_i=x_i$ 
    - For any variable $j$ depends on $i$
        - Update the variable $x_{j}=x_{j}+\dot x_j.e_j+\frac{1}{2}\ddot x_j.e_j^2$
        - For {any variable $k$ that $f_j$ depends upon}
            - Update the Quantized variable $q_{k}=q_{k}+\dot q_k.e_k$ 
        - Update the derivatives $\dot x_{j}=f_j(q,t) \; and \ddot x_j=\dot f_j(Q,t)$ 
        - Compute the next time when $q_j=x_j \; or |q_j-x_j|=2\Delta_i$ 
    - For any zero crossing function $zc$ depends on $i$
        - For {any variable $k$ that $zc$ depends upon}
            - Update the Quantized variable $q_{k}=q_{k}+\dot q_k.e_k$
        - Update $zc$
        - Compute the next event time of $zc$ 
2. If an event needs to occur
    - Recheck validity of the event
    - Execute the event and update the related quantized variables
    - For any variable $j$ depends on the event
        - Update the variable $x_{j}=x_{j}+\dot x_j.e_j+\frac{1}{2}\ddot x_j.e_j^2$
        - For any variable $k$ that $f_j$ depends upon
            - Update the Quantized variable $q_{k}=q_{k}+\dot q_k.e_k$
        - Update the derivatives $\dot x_j=f_j(Q,t) \; and \ddot x_j=\dot f_j(Q,t)$ 
        - Compute the next time when $q_j=x_j \; or |q_j-x_j|=2\Delta_j$ 
    - For any zero crossing function $zc$ depends on the event
        - For any variable $k$ that $zc$ depends upon
            - Update the Quantized variable $q_{k}=q_{k}+\dot q_k.e_k$
        - Update $zc$
        - Compute the next event time of $zc$ 





## [References](@id refs2)


[1]  Migoni, G., Kofman, E., and Cellier, F. (2012). Quantization-based new integration methods
for stiff odes. Simulation: Transactions of the Society for Modeling and Simulation International,
88(4):378–407.

[2]  G. Migoni, M. Bortolotto, E. Kofman, and F. Cellier. Linearly implicit quantization-based
integration methods for stiff ordinary differential equations. Simulation Modelling Practice
and Theory, vol.35:pp.118–136, 2013.










