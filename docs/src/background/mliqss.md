# Modified Linearly Implicit Quantized State System Methods

In LIQSS, all steps involves single variable updates depending on the future direction of $x_{i}$, as shown in bullet 4 of the **mLIQSS1 algorithm**.
While this intrinsic behavior of single cheap steps is advantageous, an issue arises from it. This asynchrony in updates can lead to unintended interactions between variables, causing them to stray from the true solution and trigger subsequent steps sooner than necessary. This behavior can result in cycles and unnecessary steps in the simulation. 
In essence, the lack of synchronization in variable updates can introduce unintended feedback loops, disrupting the simulation dynamics and potentially leading to inaccuracies or inefficiencies in the simulation results. Addressing this issue involves cycle detection mechanisms [[1]](@ref refs3).

The following conditions in Eq.(3) are used to detect cycles for order 1:



```math
    \begin{align}
     &  \dot x_j.dx_j<0  \;\; \; and \;\;\; \dot x_i.dx_i<0
    \end{align}
```



where $dx_i$ and $dx_j$ are linear approximations of the future derivatives. They are given as follows: 
```math
\begin{align}
 & dx_i=a_{ii}.q_i+a_{ij}.q_j+u_{ij}  \nonumber\\
 & dx_j=a_{jj}.q_j+a_{ji}.q_i+u_{ji}
\end{align}
```
 where $a_{ij}$ is the entry ($i$,$j$) of the Jacobian of the linearized derivative function and $u_{ij}$ is an affine coefficient.

 A modification to the original LIQSS algorithm was implemented in the following Algorithm:

**mLIQSS1 algorithm**
1. If a variable $i$ needs to change
    - Compute the elapsed time $e$ since the last update of variable $i$
    - Update its value using Taylor expansion: $x_{i}=x_{i}+\dot x_i.e$
    - Update the quantum $\Delta_i$
    - If $\dot x . (a_{ii}.(x_i+sign(\dot x).\Delta_i)+u_{ii}) >0$
        - Update the Quantized variable $q_i=x_i+sign(\dot x).\Delta_i$
    - Else
        - Update the Quantized variable $q_i=\frac{-u_{ii}}{a_{ii}}$   
    - Compute the next time when $q_i=x_i$ 
    - For any variable $j$ depends on $i$ such that $a_{ij}.a_{ji} \neq 0$
        - compute the prediction of derivatives of $i$ and $j$ using Eq.(2)
        - If the conditions in Eq.(1) are met:
            - perform a Backward Euler step as shown in Eq.(4) using **The Iteration approach in the simultaneous update order 1**
            - For any variable $k$ that depends on $j$
                -  update $x_k$, $\dot x_k$, $a_{kj}$, and next time of change.
    - For any variable $j$ depends on $i$
        - Update the variable $x_{j}=x_{j}+\dot x_j.e_j$
        - Update the derivative $\dot x_{j}=f_j(q,t)$ 
        - Compute the next time when $q_j=x_j \; or |q_j-x_j|=2\Delta_i$ 
    - For any zero crossing function $zc$ depends on $i$
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







**The Iteration approach in the simultaneous update order 1**

  - Update the step size $h=Final\;Simulation\;Time-current\;Sim\;Time$; 
  -   $q_i=x_i+ds_i(h)$ \; $q_j=x_j+ds_j(h)$
  - If ($|q_i-x_i|>\Delta_i$ or $|q_j-x_j|>\Delta_j$)
     - Update the step size $h_i=\frac{\Delta_i}{|\dot x_i|}$;  $h_j=\frac{\Delta_j}{|\dot x_j|}$; $h=min(h_i,h_j)$
     - Update the Quantized variable $q_i=x_i+ds_i(h)$ ; $q_j=x_j+ds_j(h)$
  - While ($|q_i-x_i|>\Delta_i$ or $|q_j-x_j|>\Delta_j$) and $AllowedIters>0$
     - Update $AllowedIters=AllowedIters-1$
     - Update the step size $h_i=h_i.\sqrt{\frac{\Delta_i}{|q_i-x_i|}}$; $h_j=h_j.\sqrt{\frac{\Delta_j}{|q_j-x_j|}}$; $h=min(h_i,h_j)$
     - Update the Quantized variable $q_i=x_i+ds_i(h)$ ; $q_j=x_j+ds_j(h)$
  - If $AllowedIters=0$
     - cancel the simultnaeous update and perform a LIQSS step


  
When a cycle is detected, a simultaneous update must occur. 
 In order to simultaneously update $q_i$ and $q_j$ as the future values of $x_i$ and $x_j$, the following equation is used:
```math
\begin{align}
 & q_i=x_i+h.dx_i   \nonumber\\
 & q_j=x_j+h.dx_j 
\end{align}
```
Using this equation and Eq.(2) is equivalent to performing a Backward Euler step, and it can formulated as follows:
```math
\begin{equation}
Q_{ij}=(I-h.A_{ij})^{-1}.(X_{ij}+h.U_{ij})
\end{equation}
```
where 
$Q_{ij}=$ $\begin{pmatrix}
  q_i \\ 
  q_j 
\end{pmatrix}$; $A_{ij}=$  $\begin{pmatrix}
  a_{ii} & a_{ij}\\ 
  a_{ji} & a_{jj}
\end{pmatrix}$; $X_{ij}=$ $\begin{pmatrix}
  x_i \\ 
  x_j 
\end{pmatrix}$; and $U_{ij}=$ $\begin{pmatrix}
  u_{ij} \\ 
  u_{ji} 
\end{pmatrix}$ 






The value of $h$ is computed analytically or via iterations as the maximum value that satisfies Eq.(5).
```math
\begin{align}
|q_i-x_i| \leq \Delta_i ;\ |q_j-x_j| \leq \Delta_j
\end{align}
```
Where $\Delta_i $ and $\Delta_j $ are the quantums of the variables $x_i$ and $x_j$ respectively.






## mLIQSS2
While higher-order methods offer improved accuracy in tracking the system's dynamics, they also introduce additional challenges. One of the key issues arises from the increased complexity in variable interactions, which can lead to discrepancies between the predicted and actual trajectories of the variables, potentially causing premature or unnecessary steps in the simulation.To address these challenges, conducting the cycle detection mechanism and updating the quatized states consider higher order terms.

Cycle detection mechanisms for Order 2 must consider both the first and second derivatives of the variables. The following conditions in Eq.(6) are used to detect cycles for Order 2:
```math
\begin{align}
  &  |\dot x_j-dx_j|>\frac{|\dot x_j+dx_j|}{2} \; or \; |\ddot x_j-ddx_j|>\frac{|\ddot x_j+ddx_j|}{2} \nonumber\\
  & \qquad \qquad \qquad  \qquad \qquad and \nonumber\\
  &  |\dot x_i-dx_i|>\frac{|\dot x_i+dx_i|}{2}  \; or \; |\ddot x_i-ddx_i|>\frac{|\ddot x_i+ddx_i|}{2}

  \end{align}
```
where $ddx_i$ and $ddx_j$ are linear approximations of the future second derivatives. They are given as shown in Eq.(7). 
```math
\begin{align}
 & ddx_i=a_{ii}.\dot q_i+a_{ij}.\dot q_j+\dot u_{ij}  \nonumber\\
 & ddx_j=a_{jj}.\dot q_j+a_{ji}.\dot q_i+\dot u_{ji}
\end{align}
```
 where $\dot u_{ij}$ is an affine coefficient.  

To account for these second-order effects, modifications to the original LIQSS algorithm were implemented, as shown in the **mLIQSS2 algorithm**. The algorithm ensures that the higher-order derivatives are properly integrated into the update rules.

The mLIQSS order 2 method is defined in the following Algorithm.
**mLIQSS2 algorithm** 
1. If a variable $i$ needs to change
    - Compute the elapsed time $e$ since the last update of variable $i$
    - Update the variable $x_{i}=x_{i}+\dot x_i.e_i+\frac{1}{2}\ddot x_i.e_i^2$
    - Update the quantum $\Delta_i$
    - Update the derivative $\dot x_{i}=\dot x_{i}+\ddot x_i.e$
    - Update the affine coefficient $u_{ii} = \dot x_i - a_{ii} . (q_{i}+\dot q_i.e_i)$
    - Update the derivative of the affine coefficient $\dot u_{ii} = \ddot x_i - a_{ii} . \dot q_i$
    - The Quantized variable $q_i(h)=\frac{(x_i-h.a_{ii}.x_i-h^2.(a.u_{ii}+\dot u_{ii})/2)}{(1 - h . a_{ii} + h^2.a_{ii}^2 / 2)}$
    - Update the step size $h=Final\;Simulation\;Time-current\;Sim\;Time$;  $q_i=q_i(h)$
    - If $|q_i-x_i|>\Delta_i$
        - Update the step size $h=\sqrt{\frac{2\Delta_i}{|ddx_i|}}$;  $q_i=q_i(h)$
    - While $|q_i-x_i|>\Delta_i$
        - Update the step size $h=h\sqrt{\frac{\Delta_i}{|q_i-x_i|}}$;  $q_i=q_i(h)$
    - Update the Quantized derivative $\dot q_i=\frac{a_{ii}.q_i+u_{ii}+h.\dot u_{ii}}{1-h.a_{ii}}$
    - Compute the next time when $q_i=x_i$ 
    - For any variable $j$ depends on $i$ such that $a_{ij}.a_{ji} \neq 0$
        -  compute the prediction of derivatives of $i$ and $j$ using Eq.(7)
        - If the conditions in Eq.(6) are met:
            -  perform a simultaneous update as shown in Eq.(9) using **The Iteration approach in the simultaneous update order 2**
            - For any variable $k$ that depends on $j$
                -  update $x_k$, $\dot x_k$, $\ddot x_k$, $a_{kj}$, and next time of change.
    - For any variable $j$ depends on $i$
        - Update the variable $x_{j}=x_{j}+\dot x_j.e_j+\frac{1}{2}\ddot x_j.e_j^2$
        - For any variable $k$ that $f_j$ depends upon
            - Update the Quantized variable $q_{k}=q_{k}+\dot q_k.e_k$ 
        - Update the derivatives $\dot x_{j}=f_j(q,t) \; and \ddot x_j=\dot f_j(Q,t)$ 
        -  Compute the next time when $q_j=x_j \; or |q_j-x_j|=2\Delta_i$ 
    - For any zero crossing function $zc$ depends on $i$
        - For any variable $k$ that $zc$ depends upon
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


   
 
  


**The Iteration approach in the simultaneous update order 2**
  - Update the step size $h=Final\;Simulation\;Time-current\;Sim\;Time$; 
  - Calculate $q_i$ and $q_j$ using Eq.(10)
  - If ($|q_i-x_i|>\Delta_i$ or $|q_j-x_j|>\Delta_j$)
     - Update the step size $h_i=\sqrt{\frac{2\Delta_i}{|ddx_i|}}$;  $h_j=\sqrt{\frac{2\Delta_j}{|ddx_j|}}$; $h=min(h_i,h_j)$
     - calculate $q_i$ and $q_j$ using Eq.(10)

 - While ($|q_i-x_i|>\Delta_i$ or $|q_j-x_j|>\Delta_j$) and $AllowedIters>0$
     - Update $AllowedIters=AllowedIters-1$
     - Update the step size $h_i=h_i.\sqrt{\frac{\Delta_i}{|q_i-x_i|}}$; $h_j=h_j.\sqrt{\frac{\Delta_j}{|q_j-x_j|}}$; $h=min(h_i,h_j)$
     - Calculate $q_i$ and $q_j$ using Eq.(10)
 - If $AllowedIters=0$
     -  Cancel the simultnaeous update and perform a LIQSS step
 -  Calculate $\dot q_i$ and $\dot q_j$ using Eq.(9)





When cycles are detected in Order 2, simultaneous updates must account for both the positions and second-order derivatives of the involved variables. This requires solving the following system of equations for $q_i$ and $q_j$:
```math
\begin{align}
  & q_i+h.\dot q_i=x_i+h.dx_i+h^2.ddx_i/2  \nonumber \\
  & \dot q_i=dx_i+h.ddx_i \nonumber\\
  & q_j+h.\dot q_j=x_j+h.dx_j+h^2.ddx_j/2 \nonumber\\
  & \dot q_j=dx_j+h.ddx_j 
 \end{align}
```
 Using Eq.(8), Eq.(2) and Eq.(7) is equivalent to the following:
```math
 \begin{align}
 \dot Q_{ij}=(I-h.A_{ij})^{-1}.(A_{ij}.Q_{ij}+U_{ij}+h.\dot U_{ij}) \nonumber\\
  Q_{ij}+h.\dot Q_{ij}=(X_{ij}+h.(A_{ij}.Q_{ij}+U_{ij})+h^2.(A_{ij}.\dot Q_{ij}+\dot U_{ij})/2)
 \end{align}
```
```math
 \begin{align}
   \alpha=X_{ij}+h.U_{ij}+h^2.\dot U_{ij}/2 \nonumber\\
   \beta=-(h.I-h^2.A_{ij}/2).(I-h.A_{ij})^{-1}.(U_{ij}+h.\dot U_{ij})+\alpha \nonumber\\
   \gamma=(I-h.A_{ij})+(h.I-h^2.A_{ij}/2).(I-h.A_{ij})^{-1}.A_{ij} \nonumber\\
   Q_{ij}= \gamma^{-1}.\beta
  \end{align}
```
  
 where 
 $\dot Q_{ij}=$ $\begin{pmatrix}
  \dot q_i \\ 
  \dot q_j 
 \end{pmatrix}$ and $\dot U_{ij}=$ $\begin{pmatrix}
  \dot u_{ij} \\ 
  \dot u_{ji} 
 \end{pmatrix}$ 




The value of $h$ is computed analytically or via iterations as the maximum value that satisfies Eq.(5).








## [References](@id refs3)


[1]  F. Pietro, G. Migoni, and E. Kofman. Improving linearly implicit quantized state system
methods. Simulation: Transactions of the Society for Modeling and Simulation International,
vol.95(no.2):pp.127–144, 2019.





