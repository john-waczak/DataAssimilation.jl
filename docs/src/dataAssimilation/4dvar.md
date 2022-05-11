# 4D-Var

The **3D-Var** algorithm attempts to optimize a cost function to obtain the ideal analysis *for each point where we have observation data*. This can become computationally expensive as we require model evaluations *and* an optimization routine for every observation point. An alternative approach is to simultaneously optimize accross all observations in order to obtain the ideal *initial condition* that acheive the best model fit. This approach is similar to sensativity analysis which seeks to fit a model's parameters to data. 

To begin, we construct the 4d-var cost function
```math
\begin{aligned}
    J(u_0) &= \frac{1}{2}\left( u_0 - u_0^{(b)} \right)^TB^{-1}\left( u_0 - u_0^{(b)} \right) + \frac{1}{2}\sum_k\left(w_k - h(u_k) \right)^TR_k^{-1}\left(w_k - h(u_k) \right) \\ 
           &= J_b(u_0) + J_m(u_0)
\end{aligned}
```
The first term is usefull if we already have an initial guess ``u_0^{(b)}`` for the inital condition in mind. If we do not have one, we may ommit this term. 


As before, we now want to optimize this cost function. To do so, we first observe that 
```math
    u_k = \mathcal{M}^{(k)}(u_0; \theta)
```

It is easy to obtain the gradient of ``J_0`` so we shall focus on the second term. We find that 
```math
\begin{aligned}
    \nabla_{u_0}J_m &= \nabla_{u_0}\Big\{ \sum_k \frac{1}{2}  \left(w_k - h(u_k) \right)^TR_k^{-1}\left(w_k - h(u_k) \right) \Big\}\\
                    &= - \sum_k \left[\dfrac{\partial }{\partial u_0}h\left(\mathcal{M}^{(k-1)}(u_0)\right) \right]^T R_k^{-1}\left(w_k - h(u_k) \right)\\
                    &= - \sum_k \left[D_h(u_k)D_M(u_{k-1})D_M(u_{k-2})\cdots D_M(u_0) \right]^T R_k^{-1}\left(w_k - h(u_k) \right)\\
                    &= - \sum_k \left[D_M^T(u_0)D_M^T(u_1)\cdots D_M^T(u_{k-1})D_h^T(u_k) \right] R_k^{-1}\left(w_k - h(u_k) \right)\\
\end{aligned}
```
