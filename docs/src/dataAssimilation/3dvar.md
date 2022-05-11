# 3D-Var 
For the **Kalman Filter** and the **EKF**, we derived the optimal way to combine observation with simulation so as to minimize the trace of the analysis error covariance matrix, ``P_k``. An alternative approach is to recast the problem as a pure optimzation problem where rather than finding a filter ``K_k`` that will add an innovation to ``u_k^{(b)}`` to obtain the analysis ``u_k^{(a)}``, we obtain the analysis by optimizing the following cost function
```math
J(u) = \frac{1}{2}\left(w - h(u) \right)^TR^{-1}\left(w - h(u) \right) + \frac{1}{2}\left(u - u^{(b)} \right)^TB^{-1}\frac{1}{2}\left(u - u^{(b)} \right)
```
which we can justify as coming from the joint probability distribution assuming Gaussian errors
```math
\mathcal{P}(u|w) = C\exp\left(- \frac{1}{2}\left(u - u^{(b)} \right)^TB^{-1}\frac{1}{2}\left(u - u^{(b)} \right) \right)\cdot\exp\left(-  \frac{1}{2}\left(w - h(u) \right)^TR^{-1}\left(w - h(u) \right) \right)
```
with model error covariance ``B`` and measurement error covariance ``R`` as before. This is clearly a *very strong assumption*.

To optimize ``J(u)``, we begin by taking it's gradient.
```math
\nabla_uJ(u) &= -D_h^TR^{-1}(w-h(u)) + B^{-1}(u-u^{(b)})
```
Thus, finding the analysis ``u^{(a)}`` ammounts to solving the system 
```math
a-D_h^TR^{-1}(w-h(u^{(a)})) + B^{-1}(u^{(a)}-u^{(b)}) = 0
```

As for Kalman filtering, let's begin with the assumption that our model and observation function are linear. 

## Linear Case
Suppose that we have ``h(u) = Hu`` so that ``D_h(u) = H``. Then, we have 
```math 
\begin{aligned}
    D_h^TR^{-1}(w-Hu^{(a)}) &= B^{-1}(u^{(a)}-u^{(b)}) \\ 
    D_h^TR^{-1}w - D_h^TR^{-1}Hu^{(a)} &= B^{-1}u^{(a)} - B^{-1}u^{(b)} \\ 
    \left(D_h^TR^{-1}H + B^{-1} \right)u^{(a)} &= D_h^TR^{-1} + B^{-1}u^{(b)} \\
    \left(H^TR^{-1}H + B^{-1} \right)u^{(a)} &= H^TR^{-1} + B^{-1}u^{(b)}
\end{aligned}
```
