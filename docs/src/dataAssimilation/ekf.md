# Extended Kalman Filter
Given the nonlinear nature of many scientific models it is desirable to extend the **Kalman Filter** to be able to handle nonlinear models ``f(\cdot)`` (and by extension, their update function ``\mathcal{M}(\cdot)``), and nonlinear observation functions ``h(\cdot)``. This can be accomplished so long as these functions are sufficiently smooth (`C^1` to be precises) so as to admit valid Taylor approximations to first order. That is, 
```math
\begin{aligned}
    \mathcal{M}(u_{k}) &\approx \mathcal{M}(u_k^{(a)}) + D_{M}(u_k^{(a)})\xi_k^{(a)} & h(u_k) &\approx h(u_k^{(b)}) + D_h(u_k^{(b)})\xi_k^{(b)} \\ 
    D_{M} &:= \left[\dfrac{\partial \mathcal{M}_i}{\partial u_j} \right] & D_h &:= \left[ \dfrac{\partial h_i}{\partial u_j}\right]
\end{aligned}
```
where ``\mathcal{M}_i`` and ``h_i``  denote the ith component functions of ``\mathcal{M}`` and ``h``. 


Using these substitutions for the previously linear functions ``M_k`` and ``H_k``, we may follow the same derivation to obtain the following procedure. 

### 0. Initialization 
To begin we must choose values for ``u_0^{(a)}`` and ``P_0``. We must also provide models for ``Q_k`` and ``R_k``. 

### 1. Forecast Step 
```math
\begin{aligned}
    u_k^{(b)} &= \mathcal{M}(u_{k-1}^{(a)}) \\ 
    B_k &= D_M(u_{k-1}^{(a)})P_{k-1}D_M^T(u_{k-1}^{(a)}) + Q_k
\end{aligned}
```

### 2. Assimilation Step 
```math
\begin{aligned}
    K_k &= B_kD_M^T(u_k^{(b)})\left[ D_h(u_k^{(b)})B_kD_h^T(u_k^{(b)}) + R_k \right]^{-1}\\
    u_k^{(a)} &= u_k^{(b)} + K_k(w_k - h(u_k^{(b)})) \\ 
    P_k &= \left( I - K_kD_h(u_k^{(b)}) \right)B_k
\end{aligned}
```
