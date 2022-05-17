# Ensemble Kalman Filtering 

Data assimilation sees many applications including the use in systems governed by PDEs. In this case a discretization must be chosen to transform the PDEs into a set of coupled ODEs defined on each of the simulation grid points. As a consequence, performing data assimilation on such systems leads to state vectors with many millions or even billions of elements. This makes it computationally infeasable to determine the gargantuan covariance matrices and necesetates a revised approach. 

Previously, we had defined the background and analysis covariance matrices as 
```math
\begin{aligned}
    B &:= \E[\xi_b(\xi_b)^T] = \E\left[(u^{(t)-u^(b)})(u^({t})-u^{(b)})^T \right] \\
    P &:= \E[\xi_a(\xi_a)^T] = \E\left[(u^{(t)-u^(a)})(u^({t})-u^{(a)})^T \right] \\
\end{aligned}
```

In EnKF, we instead consider an ensemble of ``N`` copies of the system and use the following approximation to determine the ensemble covariances 
```math
\begin{aligned}
    B &\approx \frac{1}{N-1}\sum_{i=1}^{N}\left(u_i^{(b)} - \bar{u}^{(b)} \right)\left(u_i^{(b)} - \bar{u}^{(b)} \right)^T \\ 
    P &\approx \frac{1}{N-1}\sum_{i=1}^{N}\left(u_i^{(a)} - \bar{u}^{(a)} \right)\left(u_i^{(a)} - \bar{u}^{(a)} \right)^T \\ 
\end{aligned}
```
where the overbar denotes the mean, ``\bar{u}=\frac{1}{N}\sum_i u_i``. 

Naturally, the larger the ensemble, the better the approximation to the true state and covariances. 

##  The Algorithm 
1. Create an ensemble ``u_i^{(a)}(t_k)`` for ``i\in 1,...,N`` drawn from the distribution defined by the covariance matrix ``P_k``. 
2. Apply the forecast step to each ensemble member, i.e. 
```math
u_i^{(b)}(t_{k+1}) = \mathcal{M}(u_i^{(a)}(t_k); \theta) + \xi_p^{(i)}(t_{k+1})
```
3. Comptue the forecast covariance matrix 

```math
\begin{aligned}
    \bar{u}^{(b)}(t_{k+1}) &= \frac{1}{N}\sum_i^N u_i^{(b)}(t_{k+1}) \\ 
    \xi_i^{(b)}(t_{k+1}) &= u_i^{(b)}(t_{k+1}) - \bar{u}^{(b)}(t_{k+1}) \\ 
    B_{k+1} &= \frac{1}{N-1}\sum_i^N \xi_i^{(b)}(t_{k+1})\left(\xi_i^{(b)}(t_{k+1})\right)^T
\end{aligned}
```
4. Given an observation ``w(t_{k+1})``, generate an ensemble of *virtual observations* by sampling a Gaussian distribution with error covariance matrix ``R_{k+1}`` to obtain ``w_i(t_{k+1})``
5. Compute the Kalman gain matrix as before 
```math
K_{k+1} = B_{k+1}D_h(\bar{u}^{(b)}(t_{k+1}))^T\left[D_h(\bar{u}^{(b)}(t_{k+1}))B_{k+1}D_h(\bar{u}^{(b)}(t_{k+1}))^T + R_{k+1} \right]^{-1}
```
6. perform the analysis step 
```math
u_i^{(a)} = u_i^{(b)} + K_{k+1}\left(w_i(t_{k+1}) - h(u_i^{(b)}(t_{k+1})) \right)
```
7. Compute ensemble error covariance matrix for the analysis 
```math
\begin{aligned}
    u^{(a)}(t_{k+1}) &\approx \bar{u}^{(a)}(t_{k+1}) = \frac{1}{N}\sum_i^N u_i^{(a)}(t_{k+1}) \\ 
    \xi_i^{(a)} &= u_i^{(a)}(t_{k+1}) - \bar{u}^{(a)}(t_{k+1}) \\ 
    P_{k+1} &= \frac{1}{N-1}\sum_i^N \xi_i^{(a)}(t_{k+1})\left(\xi_i^{(a)}(t_{k+1}) \right)^T
\end{aligned}
```

We note that since each ensemble member may be propagated forward independently, it is easy to parallelize this algorithm to acheive good performance. 
