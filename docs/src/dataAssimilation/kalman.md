# Kalman Filtering 
Given some model for the error covariance matrices ``Q_k`` and ``R_k``, we would like a method that propagates *both* our model **and** the errors forward. This way we may guarantee that the accuracy of our analysis doesn't come at the cost of higher uncertainty. 

The original implementation of the Kalman filter was for strictly linear systems. We will first develop the analysis for this simplified case adn then will generalize to the **Extended Kalman Filter** (EKF) that can handle fully nonlinear situations.

In the linear case, our system may be written as 
```math
\begin{aligned}
    u_{k+1}^{(t)} &= M_ku_k^{(t)} + \xi_{k+1}^{(p)} \\ 
    w_k &= H_ku_k^{(t)} + \xi_k^{(m)}
\end{aligned}
```
where ``M_k`` and ``H_k`` are now matrices defining the linear problem. 

The goal of the Kalman filter is to derive the analysis ``u^{(a)}`` which optimizes the trace of the analysis error covariance matrix (i.e. sum of squared errors): 
```math
\mathrm{Tr}\left( P_k\right) := \E[(u_k^{(t)}-u_k^{(a)})^T(u_k^{(t)}-u_k^{(a)})]
```
Finding the analysis consists of two steps: the forecast step and the assimilation step.


## Forecast Step
Assume we have the analysis at time ``t_k`` denoted ``u_k^{(a)}``. Then the forecast for time ``t_{k+1}`` is
```math
    u_{k+1}^{(b)} = M_ku_k^{(a)}
```
The background error is therefore 
```math
\begin{aligned}
    \xi_{k+1}^{(b)} &= u_{k+1}^{(t)} - u_{k+1}^{(b)} \\ 
    &= M_ku_k^{(t)}+\xi_{k+1}^{(p)} - M_{k}u_k^{(a)} \\ 
    &= M_k\left(u_k^{(t)}-u_k^{(a)} \right) + \xi_{k+1}^{(p)} \\ 
    &= M_k\xi_k^{(a)} + \xi_{k+1}^{(p)}
\end{aligned}
```
We may now evaluate the covariance matrix of our background estimate as: 
```math
\begin{aligned}
    B_{k+1} &= \E[\xi_{k+1}^{(b)}(\xi_{k+1}^{(b)})^T] \\ 
    &= \E\left[\left(M_k\xi_k^{(a)} + \xi_{k+1}^p \right) \left(M_k\xi_k^{(a)} + \xi_{k+1}^p \right)^T \right] \\ 
\end{aligned}
```
If we presume that ``\E[\xi_k^{(b)}(\xi_{k+1}^{(p)})^T] = 0``, then the cross terms vanish and we are left with 
```math
\boxed{B_{k+1} = M_kP_kM_k^T + Q_{k+1}}
```
