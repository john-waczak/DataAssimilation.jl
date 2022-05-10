# Kalman Filtering 
Given some model for the error covariance matrices ``Q_k`` and ``R_k``, we would like a method that propagates *both* our model **and** the errors forward. This way we may guarantee that the accuracy of our analysis doesn't come at the cost of higher uncertainty. 

The original implementation of the Kalman filter was for strictly linear systems. We will first develop the analysis for this simplified case adn then will generalize to the **Extended Kalman Filter** (EKF) that can handle fully nonlinear situations.

In the linear case, our system may be written as 
```math
\begin{aligned}
    u_{k+1}^{(t)} &= M_ku_k^{(t)} + \xi_{k+1}^{(b)} \\ 
    w_k &= H_ku_k^{(t)} + \xi_k^{(m)}
\end{aligned}
```
where ``M_k`` and ``H_k`` are now matrices defining the linear problem. 

The goal of the Kalman filter is to derive the analysis ``u^{(a)}`` which optimizes the trace of the analysis error covariance matrix (i.e. sum of squared errors): 
```math
\mathrm{Tr}\left( P_k^{(a)}\right) := \E[(u_k^{(t)}-u_k^{(a)})^T(u_k^{(t)}-u_k^{(a)})]
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
    &= M_ku_k^{(t)}+\xi_
\end{aligned}
```
