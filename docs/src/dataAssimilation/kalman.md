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

If we presume that $\E[\xi_k^{(b)}(\xi_{k+1}^{(p)})^T] = 0$, then the cross terms vanish and we are left with 

```math
\boxed{B_{k+1} = M_kP_kM_k^T + Q_{k+1}}
```

Thus we now have the background (i.e forecast) estimate of the state at $t_{k+1}$ and its covariance matrix. Given a measurement $w_{k+1}$ at the same time with covariance matrix $R_{k+1}$, then we may now perform the assimilation step where we fuse the two sources of information to obtain $u_{k+1}^{(a)}$ and $P_{k+1}$.

## Data Assimilation Step
Let's suppose that the analysis has the form 

```math
u_{k+1}^{(a)} = \nu + K_{k+1}w_{k+1}
```

for some vector $\nu\in\R^n$ and matrix $K_{k+1}\in\R^{m\times n}$. In a perfect world, we would have $\E[u_{k}^{(t)}-u_{k}^{(a)}] = 0$. Therefore, 

```math
\begin{aligned}
    0 &= \E[u_k^{(t)} - u_k^{(a)}] \\ 
    &= \E[(u_k^{(b)} + \xi_k^{(b)}) - (\nu + K_kw_k)] \\ 
    &= \E[(u_k^{(b)} + \xi_k^{(b)}) - (\nu + K_kH_ku_k^{(t)} + K_k\xi_k^{(m)})] \\ 
    &= \E[u_k^{(b)}] + \E[\xi_k^{(b)}] - \E[\nu] -K_kH_k\E[u_k^{(t)}] - K_k\E[\xi_k^{(m)}]\\
    &= u_k^{(b)} + 0 - \nu - K_kH_ku_k^{(b)} - 0 \\ 
    &= u_k^{(b)} - \nu - K_kH_ku_k^{(b)} \\ 
    \Rightarrow \nu &= u_k^{(b)} - K_kH_ku_k^{(b)}
\end{aligned}
```

which we now substitute to obtain 

```math
\boxed{u_k^{(a)} = u_k^{(b)} + K_k(w_k - H_ku_k^{(b)})}
```

Now that we know the form for the analysis we may derive the optimal matrix $K_k$ by optimization of $P_k$. We have

```math
\begin{aligned}
	\xi_k^{(a)} &= u_k^{(t)} - u_k^{(a)} \\ 
                &= M_{k-1}u_{k-1}^{(t)} + \xi_{k}^{(p)} - u_k^{(b)} - K_k\left(w_k - H_ku_k^{(b)} \right) \\ 
                &= M_{k-1}u_{k-1}^{(t)} + \xi_{k}^{(p)} - M_{k-1}u_{k-1}^{(a)} - K_k\left(H_ku_k^{(t)} + \xi_k^{(m)} - H_ku_k^{(b)} \right) \\ 
                &= M_{k-1}u_{k-1}^{(t)} + \xi_{k}^{(p)} - M_{k-1}u_{k-1}^{(a)} - K_kH_ku_k^{(t)} - K_k\xi_k^{(m)} + K_kH_ku_k^{(b)} \\ 
                &= M_{k-1}u_{k-1}^{(t)} + \xi_{k}^{(p)} - M_{k-1}u_{k-1}^{(a)} - K_kH_ku_k^{(t)} - K_k\xi_k^{(m)} + K_kH_ku_k^{(b)} \\ 
                &= \Big\{ M_{k-1}(\xi_{k-1}^{(a)}+u_{k-1}^{(a)}) + \xi_{k}^{(p)} - M_{k-1}u_{k-1}^{(a)} \Big\} - K_kH_ku_k^{(t)} - K_k\xi_k^{(m)} + K_kH_ku_k^{(b)} \\ 
                &= \Big\{ M_{k-1}\xi_{k-1}^{(a)} + \xi_{k}^{(p)} \Big\} - K_kH_ku_k^{(t)} + K_kH_ku_k^{(b)} - K_k\xi_k^{(m)}\\ 
                &= M_{k-1}\xi_{k-1}^{(a)} + \xi_{k}^{(p)} - K_kH_k(M_{k-1}u_{k-1}^{(t)} + \xi_k^{(b)}) + K_kH_ku_k^{(b)} - K_k\xi_k^{(m)}\\ 
                &= M_{k-1}\xi_{k-1}^{(a)} + \xi_{k}^{(p)} - K_kH_kM_{k-1}(\xi_{k-1}^{(a)} + u_{k-1}^a) - K_kH_k\xi_k^{(b)} + K_kH_ku_k^{(b)} - K_k\xi_k^{(m)}\\ 
                &= M_{k-1}\xi_{k-1}^{(a)} + \xi_{k}^{(p)} - K_kH_kM_{k-1}(\xi_{k-1}^{(a)} + u_{k-1}^a) - K_kH_k\xi_k^{(b)} + K_kH_kM_{k-1}u_{k-1}^{(a)} - K_k\xi_k^{(m)}\\ 
                &= M_{k-1}\xi_{k-1}^{(a)} + \xi_{k}^{(p)} - K_kH_kM_{k-1}\xi_{k-1}^{(a)} - K_kH_k\xi_k^{(b)} - K_k\xi_k^{(m)}\\ 
                &= \big(I-K_kH_k \big)(M_{k-1}\xi_{k-1}^{(a)} - \xi_{k}^p) - K_k\xi_k^{(m)}\\
\end{aligned}
```

and therefore the covariance matrix is 

```math
\begin{aligned}
    P_k &= \E[\xi_{k}^{(a)}(\xi_{k}^{(a)})^T] \\ 
        &= \E\Big[\left(\big(I-K_kH_k \big)(M_{k-1}\xi_{k-1}^{(a)} - \xi_{k}^p) - K_k\xi_k^{(m)} \right) \left(\big(I-K_kH_k \big)(M_{k-1}\xi_{k-1}^{(a)} - \xi_{k}^p) - K_k\xi_k^{(m)} \right)^T \Big] \\ 
        &= \big(I-K_kH_k \big)M_{k-1}\E[\xi_{k-1}^{(a)}(\xi_{k-1}^{(a)})^T]M_{k-1}^T\big(I-K_kH_k \big)^T + \big(I-K_kH_k \big)\E[\xi_k^{(p)}(\xi_k^{(p)})^T]\big(I-K_kH_k \big)^T - K_k\E[\xi_{k}^{(m)}(\xi_k^{(m)})^T]K_k^T \\ 
        &= \big(I-K_kH_k \big)B_k\big(I-K_kH_k \big)^T - K_kR_kK_k^T
\end{aligned}
```

##  Deriving ``K_k`` 
The Kalman filter is defined at that ``K_k`` which which minimizes the sum of squared analysis errors, i.e. the trace of the analysis error covariance matrix. The following identies will be useful: 
```math
\begin{aligned}
    \mathop{\nabla}_{A}\text{tr}(AB) &= B^T \\ 
    \mathop{\nabla}_{A}\text{tr}(BA^T) &= B \\ 
    \mathop{\nabla}_{A}\text{tr}(ABA^T) &= AB^T + AB  \\ 
\end{aligned}
```
from which we obtain 
```math
\begin{aligned}
    0 &= \mathop{\nabla}_{K_k}\text{tr}(P_k) \\ 
      &= \mathop{\nabla}_{K_k}\Big\{ B_k -B_kH_k^TK_k^T - K_kH_kB_k  + K_kH_kB_kH_k^TB_k^T - K_kR_kK_k^T \Big\} \\ 
      &= -B_kH_k^T - (H_kB_k)^T + K_k\left[H_kB_kH_k^T + (H_kB_kH_k^T)^T - R_k+R_k^T\right] \\
      &= -2B_kH_k^T + 2K_k\left(H_kB_kH_k^2 - R_k \right) \\ 
  \Rightarrow K_k &= B_kH_k^T\Big[ H_kB_kH_k^T - R_k \Big]^{-1}
\end{aligned}
```
we now substitute this result to obtain a simplified form for ``P_k``. 

```math
\begin{aligned}
    P_k &= \left(I - K_kH_k \right)B_k\left(I - K_kH_k \right)^T + K_kR_kK_k^T \\ 
        &= \left(I - K_kH_k \right)B_k - \left(I - K_kH_k \right)B_k\left(K_kH_k \right)^T + K_kR_kK_k^T \\
        &= \left(I - K_kH_k \right)B_k -\left\{ \left(I - K_kH_k \right)B_k\left(K_kH_k \right)^T + K_kR_kK_k^T \right\} \\
        &= \left(I - K_kH_k \right)B_k -\left\{ \left(I - K_kH_k \right)B_kH_k^TK_k^T + K_kR_kK_k^T \right\} \\
        &= \left(I - K_kH_k \right)B_k -\left\{ \left(I - K_kH_k \right)B_kH_k^T + K_kR_k \right\}K_k^T \\
        &= \left(I - K_kH_k \right)B_k -\left\{ B_kH_k^T - K_k\left( H_kB_kH_k^T + R_k \right)  \right\}K_k^T \\
        &= \left(I - K_kH_k \right)B_k -\left\{ B_kH_k^T - B_kH_k^T \right\}K_k^T \\
        &= \left(I - K_kH_k \right)B_k 
\end{aligned}
```
**NOTE**: we have used the fact that covariance matrices are symmetric. 

## Summary 
Let's summarize the whole process. We have 

### 0. Initialization 
We must set the system to some initial condition. This means we must define ``u_0^a`` and ``P_0``. We must also come up with a model for the process noise covariance ``Q_k`` and measurement error covariance ``R_k``.

### 1. Forecast Step 
```math
\begin{aligned}
    u_k^{(b)} &= M_{k-1}u_{k-1}^{(a)} \\ 
    B_k &= M_{k-1}P_{k-1}M_{k-1}^T + Q_{k}
\end{aligned}
```

### 2. Assimilation Step 
```math
\begin{aligned}
    K_k &= B_kH_k^T\Big[ H_kB_kH_k^T - R_k \Big]^{-1}  \\ 
    u_k^{(a)} &= u_k^{(b)} + K_k(w_k - H_ku_k^{(b)})\\ 
    P_k &= \left(I - K_kH_k \right)B_k 
\end{aligned}
```
