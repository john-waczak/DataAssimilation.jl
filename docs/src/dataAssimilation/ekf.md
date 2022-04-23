# Extended Kalman Filter

For the **Extended Kalman Filter** we seek to find the ``x_k^a`` that minimizes ``e_k = x_k - x_k^a`` in the least-squares sense. 

## Initialization

To begin, we apply the *use-what-you-have* strategy and set 
```math
\begin{aligned}
    x_0^a &= \mu_0 \\ 
    P_0 &= \E[(x_0-x_0^a)(x_0-x_0^a)^T]
\end{aligned}
```

## Deriving the Forecast Covariance 
Consider now the forecast error produced by our model ``f`` at time index ``k``. We have 
```math
\begin{aligned}
    e_k^f &:= x_k - x_k^f \\ 
        &= f(x_{k-1}) + w_{k-1} - f(x_{k-1}^a) \\ 
\end{aligned}
```
If our functions ``f`` and ``h`` are *smooth enough* (in this case, ``C^1``), we may expand in a Taylor series and determine their linearization via the Jacobian about the analysis vector ``x_{k-1}^a``. That is, 
```math
\begin{equation}
    f(x_{k-1}) \appox f(x_{k-1}^a) + J_f(x_{k-1}^a)(x_{k-1}-x_{k-1}^a) + \text{ Higher Order Terms}
\end{equation}
```

where 
```math
\begin{equation}
    J_f := \left[ \dfrac{\partial f_i}{\partial x_j}\right] = \begin{bmatrix}
        \frac{\partial f_1}{\partial x_1} & \cdots & \frac{\partial f_1}{\partial x_n} \\ 
        \vdots & \ddots & & \\ 
        \frac{\partial f_n}{\partial x_1} & \cdots & \frac{\partail f_n}{\partial x_n}
    \end{bmatrix}
\end{equation}
```
