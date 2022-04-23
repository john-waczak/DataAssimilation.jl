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
f(x_{k-1}) \approx f(x_{k-1}^a) + J_f(x_{k-1}^a)(x_{k-1}-x_{k-1}^a) + \text{ Higher Order Terms}
```

where 
```math
J_f := \left[ \dfrac{\partial f_i}{\partial x_j}\right] 
```

Therefore we can acheive the approximation 
```math
\begin{aligned}
    e_k^f &\approx f(x_{k-1}^a) + J_f(x_{k-1}^a)e_{k-1} + w_{k-1} - f(x_{k-1}^a) \\ 
        &= J_f(x_{k-1}^a)e_{k-1} + w_{k-1}
\end{aligned}
```
We can now use this error vector to form the forecast error covariance matrix
```math
\begin{aligned}
    P_k^f &:= \E[e_k^f(e_k^f)^T] \\ 
        &\approx \E[(J_f e_{k-1} + w_{k-1})(J_f e_{k-1} + w_{k-1})^T] \\ 
        &\approx \E[(J_f e_{k-1} + w_{k-1})(e_{k-1}^T J_f^T  + w_{k-1}^T)] \\ 
        &\approx J_f\E[e_{k-1}e_{k-1}^T]J_f^T + J_f\E[e_{k-1}w_{k-1}^T] + \E[w_{k-1}e_{k-1}^T]J_f^T + \E[w_{k-1}w_{k-1}^T] \\ 
        &= J_f(x_{k-1}^a)P_{k-1}J_f^T(x_{k-1}^a) + Q_{k-1}
\end{aligned}
```

## Deriving Data Assimilation Step
At time index ``k`` we have ``x_k^f``, ``P_k^f``, ``z_k``, and ``R_k``. We will now use these to derive the optimal analysis ``x_k^a``. 

Let's assume that the result is of the form
```math
x_k^a = a + K_kz_k \;\; \in \R^n
```
with ``a\in\R^n`` and ``K_k\in\R^{m\times n}``. We desire that ``\E[x_k - x_k^a] = 0``. Thus, 
```math
\begin{aligned}
    0 &= \E[x_k - x_k^a] \\ 
    &= \E[(x_k^f + e_k^f) - (a+K_k z_k)] \\ 
    &= \E[(x_k^f + e_k^f) - (a+K_k h(x_k)  +K_kv_k)] \\ 
    &= \E_{x_k}[x_k^f] + \E_{x_k}[e_k^f] - \E_{x_k}[a] - K_k\E_{x_k}[h(x_k)] - K_k\E_{x_k}[v_k] \\ 
    &= x_k^f + 0 - a - K_k\E_{x_k}[h(x_k)]  - 0 \\ 
    &= x_k^f - a - K_k\E[h(x_k)] \\ 
    \Rightarrow a &= x_k^f - K_k\E[h(x_k)]
\end{aligned}
```
We now substitute to obtain

```math
x_k^a - x_k^f - K_k\E[h(x_k)] + K_kz_k
```
Since ``h`` is smooth enough, we may expand it about ``x_k^f`` to find that 
```math
h(x_k) \approx h(x_k^f) + J_h(x_k^f)(x_k-x_k^f)
```
where 
```math
J_h := \left[ \dfrac{\partial h_i}{\partial x_j} \right]
```
This leads to the approximation 
```math
\begin{aligned}
    \E[h(x_k)] &\approx \E[h(x_k^f) + J_h(x_k^f)e_k^f] \\ 
        &= \E[h(x_k^f)] + J_h\E[e_k^f] \\ 
        &= \E[h(x_k^f)] = h(x_k^f)
\end{aligned}
```
The analysis then becomes 
```math
\begin{aligned}
    x_k^a &= x_k^f - K_k\E[h(x_k)] + K_k z_k \\ 
        &\approx x_k^f - K_kh(x_k^f) + K_kz_k \\ 
        &= x_k^f + K_k(z_k - h(x_k^f))
\end{aligned}
```


## Determining the Filter Matrix ``K_k``
To determine the ``K_k`` that will give us the optimal ``x_k^a``, we will first derive the analysis error covariance matrix and then optimize its trace (i.e. the sum of squared errors) with respect to ``K_k``. Following the same procedure as before, we begin by infestigating the error ``e_k=x_k-x_k^a``.
```math
\begin{aligned}
    e_k &= x_k - x_k^a \\ 
    &= f(x_{k-1}) + w_{k-1} - x_k^f - K_k(z_k-h(x_k^f)) \\ 
    &\approx f(x_{k-1}) + w_{k-1} - f(x_{k-1}^a) - K_k(h(x_k) - h(x_k^f) + v_k) \\ 
    &\approx J_f(x_{k-1}^a)e_{k-1} + w_{k-1} - K_kJ_h(x_k^f)(J_f(x_{k-1}^a)e_{k-1} + w_{k-1}) - K_kv_k \\ 
    &= J_f(x_{k-1}^a) e_{k-1} - K_k J_h(x_k^f)J_f(x_{k-1}^a)e_{k-1} + w_{k-1} - K_kJ_h(x_k^f) w_{k-1} - K_kv_k \\ 
    &= \left(I-K_kJ_h(x_k^f) \right)J_f(x_{k-1}^a)e_{k-1}  + \left(I - K_kJ_h(x_k^f)\right)w_{k-1} - K_kv_k
\end{aligned}
```
We now form the covariance matrix
```math
\begin{aligned}
    P_k &:= \E[e_ke_k^T] \\ 
    &= \left( I - K_k J_h(x_k^f)\right)J_f(x_{k-1}^a)P_{k-1}J_f^T(x_{k-1}^a)\left( I - K_k J_h(x_k^f)\right)^T + \left( I - K_k J_h(x_k^f)\right)Q_{k-1}\left( I - K_k J_h(x_k^f)\right)^T + K_kR_kK_k^T \\ 
    &= \left( I - K_k J_h(x_k^f)\right)\left[ J_f(x_{k-1}^a)P_{k-1} J_f^T(x_{k-1}^a) + Q_{k-1}  \right]\left( I - K_k J_h(x_k^f)\right)^T + K_kR_kK_k^T \\ 
    &= \left( I - K_k J_h(x_k^f)\right)P_k^f \left( I - K_k J_h(x_k^f)\right)^T + K_kR_kK_k^T
\end{aligned}
```
We now seek to minimize ``\text{tr}(P_k)`` with respect to ``K_k``. The following identities will be helpful: 
```math
\begin{aligned}
    \nabla_{A}\text{tr}(AB) &= B^T \\ 
    \nabla_{A}\text{tr}(BA) &= B \\ 
    \nabla_{A}\text{tr}(ABA) &= AB^T + AB  \\ 
\end{aligned}
```
From these, we obtain 
```math
\begin{aligned}
    0 &= \nabla_{K_k}\text{tr}(P_k) \\ 
        &= -\left(J_h(x_k^f) P_k^f \right)^T - P_k^f\left( J_h(x_k^f) \right)^T + K_k\left( \left\{ J_h(x_k^f)P_k^f J_h^T(x_k^f)\right\}^T  + J_h(x_k^f P_k^f J_h^T(x_k^f))\right) + K_k(R_k^T + R_k) \\ 
        &= -2P_k^f J_h^T(x_k^f) + 2K_k\left[ J_h(x_k^f) P_k^f J_h^T(x_k^f) + R_k \right] \\ 
    K_k &= P_k^f J_h^T(x_k^f)\left[ J_h(x_k^f)P_k^f J_h^T(x_k^f) + R_k \right]^{-1}
\end{aligned}
```
