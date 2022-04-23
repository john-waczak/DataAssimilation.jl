# Data Assimilation Overview

## Overview

The proper application of scientific models to make real-world predictions requires that we commit ourselves to a full accounting of all possible sources of uncertainty when reporting results. Further, the explosion of *big data* across scientific fields provieds a plethora observational data that our models are typically unequipped to incorporate. The field of **Data Assimilation** addresses this problem by providing a family of techniques engineered to conbine model output together with observational data whilst enabling a complete accounting the sources of uncertainty. For chaotic systems in particular, data assimilation enables integration on long time scales that would be impossible via models alone. 

## Problem Framing
Data Assimilation can be easily understood in terms of discrete dynamical systems that involve both a model function ``f(\cdot)`` as well as observations ``z``. This can be summarized by: 
```math
\begin{aligned}
    x_k &= f(x_{k-1}) + w_{k-1} \\ 
    z_k &= h(x_k) + v_k
\end{aligned}
```
where 
```math
\begin{aligned}
    & x_k\in \R^n &\text{the state vector} \\ 
    & w_k \in \R^n &\text{model error} \\ 
    & z_k \in \R^m &\text{observation vector} \\ 
    & v_k \in \R^m &\text{observation error}\\ 
    & f:\R^n \to \R^n &\text{Model function (aka forecast function)} \\ 
    & h:\R^n \to \R^m &\text{Observation function} \\ 
    & x_k^a \in \R^n &\text{the analysis vector; our best estimate} \\ 
    & x_k^f := f(x_k^a) \in \R^n &\text{model output; forecast vector}\\ 
\end{aligned}
```

In words, we say that the current state of the system is given by our model acting on the previous state *plus* some uncertainty in the previous state captured by ``w_{k-1}``. At the same time we may make a set of observations captured by ``z_k`` which come directly from the model state via a function ``h(\cdot)`` and an associated uncertainty ``v_k``. 

## Initial Conditions 
To simulate the model, we must supply a vector of initial conditions, say ``x_0``. We also supply the initial mean and covariance matrix, denoted
```math
\begin{aligned}
    \mu_0 &= \E[x_0] & \in R^n \\
    P_0 &= \E\left[(x_0-\mu_0)(x_0-\mu_0)^T\right] & \in \R^{n\times n}
\end{aligned}
```

Here ``\E[\cdot]`` denotes the *expectation value*. 

## Assumptions 
We assume that the model and observations uncertainties have mean zero. Further, we assume that they are neither correleted with each other nor with the initial state vector. Finally, we make the (reasonable?) assumption that uncertainties have are not correlated between time steps, i.e 
```math
\begin{aligned}
    &\E[w_k] = 0 & &\E[w_kw_j^T] = 0 \text{ for } k\neq j \\ 
    &\E[v_k] = 0 & &\E[v_kv_j^T] = 0 \text{ for } k\neq j \\ 
    &\E[v_kx_0^T] = 0 & &\E[v_k x_0^T] = 0 \\ 
    &\E[w_kv_j^T] = 0 & & 
\end{aligned}
```

We also define the error covariance matrices 
```math
\begin{aligned}
    Q_k &:= \E[w_kw_k^T] \\ 
    R_k &:= \E[v_kv_k^T]
\end{aligned}
```
which we will use in our consideration of the final error of our analysis. 


## Goal 
The goal of data assimilation is to optimally combine our model forecast ``x_k^f`` with our observations ``z_k`` to acheive an analysis ``x_k^a`` which minimizes the error ``e_k = x_k - x_k^a`` between the true value ``x_k`` and our prediction ``x_k^a``. 

In the next section we will examine a popular approach called the **Extended Kalman Filter**.
