# Data Assimilation Overview

## Overview

The proper application of scientific models to make real-world predictions requires that we commit ourselves to a full accounting of all possible sources of uncertainty when reporting results. Further, the explosion of *big data* across scientific fields provieds a plethora observational data that our models are typically unequipped to incorporate when making predictions. The field of **Data Assimilation** addresses this problem by providing a family of techniques engineered to combine model output together with observational data whilst enabling a complete accounting the sources of uncertainty. For chaotic systems in particular, data assimilation enables integration on long time scales that would be impossible via models alone. 


In this overview, we will follow the examples from [this nice paper](https://www.mdpi.com/2311-5521/5/4/225/htm). 

## Framing the Problem
Data assimilation can be understood most generally in terms of dyscrete dynamical systems. This enables us to apply the methods to most mathematical models from gridded PDE solvers to systems of ordinary differential equations. Our goal is to find the best prediction for the system state vector ``u`` that combines our model predictions, also known as forecasts, with observational data. Model predictions are summarized via the discrete update equation: 

```math
u_{k+1} = \mathcal{M}(u_k; \theta)
```

For ODE systems, ``\mathcal{M}`` represents the time integration scheme for a system of ODEs like 

```math
\dfrac{du}{dt} = f(u, t; \theta)
```

To measure the performance of our assimilation scheme, we denote the *true* value of the state vector asd ``u^{(t)}``. The output of our model is denoted ``u^{(b)}`` (*b* subscript for *background*). The discrepancy between the true value and our forecast is denoted ``\xi^{(b)} = u^{(t)} - u^{(b)}`` characterizing the extent to which our model prediction is imperfect. 

The observations of our system are denoted by ``w_k = w(t_k)``. These observations do not neccessarily need to be components of the state vector ``u``, but rather, are related to it via the *observation function* ``h``. For example, one may attempt to predict sea surface temperature using data assimilation with data from satellite observations. The function ``h`` would then be the Stefan-Boltzmann law. However, real world data is noisy, which we must take into account. We write 
```math
w_k = h(u_k) + \xi_k^{(m)}
```
where ``\xi_k^{(m)}`` denotes this measurement noise. 


Given our model predictions ``u_{k}^{(b)}`` and observations ``w_k``, we seek to obtain the *optimal* or best-possible prediction called the **analysis**, ``u^{(a)}``. This analysis will still not be perfect, so we further specify the analysis error via 
```math
\xi^{(a)} = u^{(t)} - u^{(a)}
```

## Summary 

```math
\begin{aligned}
    &u_k^{(t)} \in \R^n &\text{the true state vector} \\ 
    &u_k^{(b)} \in \R^n &\text{the k^{th} model forecast} \\ 
    &u_k^{(a)} \in \R^n &\text{the analysis} \\ 
    &w_k \in \R^m &\text{the k^{th} observation vector} \\ 
    &\xi^{(b)} \in \R^n &\text{the model forecast error}\\
    &\xi^{(m)} \in \R^m &\text{the observation noise vector}\\ 
    &\xi^{(a)} \in \R^n &\text{the analysis error}\\
    &\xi^{(p)} \in \R^n &\text{the process noise if we used our model on the true state}\\
    &\mathcal{M}:\R^n\to\R^n &\text{the model update function}\\
    &f:\R^n\to\R^n &\text{differential equation model}\\ 
    &h:\R^n\to\R^m  &\text{observation function}
\end{aligned}
```

## Assumptions
To make possible the derivation of a *unique* analysis ``u^{(a)}``, the following assumptions are in order.
```math
\begin{aligned}
    &\E[\xi_k^{(b)}] = 0 & &\E[\xi_k^{(b)}(\xi_j^{(b)})^T] = 0 \text{ for } k\neq j\\
    &\E[\xi_k^{(m)}] = 0 & &\E[\xi_k^{(m)}(\xi_j^{(m)})^T] = 0 \text{ for } k\neq j\\
    &\E[\xi_k^{(b)}(u_0)^T] = 0 & &\E[\xi_k^{(m)}(u_0)^T] = 0\\
    &\E[\xi_k^{(b)}\xi_j^{(m)}] = 0 & &  \\ 
    &\E[u_k^{(t)}] = u_k^{(b)} & &
\end{aligned}
```

We also define the error covariance matrices
```math
\begin{aligned}
    Q_k &:= \E[\xi_k^{(p)}(\xi_k^{(p)})^T] \\ 
    R_k &:= \E[\xi_k^{(m)}(\xi_k^{(m)})^T] \\ 
    B_k &:= \E[\xi_k^{(b)}(\xi_k^{(b)})^T] 
\end{aligned}
```

which we will use in our consideration of the final error of our analysis.
