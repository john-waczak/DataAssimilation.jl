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
    & x_k\in \R^n & \text{ the state vector} \\ 
    & w_k \in \R^n & \text{ model uncertainty} \\ 
    & z_k \in \R^m & \text{ observation vector} \\ 
    & v_k \in \R^m & \text{ observation uncertainty}\\ 
    & f:\R^n \to \R^n & \text{ Model function (aka forecast function)} \\ 
    & h:\R^n \to \R^m & \text{ Observation function}
\end{aligned}
```
