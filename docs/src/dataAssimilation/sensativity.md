# Sensativity Analysis for Differential Equations

Provided some model for a physical system in the form of a set of differential equations, a natural question is: How can we select the parameters for our model in order to get the best possible fit to some experimental data. Similarly, one may wonder what would happen to the prediction of your model if you were to slightly change the values of some parameters. In other words, how *sensative* is the output of our model to your choice of parameter values? 

In the most general sense, we may frame the problem as follows. Suppose we have a model of the form 
```math
\dfrac{du}{dt} = f(u,t,\theta)
```
Given this model, our goal is to optimize a cost function 
```math
J(u; \theta) := \int_0^T g(u;\theta)dt
```
where ``g(u;\theta)`` is usually taken to be some *quadratic form*. 

As an example, we might consider ``g(u\; \theta) = (u(t)-w(t))^T(u(t)-w(t))`` where ``w(t)`` denotes the vector of observations at time ``t``. 


Our goal then is to find out how ``J`` depends on the parameters ``\theta``, in other words, to find ``\partial J / \partial \theta``. To do this, we will use the method of Lagrange multipliers to generate a so called *adjoint equation* that enables us to find this derivative in a way that minimizes computational cost. As always, this method begins by adding a term that evaluates to 0 into our cost function: 
```math
    \mathcal{L} := \int_0^T \left[ g(u;\theta) + \lambda^T(t)\left(f-\dfrac{du}{dt}\right) \right] dt
```
From this, we find 
```math
\begin{aligned}
    \dfrac{\partial \mathcal{L}}{\partial \theta} &:= \int_0^T\left[ \frac{\partial g}{\partial \theta} + \frac{\partial g}{u}\frac{\partial u}{\partial \theta} + \lambda^T(t)\left( \frac{\partial f}{\partial \theta} + \frac{\partail f}{\partial u}\frac{\partail u}{\partial \theta} - \frac{d}{dt}\frac{\partail u}{\partial \theta} \right)\right]dt \\ 
    &= \int_0^T \left[ \frac{\partial g}{\partial \theta} + \lambda^T(t)\frac{\partial f}{\partial \theta} + \left( \frac{\partial g}{\partial u} + \lambda^T(t)\frac{\partial f}{\partial u} - \lambda^T(t)\frac{d}{dt} \right)\frac{\partial u}{\partial \theta} \right]dt
\end{aligned}
```
This reorganization is nice because the term ``\partial u/\partial \theta`` is the one thats *hard* to compute. Therefore, if we can make the terms in the paretheses evaluate to 0, we will be able to remove this pesky term. Let's use integration by parts to further rearrange by moving the ``d/dt``.

```math
\begin{aligned}
    \int_0^T-\lambda^T(t)\frac{d}{dt}\frac{\partial u}{\partial \theta} dt &= \left[-\lambda^T(t)\frac{\partial u}{\partial \theta} \right]_0^T_+ \int_0^T \frac{d\lambda^T(t)}{dt}\frac{\partial u}{\partial \theta}dt \\ 
    &= \lambda^T(0)\frac{\partial u_0}{\partial \theta} - \lambda^T(T)\frac{\partial u(T)}{\partial \theta} + \int_0^T \left[ \frac{d\lambda}{dt} \right]^T\frac{\partial u}{\partial \theta}dt
\end{aligned}
```
so that plugging this back into our expression for ``\partial \mathcal{L}//\partial \theta``, we obtain
```math
\frac{\partial \mathcal{L}}{\partial \theta} = \int_0^T \left[ \frac{\partial g}{\partial \theta} + \lambda^T\frac{\partial f}{\partial \theta} + \left( \frac{\partial g}{\partial u} + \lambda^T\frac{\partial f}{\partial u} + \left[\frac{d\lambda}{dt}\right]^T \right)\frac{\partial u}{\partial \theta}\right]dt + \lambda^T(0)\frac{\partial u_0}{\partial \theta} - \lambda^T(T)\frac{\partial u(T)}{\partial \theta}
```
Thus, forcing the nasty terms to dissappear is equivalent find the ``\lambda(t)`` subject to the differential equations 
```math
\begin{aligned}
   \frac{\partial g}{\partial u} + \lambda^T(t)\frac{\patial f}{\partial u} + \frac{d\lambda^T(t)}{dt} &= 0 \\ 
   \lambda^T(T) &= 0
\end{aligned}
```
or by taking the transpose: 
```math
\begin{aligned}
    \frac{d}{dt}\lambda &= - \left[ \frac{\partial g}{\patial u} \right]^T - \left[ \frac{\partial f}{\partial u} \right]^T\lambda  \\ 
    \lambda(T) &= 0
\end{aligned}
```

## Summary
To find the sensativities ``\partial J/\partial \theta``, we perform the following: 
1. Integrate the model ``du/dt = f(u,t,\theta)`` forward to obtain ``u(t)``.
2. Integrate the adjoint model ``d\lambda/dt = -(\partial f/ \patial u)^T\lambda - (\partial g / partial u)^T `` backwards in time from ``T`` to ``0`` to obtain ``\lambda(t)``.
3. Evaluate ``\partial J / \partial \theta = \int_0^T\left( \partial g/ partial \theta + \lambda^T \partial f/\partial \theta\right)dt + \lambda^T(0)\partial u_0/\partial \theta  ``
