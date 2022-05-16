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
