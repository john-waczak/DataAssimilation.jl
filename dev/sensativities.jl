using DiffEqSensitivity
using DifferentialEquations
using ForwardDiff
using DiffResults
using Zygote
using Plots
using FiniteDifferences





# define function for RHS of lorenz system
function lorenz!(du, u, p, t)
    x,y,z = u
    σ,ρ,β = p

    du[1] = σ*(y - x)
    du[2] = x*(ρ - z) - y
    du[3] = x*y - β*z
end

# set up initial conditions for truth
u0 = [1.0; 0.0; 0.0]
p = [10.0, 28.0, (8/3)]
# tspan = (0.0, 100.0)
tspan = (0.0, 0.5)

# set up the ODE problem
prob = ODEProblem(lorenz!, u0, tspan, p)


# solve and plot
sol = solve(prob, Tsit5(), reltol=1e-6,abstol=1e-6)


# integrate model from k to k+1
function model_forward(uₖ)
    _prob = remake(prob, u0=uₖ)
    sol = solve(_prob, Tsit5(), reltol=1e-6,abstol=1e-6)
    return sol[:, end]
end




# # test that it works
J_forwarddiff = ForwardDiff.jacobian(model_forward, u0)

u_test = [1.0, 2.0, 2.0]

res = DiffResults.JacobianResult(u_test)
ForwardDiff.jacobian!(res, model_forward, u_test)
DiffResults.value(res)
DiffResults.jacobian(res)



# ----- Now try with zygote and compare
J_reversediff = Zygote.jacobian(model_forward, u0)[1]

res = Zygote.withjacobian(model_forward, u_test)
res.val
res.grad[1]


# ---- Let's verify with ForwardDiff
J_finitediff = FiniteDifferences.jacobian(central_fdm(5,1), model_forward, u0)[1]

J_forwarddiff ≈ J_reversediff

