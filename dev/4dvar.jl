using Flux
using DiffEqSensitivity
using DifferentialEquations
using LinearAlgebra
using Statistics
using Distributions
using Plots
using Random

# set random seed for reproducible results
Random.seed!(1)



## ----- setup true solution ------------------------------------


# define function for RHS of lorenz system
function lorenz!(du, u, p, t)
    x,y,z = u
    σ,ρ,β = p

    du[1] = σ*(y - x)
    du[2] = x*(ρ - z) - y
    du[3] = x*y - β*z
end

θ = [10.0, 28.0, (8/3)]  # true parameter values
dt = 0.01  # time step for simulation (this may be unnecessary)
tend = 10.0  # max time for total simulation
ts_sim = (0.0, tend)
u0True = [1.0, 1.0, 1.0]  # true initial condition

# integrate to obtain the true solution
prob_true = ODEProblem(lorenz!, u0True, ts_sim, θ, saveat=dt)
sol_true = solve(prob_true)



## ----- setup observation data ------------------------------------

# set up observation uncertainties
σ_m = 0.15  # i.e. observation std
# set up observation error covariance matrix
R = σ_m^2 .* I(3)

# set up details for measurement points
dt_m = 0.2
tm_end = 2.0  # final observation time
ts_m = dt_m:dt_m:tm_end

# define observation operator
function h(u)
    w = u
    return w
end


# populate observation vector
w = hcat([h(sol_true(tᵢ))+rand(MvNormal(zeros(3), R)) for tᵢ ∈ ts_m]...)



## ----- perform the assimilation ------------------------------------

# initial condition guess
u0a = [2.0, 3.0, 4.0]

# set up problem for (bad) forecast
prob = ODEProblem(lorenz!, u0a, (0.0,tend), θ)
sol_forecast = solve(prob)




function timeindex(t, Δt)
    return Int(round(t/Δt)) + 1
end



# for some reason I need to collect the matrices first
Rinv = collect(inv(R))
u0b = copy(u0a)
σ_b = 0.1  # model error std
B = σ_b^2 .* I(3)  # initial model error covariance
Binv = collect(inv(B))

function loss()
    res = Array(solve(ODEProblem(lorenz!, u0a, (0.0,tm_end), θ), saveat=dt_m:dt_m:tm_end))
    # ℓ = 0.5*(u0a-u0b)'*Binv*(u0a-u0b) + 0.5*sum((w[:,i]-h(res[:,i]))'*Rinv*(w[:,i]-h(res[:,i])) for i ∈ 1:size(res,2))
    ℓ = 0.5*sum((w[:,i]-h(res[:,i]))'*Rinv*(w[:,i]-h(res[:,i])) for i ∈ 1:size(res,2))
    return ℓ
end


u0a = [2.0, 3.0, 4.0]
loss()

# let's do the optimization
data = Iterators.repeated((), 500)
opt = ADAM(0.1)

cb = function ()
    display(loss())
end

Flux.train!(loss, Flux.params(u0a), data, opt, cb=cb)


u0a

prob = ODEProblem(lorenz!, u0a, (0.0, tend), θ)
sol_analysis = solve(prob)



var_names = ["x(t)", "y(t)", "z(t)"]
plots = []
for i ∈ 1:3
    p = plot(sol_true, vars=(0,i), color=:black, alpha=0.75, linewidth=2, label="true system")
    plot!(p, ts_m, w[i, :], seriestype=:scatter, color="light green", label="observation" )
    plot!(p, sol_forecast, vars=(0, i), linestyle=:dashdot, color="purple", label="forecast")
    plot!(p, sol_analysis, vars=(0,i), color="royal blue", alpha=1, linewidth = 1, label="analysis")
    xlabel!(p, "t")
    ylabel!(p, var_names[i])
    push!(plots, p)
end

plot(plots..., layout=(3,1))

savefig("4dvar_using_flux.svg")
savefig("4dvar_using_flux.png")


