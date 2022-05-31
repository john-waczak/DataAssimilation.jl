using DataAssimilation
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

# define observation operator jacobian. We should probably generate this symbolically
function Dh(u)
    n = size(u, 1)
    D = I(n)
    return D
end



# populate observation vector
w = hcat([h(sol_true(tᵢ))+rand(MvNormal(zeros(3), R)) for tᵢ ∈ ts_m]...)

# set up data assimilation
u0b = [2.0; 3.0; 4.0]  # guess at initial condition
σ_b = 0.1  # model error std
B = σ_b^2 .* I(3)  # initial model error covariance
Q = 0.0 * I(3)  # process noise (assumed zero)


# integrate to obtain the true solution
prob_no_analysis = ODEProblem(lorenz!, u0b, ts_sim, θ, saveat=dt)
sol_no_analysis = solve(prob_no_analysis)


fourdvar_prob = AssimilationProblem(u0b,
                                    ts_sim,
                                    lorenz!,
                                    θ,
                                    dt,
                                    Q,
                                    B,
                                    R,
                                    dt_m,
                                    h,
                                    Dh,
                                    w
                                    )
u0a = FourDVarSolve(fourdvar_prob, 500, ADAM(0.25))



prob = ODEProblem(lorenz!, u0a, (0.0, tend), θ)
sol_analysis = solve(prob)



var_names = ["x(t)", "y(t)", "z(t)"]
plots = []
for i ∈ 1:3
    p = plot(sol_true, vars=(0,i), color=:black, alpha=0.75, linewidth=2, label="true system")
    plot!(p, ts_m, w[i, :], seriestype=:scatter, color="light green", label="observation" )
    plot!(p, sol_no_analysis, vars=(0, i), linestyle=:dashdot, color="purple", label="bad forecast")
    plot!(p, sol_analysis, vars=(0,i), color="royal blue", alpha=1, linewidth = 1, label="analysis")
    xlabel!(p, "t")
    ylabel!(p, var_names[i])
    push!(plots, p)
end

final_plot = plot(plots..., layout=(3,1))

savefig("4dvar_using_flux.svg")
savefig("4dvar_using_flux.png")

display(final_plot)
