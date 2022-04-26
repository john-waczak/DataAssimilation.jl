using DifferentialEquations
using LinearAlgebra
using Statistics
using Distributions
using Plots


# eventually we want to have:
#
# AnalysisProblem(lorenz!, u0, tspan, p, Q0, R0,...)
#
# and then
#
# sol = solve(AnalysisProb)
#
# probably just pass the kwargs of AnalysProblem() through to the inner ODEProblem()



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
tspan = (0.0, 100.0)

# set up the ODE problem
prob = ODEProblem(lorenz!, u0, tspan, p)

# solve and plot
sol_true = solve(prob)
plot(sol_true, vars=(1,2,3))

savefig("true_sol_butterfly.png")


# define Jacobian ∂fᵢ/∂xⱼ
function Jf(u, p)
    x,y,z = u
    σ,ρ,β = p

    return [
        -σ σ 0;
        ρ-z -1 -x;
        y x -β
    ]
end

# test the jacobian
display(Jf(u0, p))



# define the observation function and it's Jacobian
# for simplicity, let's assume the observation
# is the full state vector (i.e. not nonlinear)
function h(u, p)
    return u
end

function Jh(u, p)
    return 1.0I
end


# test observation function and Jacobian
display(h(u0, p))
display(Jh(u0, p))



# define EKF analysis step
function EKF(xf, Pf, y, R, h, Jh, params)
    # compute linerized observation operator
    H = Jh(xf, params)

    # compute Kalman gain
    denom = H*Pf*H' + R
    K =Pf*H'*inv(denom)

    # compute analysis
    xₐ = xf + K*(y-h(xf, params))
    Pₐ = (I - K*H)*Pf

    return xₐ, Pₐ
end


# test the analysis
EKF(u0, 0.1^2 .* I(3), u0, 0.5^2 .* I(3), h, Jh, p) # should give back u0




# let's test it out

dt = 0.01
tspan = (0.0, 10.0)
u₀_true = [1.0; 1.0; 1.0]
params = [10.0, 28.0, (8/3)]


# set up observation vector
σ_m = 0.15  # measurement noise
R = σ_m^2 .* I(3)

dt_m = 0.2 # time between observations
tm_m = 2.0 # maximum time for observations
ts_obs = 0:dt_m:tm_m

y = zeros(size(ts_obs, 1))


# initalize array for true state values
u_true = zeros(3, size(tspan[1]:dt:tspan[2], 1))
u_true[:, 1] = u₀_true

# set up ode problem and solve it for true solution
prob_true = ODEProblem(lorenz!, u₀_true, tspan, params, saveat=dt)
sol_true = solve(prob_true)

# fill up the array
u_true[:, 2:end] .= sol_true[:, 2:end]

# fill up observation array
y = hcat([sol_true(t)+rand(MvNormal(zeros(3), R)) for t ∈ ts_obs]...)

# let's see how it worked
plot(tspan[1]:dt:tspan[2], u_true[1,:], linewidth=2, color="royal blue", alpha=0.75, label="truth")
plot!(ts_obs, y[1,:], seriestype=:scatter, label="obs")
xlabel!("t")
ylabel!("x(t)")


# Initialization
u₀_f = [2.0; 3.0; 4.0]
σ_f = 0.1
P₀ = σ_f^2 .* I(3)

# set up result
uₐ = zeros(3, size(tspan[1]:dt:tspan[2], 1))
uₐ[:, 1] = u₀_f

# build an array of times for the integration
times = [t for t ∈ ts_obs]
push!(times, tspan[2])


# loop through times and perform the analysis
for k ∈ 2:size(times, 1)
    Δt = (times[k-1], times[k])
    idx_i = findfirst(x->x==Δt[1], tspan[1]:dt:tspan[2])
    idx_f = findfirst(x->x==Δt[2], tspan[1]:dt:tspan[2])

    x₀ = uₐ[:, idx_i]

    prob = ODEProblem(lorenz!, x₀, Δt, params, saveat=dt)
    sol = solve(prob)

    # update uₐ with the analysis
    uₐ[:, idx_i:idx_f]  .= sol[:,:]

end


plot(tspan[1]:dt:tspan[2], u_true[1,:], linewidth=2, color="royal blue", alpha=0.75, label="truth")
plot!(tspan[1]:dt:tspan[2], uₐ[1,:], linewidth=2, color="red", alpha=0.75, label="analysis")
xlabel!("t")
ylabel!("x(t)")




