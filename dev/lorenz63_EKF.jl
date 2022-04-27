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
function lorenz!(du, u, θ, t)
    x,y,z = u
    σ,ρ,β = θ

    du[1] = σ*(y - x)
    du[2] = x*(ρ - z) - y
    du[3] = x*y - β*z
end



function EKF(xf, y, h, Jh, R, Pf, θ)
    # compute linerized observation operator
    H = Jh(xf, θ)

    # compute Kalman gain
    denom = H*Pf*H' + R
    K =Pf*H'*inv(denom)

    # compute analysis
    println("obs error:\t $(y .- h(xf, θ))")
    println("K matrix:\t $(K)")
    xₐ = xf + K*(y .- h(xf, θ))
    Pₐ = (I - K*H)*Pf

    return xₐ, Pₐ
end



# parameters
θ = [10.0, 28.0, 8/3]  # σ,ρ,β

dt = 0.01
tm = 10.0
t = 0:dt:tm


# define observation operator and the Jacobian functions
function h(u, θ)
    y = u
    return y
end


function Jh(u, θ)
    return I(size(u, 1))
end

function Jf(u, θ)
    x,y,z = u
    σ,ρ,β = θ
    jac = [
        -σ σ 0;
        ρ-z -1 -x;
        y x -β
    ]
    return jac
end



# 1. Generate true values
u0_true = [1.0; 1.0; 1.0]
prob = ODEProblem(lorenz!, u0_true, (t[1],t[end]), θ, saveat=dt)
sol_true = solve(prob)


# 2. Generate observations w/ noise
σ_obs = 0.15  # observation error
R = σ_obs^2 .* I(3)

dt_obs = 0.2 # observation dt
tm_obs = 2.0 # time of last observation
t_obs = 0:dt_obs:tm_obs

y = hcat([sol_true(tᵢ)+rand(MvNormal(zeros(3), R)) for tᵢ ∈ t_obs]...)


# 3. Data Assimilation
u0_b = [2.0, 3.0, 4.0]
σ_m = 0.1
P₀ = σ_m^2 .* I(3)
Q = 0.0 .* I(3)


# set up arrays for results
uf = zeros(3, size(t,1))
uₐ = zeros(3, size(t,1))
Pf = zeros(3, 3, size(t,1))
Pₐ = zeros(3, 3, size(t,1))

# set initial values
uₐ[:,1] = u0_b
Pₐ[:,:,1] = P₀

# set up the integrator
prob = ODEProblem(ff, u0_b, (t[1], t[end]), θ, dt=dt)
# integrator = init(prob, Euler())
integrator = init(prob, Tsit5())

# loop and perform the integration
k_obs = 1
for k ∈ 1:(size(t,1)-1)
    # step forward unitl we get to next time point
    step!(integrator, dt, true)
    uₐ[:, k+1] = integrator.u

    # compute jacobian of model
    Df = I(3) + dt*Jf(uₐ[:,k], θ)
    # Df = Jh(uₐ[:,k], θ)

    # propagate the error covariances
    Pₐ[:,:,k+1] = Df*Pₐ[:,:,k]*Df' + Q

    # if we have an observation, calc Kalman gain and update uₐ,Pₐ
    if any(integrator.t ≈ tᵢ for tᵢ ∈ t_obs)

        H = Jh(uₐ[:, k+1], θ)
        # compute Kalman gain
        denom = H*Pₐ[:,:,k+1]*H' + R
        K =Pₐ[:,:,k+1]*H'*inv(denom)

        uₐ[:,k+1] = uₐ[:,k] + K*(y[k_obs] .- h(uₐ[:,k+1], θ))
        Pₐ[:,:,k+1] = (I - K*H)*Pₐ[:,:,k+1]

        global k_obs += 1
    end

    # update integrator
    integrator.u = uₐ[:, k+1]
end


p = plot(sol_true.t, sol_true[1,:], label="true")
plot!(t_obs, y[1,:], seriestype=:scatter, label="obs")
plot!(t, uₐ[1,:], label="analysis")
#xlims!(0, 2)
ylims!(-20,20)
display(p)

