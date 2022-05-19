using DifferentialEquations
using DiffEqSensitivity
using LinearAlgebra
using Zygote
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


# visualize the data
px = plot(sol_true, vars=(0,1), color=:black, alpha=0.75, linewidth=2, label="true system")
xlabel!("t")
ylabel!("x(t)")
py = plot(sol_true, vars=(0,2), color=:black, alpha=0.75, linewidth=2, label="true system")
xlabel!("t")
ylabel!("y(t)")

pz = plot(sol_true, vars=(0,3), color=:black, alpha=0.75, linewidth=2, label="true system")
xlabel!("t")
ylabel!("z(t)")

train_plot = plot(px, py, pz, layout=(3,1))




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

# define observation operator jacobian. We should probably generate this symbolically
function Dh(u)
    n = size(u, 1)
    D = I(n)
    return D
end



# populate observation vector
w = hcat([h(sol_true(tᵢ))+rand(MvNormal(zeros(3), R)) for tᵢ ∈ ts_m]...)


# update the graphs with observation stuff
plot!(px, ts_m, w[1, :], seriestype=:scatter, color="light green", label="observation" )
plot!(py, ts_m, w[2, :], seriestype=:scatter, color="light green", label="observation" )
plot!(pz, ts_m, w[3, :], seriestype=:scatter, color="light green", label="observation" )

train_plot = plot(px, py, pz, layout=(3,1))




## ----- perform the assimilation ------------------------------------

# set up data assimilation
u0b = [2.0; 3.0; 4.0]  # guess at initial condition
σ_b = 0.1  # model error std
B = σ_b^2 .* I(3)  # initial model error covariance
Q = 0.0 * I(3)  # process noise (assumed zero)

# set up arrays for solution outputs... we should probably pre-allocate these
ua = zeros(3, size(0:dt:tend, 1))
ua[:,1] = u0b


# function for getting indices from time value
function timeindex(t, Δt)
    return Int(round(t/Δt)) + 1
end


# helper function to update the analysis result
function update!(sol, t_now, t_next, dt)
    idx_now = timeindex(t_now, dt)
    idx_next = timeindex(t_next, dt)

    ua[:, idx_now:idx_next] = hcat(sol(t_now:dt:t_next)...)
end
Zygote.@nograd update!  # tell zygote to ignore update! when computing derivatives


function model_forward!(uₖ)
    prob = ODEProblem(lorenz!, uₖ, (t_now, t_next), θ)
    sol = solve(prob, Tsit5(), reltol=1e-6,abstol=1e-6)

    update!(sol, t_now, t_next, dt)

    return sol[:, end]
end


# ------------- perform the analysis ---------------------------------------------
# times = hcat(collect(0:dt_m:tm_end)..., tend)
times = 0:dt_m:tend


t_now = times[1]
t_next = times[2]
model_forward!(u0b)

for k ∈ 1:(length(times)-1)
    t_now = times[k]
    t_next = times[k+1]


    idx_now = timeindex(t_now, dt)
    u_now = ua[:, idx_now]

    u_next, DM = model_forward!(u_now)
    DM = DM[1]

    # update error covariance matrix
    B = DM * B * DM' + Q

    if t_next ∈ ts_m
        println("doing the update!")
        n = size(ua, 1)
        idx_next = timeindex(t_next, dt)
        idx_h = timeindex(t_next, dt_m) - 1 # -1 since this doesn't include 0

        # compute the Kalman gain matrix
        DH = Dh(ua[:,idx_next])
        println("\t Dₕ=$(DH)")
        denom = DH*B*DH' + R
        println("\t denom:= $(denom)")
        K = B*DH'*inv(denom)

        println("\tK= $(K)")
        # perform the analysis

        println("\tua = $(ua[:, idx_next])")
        ua[:,idx_next] .= ua[:, idx_next] + K*(w[:, idx_h] - h(ua[:, idx_next]))

        println("\tua = $(ua[:, idx_next])")
        B = (I(n) - K*DH)*B
    end
end


plot(0:dt:tend, ua[1,:])
