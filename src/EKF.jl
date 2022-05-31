using DifferentialEquations

using LinearAlgebra

using DiffEqSensitivity
using Zygote

using Statistics
using Distributions


include("AssimilationProblem.jl")


# define function to get time index given a Δt
function timeindex(t, Δt)
    return Int(round(t/Δt)) + 1
end


# helper function to update the analysis result
function update!(sol, t_now, t_next, Δt, uₐ)
    idx_now = timeindex(t_now, Δt)
    idx_next = timeindex(t_next, Δt)

    uₐ[:, idx_now:idx_next] = hcat(sol(t_now:Δt:t_next)...)
end
Zygote.@nograd update!  # tell zygote to ignore update! when computing derivatives



function model_forward!(uₖ, uₐ, t_now, t_next, θ)
    prob = ODEProblem(lorenz!, uₖ, (t_now, t_next), θ)
    sol = solve(prob, Tsit5(), reltol=1e-6,abstol=1e-6)

    update!(sol, t_now, t_next, dt, uₐ)

    return sol[:, end]
end


function EKFsolve(prob::AssimilationProblem)
    # 1. allocate solution arrays for u(t) and B(t)
    uₐ = zeros(size(prob.u₀, 1), size(tspan[1]:prob.Δt:tspan[2], 1))
    B = zeros(size(prob.u₀,1), size(prob.u₀,1), size(tspan[1]:prob.Δt:tspan[2], 1))


    # 2. Set initial values
    uₐ[:,1] = prob.u₀
    B[:,:,1] = prob.B₀

    # define new version of function with uₐ pre-supplied
    # mdl_forward!(uₖ) = model_forward!(uₖ, uₐ)
    # mdl_forward!(uₖ) = model_forward!(uₖ, uₐ, t_now, t_next, θ)


    times = prob.tspan[1]:prob.Δt:prob.tspan[2]

    n_obs = size(prob.w, 2)
    times_obs = Δt_obs:Δt_obs:(Δt_obs * n_obs)

    for k ∈ 1:(length(times)-1)
        t_now = times[k]
        t_next = times[k+1]

        idx_now = timeindex(t_now, prob.Δt)
        u_now = uₐ[:, idx_now]


        # u_next, DM = model_forward!(u_now)
        u_next, DM = Zygote.withjacobian(model_forward!, u_now, uₐ, t_now, t_next, prob.p)
        DM = DM[1]

        # update error covariance matrix
        B[:,:,k+1] = DM * B[:,:,k] * DM' + prob.Q

    end

    return uₐ, B
end




u0b = [2.0; 3.0; 4.0]  # guess at initial condition
tspan = (0.0, 10.0)
function lorenz!(du, u, p, t)
    x,y,z = u
    σ,ρ,β = p

    du[1] = σ*(y - x)
    du[2] = x*(ρ - z) - y
    du[3] = x*y - β*z
end
θ = [10.0, 28.0, (8/3)]  # true parameter values
Δt = 0.01  # time step for simulation (this may be unnecessary)


Q = 0.0 * I(3)
B₀ = 0.1^2 .* I(3)
R = 0.15^2 .* I(3)
Δt_obs = 0.2

function h(u)
    w = u
    return w
end

function Dh(u)
    n = size(u, 1)
    D = I(n)
    return D
end

w = hcat([rand(MvNormal(zeros(3), R)) for tᵢ ∈ 0.0:Δt_obs:2.0]...)

prob = AssimilationProblem(u0b, tspan, lorenz!, θ, Δt, Q, B₀, R, Δt_obs, h, Dh, w)

sol = EKFsolve(prob)
