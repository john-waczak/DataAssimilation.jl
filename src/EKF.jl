using DifferentialEquations

using LinearAlgebra

using DiffEqSensitivity
using Zygote

using Statistics
using Distributions


include("./AssimilationProblem.jl")



"""
    timeindex(t, Δt)

For a given time value and time step, return the corresponding array index.
"""
function timeindex(t, Δt)
    return Int(round(t/Δt)) + 1
end



"""
    update!(sol, t_now, t_next, Δt, uₐ)

Update solution array uₐ with the values from the solution object sol at specified time step Δt
"""
function update!(sol, t_now, t_next, Δt, uₐ)
    idx_now = timeindex(t_now, Δt)
    idx_next = timeindex(t_next, Δt)

    uₐ[:, idx_now:idx_next] = hcat(sol(t_now:Δt:t_next)...)
end
Zygote.@nograd update!  # tell zygote to ignore update! when computing derivatives



"""
    model_forward!(uₖ, uₐ, t_now, t_next, prob)

Given value uₖ and time range t_now and t_next, integrate the ODE problem forward and update solution object.
"""
function model_forward!(uₖ, uₐ, t_now, t_next, prob)
    ode_prob = ODEProblem(prob.f, uₖ, (t_now, t_next), prob.p)
    sol = solve(ode_prob, Tsit5(), reltol=1e-6,abstol=1e-6)

    update!(sol, t_now, t_next, prob.Δt, uₐ)

    return sol[:, end]
end



"""
    EKFsolve(prob::AssimilationProblem)

Given an AssimilationProblem `prob`, obtain the analysis uₐ and error covariance matrix B using the Extended Kalman Filter (EKF) method.
"""
function EKFsolve(prob::AssimilationProblem)
    # 1. allocate solution arrays for u(t) and B(t)
    uₐ = zeros(size(prob.u₀, 1), size(prob.tspan[1]:prob.Δt:prob.tspan[2], 1))
    B = zeros(size(prob.u₀,1), size(prob.u₀,1), size(prob.tspan[1]:prob.Δt:prob.tspan[2], 1))


    # 2. Set initial values
    uₐ[:,1] = prob.u₀
    B[:,:,1] = prob.B₀


    # 3. define time ranges
    times = prob.tspan[1]:prob.Δt:prob.tspan[2]

    n_obs = size(prob.w, 2)
    times_obs = prob.Δt_obs:prob.Δt_obs:(prob.Δt_obs * (n_obs-1))

    for k ∈ 1:(length(times)-1)

        # 4. Perform the forecast to obtain the background

        t_now = times[k]
        t_next = times[k+1]

        idx_now = timeindex(t_now, prob.Δt)
        idx_next = timeindex(t_next, prob.Δt)

        u_now = uₐ[:, idx_now]

        u_next, DM = Zygote.withjacobian(model_forward!, u_now, uₐ, t_now, t_next, prob)
        DM = DM[1]

        # update error covariance matrix
        B[:,:,idx_next] = DM * B[:,:,idx_now] * DM' + prob.Q


        # 5. Perform the assimilation step
        if t_next ∈ times_obs
            n = size(uₐ, 1)
            idx_h = timeindex(t_next, prob.Δt_obs) - 1 # -1 since this doesn't include 0

            # compute the Kalman gain matrix
            DH = prob.Dh(uₐ[:,idx_next])

            denom = DH*B[:,:, idx_next]*DH' + prob.R

            K = B[:,:, idx_next]*DH'*inv(denom)

            uₐ[:,idx_next] .= uₐ[:, idx_next] + K*(prob.w[:, idx_h] - prob.h(uₐ[:, idx_next]))

            B[:,:,idx_next] = (I(n) - K*DH)*B[:,:,idx_next]
        end


    end


    # 6. Return the resulting analysis uₐ, and error covariance matrix B
    return uₐ, B
end

