using Flux
using DiffEqSensitivity
using DifferentialEquations
using LinearAlgebra
using Statistics
using Distributions
using Random

include("./AssimilationProblem.jl")









function FourDVarSolve(prob::AssimilationProblem, nsteps::Int, opt)

    nobs = size(prob.w, 2)
    tobs_end = (nobs-1)*prob.Δt_obs
    Rinv = collect(inv(prob.R))

    function loss()
        res = Array(solve(ODEProblem(prob.f, prob.u₀, prob.tspan, prob.p), saveat=prob.Δt_obs:prob.Δt_obs:tobs_end))
        # ℓ = 0.5*(u0a-u0b)'*Binv*(u0a-u0b) + 0.5*sum((w[:,i]-h(res[:,i]))'*Rinv*(w[:,i]-h(res[:,i])) for i ∈ 1:size(res,2))

        ℓ = 0.5*sum((prob.w[:,i]-prob.h(res[:,i]))'*Rinv*(prob.w[:,i]-prob.h(res[:,i])) for i ∈ 1:size(res,2))
        return ℓ
    end


    data = Iterators.repeated((), nsteps)

    cb = function ()
        display(loss())
    end

    u0a = prob.u₀

    Flux.train!(loss, Flux.params(u0a), data, opt, cb=cb)

    return u0a
end

