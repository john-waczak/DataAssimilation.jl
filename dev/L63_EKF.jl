# using DifferentialEquations
using LinearAlgebra
using Statistics
using Distributions
using Plots
using Random

function Lorenz63(u, θ)
    σ,β,ρ = θ
    x,y,z = u

    du = [
        σ * (y-x);
        x*(ρ-z) - y;
        x*y - β*z
    ]
    return du
end



function JLorenz63(u, θ)
    σ,β,ρ = θ
    x,y,z = u

    J = [
        -σ σ 0;
        (ρ-z) -1 -x;
        y x -β
    ]

    return J
end



function RK4(f, u, dt, θ)
    k1 = f(u,θ)
    k2 = f(u+k1*dt/2, θ)
    k3 = f(u+k2*dt/2, θ)
    k4 = f(u+k3*dt, θ)

    u_new = u + (dt/6)*(k1 + 2k2 + 2k3 + k4)
    return u_new
end


function JRK4(f, Jf, u, dt, θ)
    n = size(u,1)
    k1 = f(u, θ)
    k2 = f(u+k1*dt/2, θ)
    k3 = f(u+k2*dt/2, θ)
    #k4 = rhs(state+k3*dt,*args)

    dk1 = Jf(u, θ)
    dk2 = Jf(u+k1*dt/2, θ) * (I(n) + dk1*dt/2)
    dk3 = Jf(u+k2*dt/2, θ) * (I(n) + dk2*dt/2)
    dk4 = Jf(u+k3*dt, θ) * (I(n) + dk3*dt)

    return I(n) + (dt/6)*(dk1 + 2dk2 + 2dk3 + dk4)
end




function EKF(ub, w, ObsOp, JObsOp, R, B)
    n = size(ub, 1)

    Dh = JObsOp(ub)

    D = Dh*B*Dh' + R
    K = B * Dh' * inv(D)

    ua = ub + K * (w - ObsOp(ub))
    P = (I(n) - K*Dh)*B

    return ua, P
end


function h(u)
    w = u
    return w
end


function Dh(u)
    n = size(u, 1)
    D = I(n)
    return D
end



θ = [10.0, (8/3), 28.0]
dt = 0.01
tm = 10
nt = Int(tm/dt)
t = 0:dt:tm


u0True = [1.0; 1.0; 1.0]
Random.seed!(1)
σ_m = 0.15
R = σ_m^2 .* I(3)

dt_m = 0.2
tm_m = 2
nt_m = Int(tm_m/dt_m)

ind_m = (Int(dt_m/dt)+1):Int(dt_m/dt):(Int(tm_m/dt)+1)
t_m = t[ind_m]

uTrue = zeros(3,nt+1)
uTrue[:,1] = u0True
km = 1
w = zeros(3, nt_m)
for k ∈ 1:nt
    uTrue[:, k+1] = RK4(Lorenz63, uTrue[:,k], dt, θ)

    if (km <= nt_m) && (k+1==ind_m[km])
        w[:,km] = h(uTrue[:,k+1]) + rand(Normal(0,σ_m), 3)
        km = km+1
    end
end

u0b = [2.0; 3.0; 4.0]
σ_b = 0.1
B = σ_b^2 .* I(3)
Q = 0.0 * I(3)

ub = zeros(3, nt+1)
ub[:,1] = u0b

ua = zeros(3, nt+1)
ua[:,1] = u0b
km = 1

for k ∈ 1:nt
    ub[:, k+1] = RK4(Lorenz63, ub[:,k], dt, θ)
    ua[:, k+1] = RK4(Lorenz63, ua[:,k], dt, θ)

    DM = JRK4(Lorenz63, JLorenz63, ua[:,k], dt, θ)

    B = DM * B * DM' +  Q

    if (km <= nt_m) && (k+1 == ind_m[km])
        ua[:, k+1], B = EKF(ua[:, k+1], w[:, km], h, Dh, R, B)
        km = km+1
    end

end

var_names = ["x(t)", "y(t)", "z(t)"]
plots = []
for i ∈ 1:3
    p = plot(t, uTrue[i, :], color=:black, alpha=0.75, linewidth=2, label="truth")
    plot!(p, t_m, w[i, :], seriestype=:scatter, color="light green", label="observation" )
    plot!(p, t, ub[i, :], linestyle=:dashdot, color="purple", label="forecast")
    plot!(p, t, ua[i, :], color="royal blue", alpha=1, linewidth = 1, label="analysis")
    xlabel!(p, "t")
    ylabel!(p, var_names[i])
    push!(plots, p)
end

plot(plots..., layout=(3,1))
