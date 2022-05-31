struct AssimilationProblem
    """The solution initial condition, i.e. `u(tspan[1])=u₀` """
    u₀
    """The solution `u(t)` will be computed for `tspan[1] ≤ t ≤ tspan[2]`. """
    tspan
    """The underlying ODE is `du = f(u,p,t)` for out-of-place and `f(du,u,p,t)` for in-place."""
    f
    """Constant parameters to be supplied as the second argument of f"""
    p
    """Analysis time step """
    Δt
    """Process noise covariance matrix"""
    Q
    """Initial model error covariance matrix"""
    B₀
    """Observation error covariance matrix"""
    R
    """Observation time step """
    Δt_obs
    """Observation function `h:u\to w` """
    h
    """Jacobian of the observation function, i.e. `Dh = ∂ h_i / ∂ u_j` """
    Dh
    """Observation data"""
    w
end

export AssimilationProblem
