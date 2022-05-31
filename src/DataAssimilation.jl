module DataAssimilation

# Write your package code here.
include("AssimilationProblem.jl")
include("EKF.jl")
include("4dvar.jl")

export AssimilationProblem
export timeindex
export update!
export model_forward!
export EKFsolve
export FourDVarSolve
end
