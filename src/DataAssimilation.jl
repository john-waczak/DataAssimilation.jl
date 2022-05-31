module DataAssimilation

# Write your package code here.
include("AssimilationProblem.jl")
include("EKF.jl")

export AssimilationProblem
export timeindex
export update!
export model_forward!
export EKFsolve


end
