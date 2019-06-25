using HSP
using Test

RED, MinimumHansenP = minimumRED("parameter_values.csv", header=true, tol=1e-2,
                                max_iter=100)

@test typeof(RED) == Array{Float64, 1}
