using HSP, CSV
using Test

# test loading data with load_parameters function
RED, MinimumHansenP = minimumRED("parameter_values.csv", header=1, tol=1e-2,
                                max_iter=100)

@test typeof(RED) == Array{Float64, 1}

# test with pre loaded data
df = CSV.read("parameter_values.csv", header=1)
RED, MinimumHansenP = minimumRED(data=df)

@test typeof(RED) == Array{Float64, 1}
