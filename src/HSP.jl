module HSP

import Optim
using Statistics: mean
using DataFrames
import CSV
using LinearAlgebra: dot
using DelimitedFiles

export minimumRED, load_parameters

# define a struct to save the loaded parameters
mutable struct HansenParameters
    delta_d::Vector{<:Real}
    delta_p::Vector{<:Real}
    delta_h::Vector{<:Real}
    solubility::Vector{<:Real}
end

mutable struct MinimumHansenParameters
    delta_d::Real
    delta_p::Real
    delta_h::Real
    R_o::Real
end

"""
    load_parameters(path_to_parameters::String; kwargs...)

Load the table with parameters for the current problem.
"""
function load_parameters(path_to_parameters::String; kwargs...)
    # defining parameters
    header = false
    delim  = ','
    for (p, v) in kwargs
        if p == :header
            header = v
        elseif p == :delim
            delim = v
        end
    end
    # load the file with the parameters as a DataFrame
    df = CSV.read(path_to_parameters, header=header, delim=delim; kwargs...)
    # add the values to HansenParameters
    HansenParameters(convert(Array{Float64, 1}, df[:,1]),
                     convert(Array{Float64, 1}, df[:,2]),
                     convert(Array{Float64, 1}, df[:,3]),
                     convert(Array{Float64, 1}, df[:,4]))
end

"""
    QFval(targest_deltaR::Vector{Float64}, HansenP::HansenParameters)

Description missing
"""
function QFval(targets_deltaR::Vector{Float64}, HansenP::HansenParameters)
    delta_dp = targets_deltaR[1]
    delta_pp = targets_deltaR[2]
    delta_hp = targets_deltaR[3]
    R_o  = targets_deltaR[4]

    # define the Hansen Parameters
    delta_d = HansenP.delta_d
    delta_p = HansenP.delta_p
    delta_h = HansenP.delta_h
    solubility = HansenP.solubility

    # Calculate the distance between Hansen Parameters in Hansen Space
    R_a = sqrt.(
    4*(delta_d .- delta_dp).^2 + (delta_p .- delta_pp).^2 + (delta_h .- delta_hp).^2
    )

    # define the number of solvents
    solvent_number = size(HansenP.delta_p, 1)

    # initialize vector to store values depending on the solubility
    A = zeros(solvent_number)
    for i = 1:solvent_number
        if R_a[i] > R_o
            if solubility[i] == 0
                A[i] = 1
            else
                A[i] = exp(R_o - R_a[i])
            end
        else
            if solubility[i] == 0
                A[i] = exp(R_a[i] - R_o)
            else
                A[i] = 1
            end
        end
    end
    # taking the n-th root of each entry and multiplying them
    DataFit = prod(A)^(1/solvent_number)
    # return the absolute value of Datafit - 1
    return abs(DataFit - 1)
end

"""
    minimumRED(path_to_parameters::String, tol::Real=1e-2,
                 optimizer::UnionAll=Optim.NelderMead,
                 max_iter::Int=1000,
                 random_guess::Bool=false[; kwargs...])

Find minimum values of the hansen solubility parameters and RED.

# Examples
```jldoctest
# path to file
path_to_parameters = "path/to/file.csv"
RED, MinimumHansenP = minimumRED(path_to_parameters);
```
"""

function minimumRED(path_to_parameters::String=""; tol::Real=1e-2,
                       optimizer=Optim.NelderMead, max_iter::Int=1000,
                       random_guess=false, kwargs...)
    verbosity = 0
    pre_loaded = nothing
    for (p,v) in kwargs
        if p == :verbosity
            verbosity = v
        elseif p == :data
            pre_loaded = v
        end
    end
    # if the user already gives the preloaded data
    if pre_loaded != nothing
        HansenP = HansenParameters(convert(Array{Float64,1}, pre_loaded[:,1]),
                    convert(Array{Float64,1}, pre_loaded[:,2]),
                    convert(Array{Float64,1}, pre_loaded[:,3]),
                    convert(Array{Float64,1}, pre_loaded[:,4]))
    else
        # load the given hasen parameters
        HansenP = load_parameters(path_to_parameters; kwargs...)
    end
    if verbosity == 1
        println("Parameters loaded.")
    end
    # define the initial guess
    if random_guess
        initial_guess = rand(4)
    else
        # define the mean for each parameter vector
        mean_delta = [mean(getfield(HansenP, i))
                   for i in fieldnames(HansenParameters)[1:3]]
        # append a 4th entry, it depends on the previous 3 entries
        push!(mean_delta, sqrt(dot(mean_delta,mean_delta)))
        initial_guess = mean_delta
    end
    # define a function that will be minized
    QF(x) = QFval(x, HansenP)
    # start the minimization
    counter = 0
    min_val = 1
    result_optim = 0
    n_interval = floor(max_iter/5)
    while min_val > tol
        # minimize the function and output all the results
        result_optim = Optim.optimize(QF, initial_guess, optimizer(;kwargs...))
        # the value attained by the function QF in its minimum
        min_val = Optim.minimum(result_optim)
        initial_guess = Optim.minimizer(result_optim)
        # check the counter condition
        if counter > max_iter
            break
        else
            if counter % n_interval == 0
                println("Iteration: " * string(counter) * ", Tol: " * string(min_val))
            end
            counter = counter + 1
        end
    end
    # the minimum of QF
    QF_minimum = Optim.minimizer(result_optim)
    MinimumHansenP = MinimumHansenParameters(QF_minimum[1], QF_minimum[2],
                                             QF_minimum[3], QF_minimum[4])
    # the difference
    diff_params = [getfield(HansenP, i) .- getfield(MinimumHansenP, i)
                    for i in fieldnames(HansenParameters)[1:3]]
    # calculate Ra
    R_a = sqrt.(4*diff_params[1].^2 + diff_params[2].^2 + diff_params[3].^2)
    # calculate RED
    RED_values = R_a ./ MinimumHansenP.R_o
    print("\n")
    println("===================")
    print("R_o value: ")
    println(MinimumHansenP.R_o)
    println("===================")
    println("Other values")
    df = DataFrame(RED = RED_values, δ_d = HansenP.delta_d, δ_p = HansenP.delta_p,
                    δ_h = HansenP.delta_h)
    println(df)
    return RED_values, MinimumHansenP
end

end # module
