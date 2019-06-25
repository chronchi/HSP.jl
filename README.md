HSP - Hansen Solubility Parameter
=========

Package that provides functions to calculate the
[HSPs](https://en.wikipedia.org/wiki/Hansen_solubility_parameter) given a file
containing a list of vectors in Hansen Space.

# Installation
In order to install the package follow the steps below

```julia
using Pkg
Pkg.add("https://github.com/chronchi/HSP.jl")
```

# Usage

To calculate the best Hansen Solubility Parameters you just need to
give as input to the function a path to a file containing the parameters.
The output will be the RED and a structure containing the parameters and R_o.

```julia
path_to_parameters = "path/to/file"
RED, MinimumHansenP = minimumRED(path_to_parameters)
```

To access the minimum values, like R_o and δ's, do as following

```julia
δ_d = MinimumHansenP.delta_d
δ_p = MinimumHansenP.delta_p
δ_h = MinimumHansenP.delta_h
R_o = MinimumHansenP.R_o
```

If the file has a header, specify when calling the function.
You can also specify the delimiter of the file.
```julia
header = 1
delim = ','
path_to_parameters = "path/to/file"
RED, MinimumHansenP = minimumRED(path_to_parameters,
                        header=header, delim=delim)
```

If you already have the data loaded in a matrix where the
first column is δ_d, second δ_p, third δ_h, fourth the solubility, then you can pass it as input specifying by *data*
as follows

```julia
data # pre loaded data
RED, MinimumHansenP = minimumRED(data=data)
```


# Arguments

Arguments for the function `minimumRED`.

- `path_to_parameters::String`: path to parameters file.
- `random_guess::Bool`: if the initial guess is random or not
- `tol::Float64`: tolerance value for the optimization
- `max_iter::Integer`: maximum number of iterations to find the minimum
-  `optimizer::UnionAll`: optimizer used in the minimization problem

One can also specify the parameters for the optimizer as an
argument in the function `minimumRED`
