module Filtration
greet() = print("Hello World!")

# using
using DocStringExtensions
using LaTeXStrings
using Plots
using ForwardDiff
using Symbolics, Nemo

# include 
include("models/membrane_filtration_model.jl")

# export functions only for user
export membrane_filtration_model
export ismonotonic, isLfunction, isKfunction
export get_roots

"""
[`Filtration`](@ref) module.

Lists all the imported modules and packages:

$(IMPORTS)

List of all the exported names:

$(EXPORTS)
"""

end