module Filtration

# using
using DocStringExtensions
using LaTeXStrings
using Plots
using ForwardDiff
using Symbolics, Nemo

# useless function
greet() = print("Hello World!")

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