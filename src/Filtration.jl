"""
[`Filtration`](@ref) module.

Lists all the imported modules and packages:

$(IMPORTS)

List of all the exported names:

$(EXPORTS)
"""
module Filtration

# using
using DocStringExtensions
using LaTeXStrings
using Plots
using ForwardDiff
using Symbolics, Nemo
using Roots

# useless function
greet() = print("Hello World!")

# include 
include("models/membrane_filtration_model.jl")

# export functions only for user
export MembraneFiltrationModel
export ismonotonic, isLfunction, isKfunction
export get_root, get_roots_symbolic_algebraic_fraction
export get_psi

end