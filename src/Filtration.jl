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
using OptimalControl
using OrdinaryDiffEq
using MINPACK

# useless function
greet() = print("Hello World!")

# include 
include("models/membrane_filtration_model.jl")

# export functions only for user
export MembraneFiltrationModel
export ismonotonic, isLfunction, isKfunction
export singular_state, singular_control, singular_costate
export get_psi, get_phi, get_dphi
export delta_t_end

end