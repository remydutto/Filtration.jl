module Filtration

greet() = print("Hello World!")

# using
using LaTeXStrings
using Plots
using ForwardDiff
using Symbolics, Nemo

# include 

include("models/Benyahia_and_al.jl")

# export functions only for user
export membrane_filtration_model
export ismonotonic, isLfunction, isKfunction
export get_roots

end