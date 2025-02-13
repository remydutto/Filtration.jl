module Filtration

# using
using LaTeXStrings
using Plots
using ForwardDiff
using Symbolics, Nemo

# include 

include("models/Benyahia_and_al.jl")

# export functions only for user
export Benyahia
export ismonotonic, isLfunction, isKfunction
export get_roots

end