module NMTK
using Reexport
using LinearAlgebra
@reexport using ModelingToolkit
@reexport using DiffEqOperators
@reexport using DifferentialEquations
@reexport using Plots

include("diffeqopfix.jl")

export discretize
end
