using Documenter
using QuantumCircuitOpt

makedocs(
    sitename = "QuantumCircuitOpt",
    format = Documenter.HTML(),
    modules = [QuantumCircuitOpt]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
