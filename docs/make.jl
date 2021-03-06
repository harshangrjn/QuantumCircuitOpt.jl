using Documenter
using QuantumCircuitOptimization

makedocs(
    sitename = "QuantumCircuitOptimization",
    format = Documenter.HTML(),
    modules = [QuantumCircuitOptimization]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
