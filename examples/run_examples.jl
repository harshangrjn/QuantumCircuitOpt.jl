import QuantumCircuitOpt as QCOpt
using LinearAlgebra

using JuMP
using Gurobi
# using HiGHS

include("optimizers.jl")
include("2qubit_gates.jl")
include("3qubit_gates.jl")
include("4qubit_gates.jl")
include("5qubit_gates.jl")
include("parametrized_gates.jl")
include("decompose_all_gates.jl")

decompose_gates = ["parametrized_hermitian_gates"]

#----------------------------------------------#
#      Quantum Circuit Optimization model      #
#----------------------------------------------#
qcopt_optimizer = get_gurobi()

result = Dict{String,Any}()
times = zeros(length(decompose_gates), 1)

for gates = 1:length(decompose_gates)
    params = getfield(Main, Symbol(decompose_gates[gates]))()

    model_options = Dict{Symbol, Any}(
        :model_type => "compact_formulation",
        :convex_hull_gate_constraints => false,
        :idempotent_gate_constraints  => false,
        :unitary_constraints          => false,
        :fix_unitary_variables        => true,
    )

    global result = QCOpt.run_QCModel(params, qcopt_optimizer; options = model_options)
    times[gates] = result["solve_time"]
end
