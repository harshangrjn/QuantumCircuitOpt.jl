import QuantumCircuitOpt as QCOpt
using LinearAlgebra

using JuMP
using Gurobi
# using CPLEX
# using HiGHS

include("optimizers.jl")
[include("$(i)qubit_gates.jl") for i in 2:5]
include("parametrized_gates.jl")
include("decompose_all_gates.jl")

decompose_gates = ["qubit_routing_circuit"]

#----------------------------------------------#
#      Quantum Circuit Optimization model      #
#----------------------------------------------#
qcopt_optimizer = get_gurobi(solver_log = true)

result = Dict{String,Any}()
times = zeros(length(decompose_gates), 1)

for gates = 1:length(decompose_gates)
    params = getfield(Main, Symbol(decompose_gates[gates]))()

    model_options = Dict{Symbol, Any}(
        :model_type => "compact_formulation",
        :convex_hull_gate_constraints => false,
        :idempotent_gate_constraints  => false,
        :unitary_constraints          => false,
        :tight_unitary_bounds         => false,
    )

    global result = QCOpt.run_QCModel(params, qcopt_optimizer; options = model_options)
    times[gates] = result["solve_time"]
end
