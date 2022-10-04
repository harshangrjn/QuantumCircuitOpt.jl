import QuantumCircuitOpt as QCOpt
using LinearAlgebra

using JuMP
using Gurobi
# using Cbc
# using CPLEX

include("optimizers.jl")
include("2qubit_gates.jl")
include("3qubit_gates.jl")
include("4qubit_gates.jl")
include("5qubit_gates.jl")

decompose_gates = ["decompose_hadamard",
                   "decompose_controlled_Z",
                   "decompose_controlled_V",
                   "decompose_controlled_H",
                   "decompose_controlled_H_with_R",
                   "decompose_magic", 
                   "decompose_magic_using_SHCnot",
                   "decompose_S", 
                   "decompose_revcnot", 
                   "decompose_revcnot_with_U", 
                   "decompose_swap",
                   "decompose_W",
                   "decompose_W_using_HCnot",
                   "decompose_GroverDiffusion_using_Clifford",
                   "decompose_GroverDiffusion_using_U3",
                   "decompose_iSwap",
                   "decompose_qft2_using_HT",
                   "decompose_RX_on_q3"]

decompose_gates = ["decompose_CNot_41"]

#----------------------------------------------#
#      Quantum Circuit Optimization model      #
#----------------------------------------------#
result = Dict{String,Any}()

for gates in decompose_gates 
    
    params = getfield(Main, Symbol(gates))()

    model_options = Dict{Symbol, Any}(:convex_hull_gate_constraints => false,
                                      :unit_magnitude_constraints   => false,
                                      :idempotent_gate_constraints  => false,
                                      :unitary_constraints          => false)
    
    qcopt_optimizer = get_gurobi()
    global result = QCOpt.run_QCModel(params, qcopt_optimizer; options = model_options)
end
