import QuantumCircuitOpt as QCO
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
                   "decompose_iSwapGate",
                   "decompose_qft2_using_HT",
                   "decompose_RX_on_q3"]

# decompose_gates = ["decompose_GroverDiffusion_using_HX"]

#----------------------------------------------#
#      Quantum Circuit Optimization model      #
#----------------------------------------------#
result_qc = Dict{String,Any}()

for gates in decompose_gates 
    
    params = getfield(Main, Symbol(gates))()

    model_options = Dict{Symbol, Any}(:model_type => "compact_formulation",
                                      :convex_hull_gate_constraints => false)
    
    qcopt_optimizer = get_gurobi()
    global result_qc = QCO.run_QCModel(params, qcopt_optimizer)
end