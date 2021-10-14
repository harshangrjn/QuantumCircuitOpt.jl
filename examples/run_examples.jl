import QuantumCircuitOpt as QCO
using JuMP
using CPLEX
using Gurobi
using Cbc

include("solver.jl")
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
                   "decompose_GroverDiffusionGate",
                   "decompose_iSwapGate",
                   "decompose_RX_on_q3"]

# decompose_gates = ["decompose_W"]

#----------------------------------------------#
#      Quantum Circuit Optimization model      #
#----------------------------------------------#
result_qc = Dict{String,Any}()

for gates in decompose_gates 
    params = getfield(Main, Symbol(gates))()
    
    qcm_optimizer = get_solver(params)
    # qcm_optimizer = JuMP.optimizer_with_attributes(Cbc.Optimizer)
    # qcm_optimizer = JuMP.optimizer_with_attributes(CPLEX.Optimizer)
    # qcm_optimizer = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "Presolve" => 1)
    
    global result_qc = QCO.run_QCModel(params, 
                                       qcm_optimizer, 
                                       model_type = "compact_formulation",
                                       convex_hull_complex_gate_constraints = false)
end

