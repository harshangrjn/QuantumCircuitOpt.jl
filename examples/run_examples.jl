import QuantumCircuitOpt as QCO
using JuMP
using CPLEX
using LinearAlgebra
using Cbc

include("solver.jl")
include("2qubit_gates.jl")
include("3qubit_gates.jl")
include("4qubit_gates.jl")
include("5qubit_gates.jl")

decompose_gates = ["decompose_hadamard", 
                    "decompose_S", 
                    "decompose_magic_M",  
                    "decompose_revcnot", 
                    "decompose_revcnot_with_U", 
                    "decompose_controlled_Z",
                    "decompose_magic_M_using_SHCnot",
                    "decompose_swap",
                    "decompose_controlled_V",
                    "decompose_W",
                    "decompose_W_using_HCnot",
                    "decompose_GroverDiffusionGate",
                    "decompose_RX_on_q3"]

decompose_gates = ["decompose_S"]

#----------------------------------------------#
#      Quantum Circuit Optimization model      #
#----------------------------------------------#
result_qc = Dict{String,Any}()

for gates in decompose_gates 
    params = getfield(Main, Symbol(gates))()
    qcm_optimizer = get_solver(params)
    global result_qc = QCO.run_QCModel(params, 
                                       qcm_optimizer, 
                                       model_type = "compact_formulation")
end

