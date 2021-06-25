using QuantumCircuitOpt
using JuMP
using CPLEX
using LinearAlgebra
using Cbc

const QCO = QuantumCircuitOpt

include("solver.jl")
include("2qubit_gates.jl")


test_gates = ["test_hadamard", 
               "test_S", 
               "test_magic_M",  
               "test_cnot_21", 
               "test_cnot_21_with_U", 
               "test_controlled_Z",
               "test_magic_M_using_SHCnot",
               "test_swap",
               "test_controlled_V",
               "test_W",
               "test_W_using_HCnot"]

#----------------------------------------------#
#      Quantum Circuit Optimization model      #
#----------------------------------------------#
for gates in test_gates 
    params = getfield(Main, Symbol(gates))()
    qcm_optimizer = get_solver(params)  
    result_qcm = QCO.run_QCModel(params, qcm_optimizer)
end

