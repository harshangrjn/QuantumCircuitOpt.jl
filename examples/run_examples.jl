using QuantumCircuitOpt
using JuMP
using CPLEX
using LinearAlgebra
using Cbc

const QCO = QuantumCircuitOpt

include("solver.jl")
include("2qubit_gates.jl")


gate_params = [test_hadamard(), 
               test_S(), 
               test_magic_M(), 
               test_controlled_R2(), 
               test_cnot_21(), 
               test_cnot_21_with_U(), 
               test_controlled_Z(),
               test_magic_M_using_SHCnot(),
               test_swap(),
               test_controlled_V(),
               test_W(),
               test_W_using_HCnot()]

# gate_params = [test_cnot_21()]

#------------------------------#
#      Optimization model      #
#------------------------------#
function run_QCO(params)
    
    qcm_optimizer = get_solver(params)  
    result_qc = QCO.run_QCModel(params, qcm_optimizer)

    return result_qc
end

result_qcm = Dict{}
for params in gate_params 
    global result_qcm = run_QCO(params)
end

