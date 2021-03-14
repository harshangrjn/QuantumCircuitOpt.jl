using QuantumCircuitOpt
using JuMP
using CPLEX
using LinearAlgebra
#using Cbc

include("solver.jl")

# User-defined inputs
params = Dict{String, Any}(
"n_qubits" => 2, # Number of qubits
"D" => 3, # Maximum depth of the decomposition (>= 2)

# Note that, for a given input gate, say H (hadamard), user input should include the gates representations on every qubit, such as H1 and H2. 
# If you prefer to include the kronecker form of gates appearing on adjacent qubits, you can do so by mentioning H⊗H
"elementary_gates" => ["H1", "H2", "H⊗H", "cnot_12"], 
"target_gate" => "cnot_21",
"initial_gate" => "Identity",

# If you prefer to use Universal gates as inputs (like in IBM architecture), provide the discretization angles here
"U_gate_discretizations" => [], 

"optimizer" => "cplex",
"presolve" => true,
"optimizer_log" => true,                           
"binary_relax" => false,
                            
# Valid inequalities which may speed up the model run time
"cuts_1" => false, #commutative matrices
)

qcm_optimizer = get_solver(params)
data = QuantumCircuitOpt.get_data(params)

model_qc = QuantumCircuitOpt.build_QCModel(data)
results = QuantumCircuitOpt.run_QCModel(model_qc, optimizer = qcm_optimizer)


