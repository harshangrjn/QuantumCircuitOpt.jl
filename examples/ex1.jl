using QuantumCircuitOpt
using JuMP
using CPLEX
#using Cbc

include("solver.jl")

# All the user 
params = Dict{String, Any}(
"n_qubits" => 2, # Number of qubits
"D" => 5, # Maximum depth of the decomposition (>= 2)

# Note that, for a given input gate, say H (hadamard), user input should include the gates representations on every qubit, such as H1 and H2. 
# If you prefer to include the kronecker form of gates appearing on adjacent qubits, you can do so by mentioning H⊗H
"elementary_gates" => ["H1", "H2", "H⊗H", "cnot_12"], 
"target_gate" => "cnot_21",

# If you prefer to use Universal gates as inputs (like in IBM architecture), provide the discretization angles here
"U_gate_discretizations" => [], 

"optimizer" => "cplex",
"presolve" => true,
"optimizer_log" => true,                           
"lp_relax" => false,
                            
# Valid inequalities
"cuts_1" => false, #commutative matrices
)

qcm_optimizer = get_solver(params)
data = QuantumCircuitOpt.get_data(params)

# model_qc = build_QCModel(data,params)
# results = run_QCModel(model_qc, optimizer = cplex)


