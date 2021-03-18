using QuantumCircuitOpt
using JuMP
using CPLEX
using LinearAlgebra
using Cbc

include("solver.jl")

#-------------------------------#
#      User-defined inputs      #
#-------------------------------#
params = Dict{String, Any}(
"n_qubits" => 2, # Number of qubits
"D" => 5, # Maximum depth of the decomposition (>= 2)

# Note that, for a given input gate, say H (hadamard), user input should include the gates representations on every qubit, such as H1 and H2. 
# If you prefer to include the kronecker form of gates appearing on adjacent qubits, you can do so by mentioning H⊗H

# "elementary_gates" => ["H1", "H2", "H⊗H", "Identity", "cnot_12"],  
# "elementary_gates" => ["T1", "Identity", "cnot_12"],  
# "target_gate" => "S1",
# "target_gate" => "cnot_21",

# "elementary_gates" => ["R_y", "R_z", "cnot_12", "Identity"], 
"elementary_gates" => ["R_y", "cnot_12", "Identity"], 
# "target_gate" => "test_R_x_1",
"target_gate" => "controlled_Z",

# Enter discretization angles for each of the matrices which are part of the elementary_gates above. 
"R_x_discretization" => [-π/2, -π/4, π/4, π/2],
"R_y_discretization" => [-π/2, -π/4, π/4, π/2],
"R_z_discretization" => [-π/2, -π/4, π/4, π/2],
"U_θ_discretization" => [],
"U_ϕ_discretization" => [],
"U_λ_discretization" => [],

"initial_gate" => "Identity",

# Choose the objective function from one of these options:  
# "minimize_depth": Minimizes the total depth of decomposition. For this option, include "Identity" matrix in the "elementary_gates" above. 
# "minimize_cnot" : Minimizes the number of CNOT gates in the decomposition. 
"objective" => "minimize_depth", 

# Choose the type of decomposition from one of these options:
# "exact": QuantumCircuitOpt finds an exact decomposition if it exists
# "approximate": QuantumCircuitOpt finds an approximate decomposition if an exact one does not exist; otherwise it will return an exact solution.
"decomposition_type" => "exact",

"optimizer" => "cplex",
"presolve" => true,
"optimizer_log" => true,                           
"relax_integrality" => false,
                            
# Valid inequalities which may speed up the model run time
"cuts_1" => false, #commutative matrices
)

#------------------------------#
#      Optimization model      #
#------------------------------#
qcm_optimizer = get_solver(params)
data = QuantumCircuitOpt.get_data(params)

model_qc = QuantumCircuitOpt.build_QCModel(data)
result_qc = QuantumCircuitOpt.optimize_QCModel!(model_qc, optimizer = qcm_optimizer)
QuantumCircuitOpt.visualize_QCModel_solution(result_qc, data)


