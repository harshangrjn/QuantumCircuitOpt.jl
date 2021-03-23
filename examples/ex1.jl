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
"D" => 8,        # Maximum depth of the decomposition (>= 2)

# Note that, for a given input gate, say H (hadamard), user input should include the gates representations on every qubit, such as H1 and H2. 
# If you prefer to include the kronecker form of gates appearing on adjacent qubits, you can do so by mentioning H1⊗H2

# "elementary_gates" => ["H1", "H2", "H1⊗H2", "Identity", "cnot_12"],  
# "target_gate" => "cnot_21",

# "elementary_gates" => ["H1", "H2", "T1", "T2", "T1_conjugate", "cnot_12", "cnot_21"],
# "target_gate" => "controlled_V",

# "elementary_gates" => ["T1", "Identity", "cnot_12"],  


# "elementary_gates" => ["H1", "H2", "Identity", "cnot_21", "S1", "S2", "cnot_12"],
# "target_gate" => "magic_M",

# "elementary_gates" => ["R_x", "Identity"], 
# "target_gate" => "test_R_x_1",

"elementary_gates" => ["U3", "cnot_12", "cnot_21", "Identity"], 
# "target_gate" => "controlled_Z",
# "target_gate" => "controlled_H_12",

# "elementary_gates" => ["U3", "Identity"],
# "target_gate" => "test_U3_1",

# "elementary_gates" => ["R_x", "R_y", "U3", "cnot_12", "Identity"],
# "elementary_gates" => ["H1", "H2", "Identity", "cnot_21", "cnot_12", "controlled_R2"],
"target_gate" => "qft2",
# "target_gate" => "controlled_R2",
# "target_gate" => "cnot_21",

# Enter discretization angles for each of the matrices which are part of the elementary_gates above.  
# "R_x_discretization" => [π/2, 0], 
# "R_y_discretization" => [-π/4, π/4, π/2, -π/2], 
# "R_z_discretization" => [-π/2, π/2, π/4, -π/4], 
"U_θ_discretization" => [π/2, -π/2],
"U_ϕ_discretization" => [-π/2, π/2],
# "U_λ_discretization" => [-π/4, π/4, π/8, -π/8, π/2, π],
"U_λ_discretization" => [π/2, π, π/8],

"initial_gate" => "Identity", 

# Choose the objective function from one of these options:  
# "minimize_depth": Minimizes the total depth of decomposition. For this option, include "Identity" matrix in the "elementary_gates" above. 
# "minimize_cnot" : Minimizes the number of CNOT gates in the decomposition. 
"objective" => "minimize_depth", 

# Choose the type of decomposition from one of these options:
# "exact": QuantumCircuitOpt finds an exact decomposition if it exists
# "approximate": QuantumCircuitOpt finds an approximate decomposition if an exact one does not exist; otherwise it will return an exact solution.
"decomposition_type" => "exact",

# Choose the mixed-integer optimizer here. Note that CPLEX or Gurobi will be ideal for faster run times.
"optimizer" => "cplex",
"presolve" => true,
"optimizer_log" => true, 
"relax_integrality" => false,
                            
)

#------------------------------#
#      Optimization model      #
#------------------------------#
qcm_optimizer = get_solver(params)
data = QuantumCircuitOpt.get_data(params)

model_qc  = QuantumCircuitOpt.build_QCModel(data, model_type = "compact_formulation", commute_matrix_cuts = false)
result_qc = QuantumCircuitOpt.optimize_QCModel!(model_qc, optimizer = qcm_optimizer)
QuantumCircuitOpt.visualize_QCModel_solution(result_qc, data)
