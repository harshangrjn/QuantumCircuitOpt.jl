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

# "elementary_gates" => ["H1", "H2", "H⊗H", "Identity", "cnot_21"],  
# "target_gate" => "controlled_Z",

# "elementary_gates" => ["R_x", "R_z", "cnot_12", "Identity"], 
"elementary_gates" => ["R_x", "Identity"], 
"target_gate" => "X1",

"R_x_discretization" => [π/2],
"R_y_discretization" => [-π/4, π/4],
"R_z_discretization" => [π/4, π/2, π, 5*π/4, 3*π/2, 2*π],
"U_θ_discretization" => [],
"U_ϕ_discretization" => [],
"U_λ_discretization" => [],

"initial_gate" => "Identity",

# Choose the objective function, which is either to minimize the number of CNOT gates or 
# the total depth of decomposition. Specify the exact type of cnot gate which needs to be minimized.
# If depth is to be minimized, use "depth" as the option.
"objective" => "depth",

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


