using QuantumCircuitOpt
using JuMP
using CPLEX
using LinearAlgebra

include("solver.jl")

#-------------------------------#
#      User-defined inputs      #
#-------------------------------#
params = Dict{String, Any}(
"n_qubits" => 2, 
"depth" => 5,  

# "elementary_gates" => ["R_x", "Identity"], 
# "target_gate" => "test_R_x_1",

# "elementary_gates" => ["U3", "Identity"],
# "target_gate" => "test_U3_1",

"R_x_discretization" => [π/4], 
"R_y_discretization" => [-π/4, π/4, π/2, -π/2, -π], 
"R_z_discretization" => [-π/2, π/2, π/4, -π/4, -π], 
"U_θ_discretization" => [-π/2, π/2],
"U_ϕ_discretization" => [0, π/2],
"U_λ_discretization" => [0, π/2],

"initial_gate" => "Identity", 
"objective" => "minimize_depth", 
"decomposition_type" => "exact",
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
