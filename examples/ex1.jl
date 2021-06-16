using QuantumCircuitOpt
using JuMP
using CPLEX
using LinearAlgebra

const QCO = QuantumCircuitOpt

include("solver.jl")

#-------------------------------#
#      User-defined inputs      #
#-------------------------------#
params = Dict{String, Any}(
"num_qubits" => 2, 
"depth" => 5,  

# "elementary_gates" => ["R_x", "Identity"], 
# "target_gate" => "test_R_x_1",

"elementary_gates" => ["U3", "Identity"],
"target_gate" => "test_U3_1",

"R_x_discretization" => [π/4], 
"R_y_discretization" => [-π/4, π/4, π/2, -π/2, -π], 
"R_z_discretization" => [-π/2, π/2, π/4, -π/4, -π], 
"U_θ_discretization" => [0, π/2],
"U_ϕ_discretization" => [0, π/2],
"U_λ_discretization" => [0, π/4],

"objective" => "minimize_depth", 
"decomposition_type" => "exact",
"optimizer" => "cplex"
                            
)

#------------------------------#
#      Optimization model      #
#------------------------------#
qcm_optimizer = get_solver(params)

result_qc = QCO.run_QCModel(params, qcm_optimizer, model_type = "balas_formulation")