using QuantumCircuitOptimization
using JuMP
using CPLEX
using Cbc

include("solver.jl")

# All the user 
params = Dict{String, Any}(
"n" => 2, # Number of qubits
"D" => 3, # Minimum depth of the decomposition (>= 2)
"instance" => "ibm", 

"optimizer" => "cplex",
"presolve" => true,
"optimizer_log" => true,                           
"lp_relax" => false,
                            
# Valid inequalities
"valid_1" => false, #commutative matrices
)

qcm_optimizer = get_solver(params)
model_qc = build_QCModel(data,params)
results = run_QCModel(model_qc, optimizer = cplex)


