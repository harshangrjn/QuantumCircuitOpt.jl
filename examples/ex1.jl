import QuantumCircuitOpt as QCO
using JuMP
using CPLEX
using Cbc

include("solver.jl")

#-------------------------------#
#      User-defined inputs      #
#-------------------------------#
function target_gate(gate::Int)
    if gate == 1
        return QCO.kron_single_gate(2, QCO.U3Gate(0,0,π/4), "q1")
    elseif gate == 2
        return QCO.kron_single_gate(2, QCO.RXGate(π/4), "q1") * QCO.kron_single_gate(2, QCO.RYGate(π/4), "q2") * QCO.kron_single_gate(2, QCO.RZGate(π/4), "q1")
    elseif gate == 3
        return QCO.SwapGate()
    elseif gate == 4
        return QCO.kron_single_gate(3, QCO.RXGate(π/4), "q3")
    end
end

params = Dict{String, Any}(
"num_qubits" => 3,
"depth" => 2,

# "elementary_gates" => ["RX", "RY", "RZ", "Identity"],
"elementary_gates" => ["U3", "Identity"],
# "elementary_gates" => ["cnot_12", "cnot_21", "Identity"],
# "elementary_gates" => ["H1", "H2", "T1", "T2", "Tdagger1", "Tdagger2", "cnot_12", "Identity"],  

"target_gate" => target_gate(4),

"RX_discretization" => [0, π/4],
"RY_discretization" => [π/4],
"RZ_discretization" => [π/2, π/4],

"U_θ_discretization" => [0, π/4],
"U_ϕ_discretization" => [0, -π/2],
"U_λ_discretization" => [0, π/2],    

"objective" => "minimize_cnot", 
"decomposition_type" => "exact",
"optimizer" => "cbc"
                            
)

#------------------------------#
#      Optimization model      #
#------------------------------#
qcm_optimizer = get_solver(params)
result_qc = QCO.run_QCModel(params, qcm_optimizer, model_type = "balas_formulation", eliminate_identical_gates = true)
