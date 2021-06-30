import QuantumCircuitOpt as QCO
using JuMP
using CPLEX
# using Cbc

include("solver.jl")

#-------------------------------#
#      User-defined inputs      #
#-------------------------------#
function target_gate(gate::Int)
    if gate == 1
        return QCO.kron_single_gate(2, QCO.U3Gate(0,0,π/4), "q1")
    elseif gate == 2
        return QCO.kron_single_gate(2, QCO.RZGate(π/4), "q1")
    elseif gate == 3
        return QCO.C2SXGate()
    elseif gate == 4
        return QCO.kron_single_gate(3, QCO.RXGate(π/4), "q3")
    end
end

params = Dict{String, Any}(
"num_qubits" => 3,
"depth" => 4,

# "elementary_gates" => ["RX", "RY", "RZ", "Identity"],
# "elementary_gates" => ["U3", "Identity", "cnot_12"],
"elementary_gates" => ["U3", "Identity"],
# "elementary_gates" => ["H1", "H2", "T1", "T2", "Tdagger1", "Tdagger2", "cnot_12", "Identity"],  

"target_gate" => target_gate(4),

"RX_discretization" => [π/4],
"RY_discretization" => [-π/4, π/4, π/2, -π/2, -π],
"RZ_discretization" => [-π/2, π/2, π/4, -π/2, -π],

"U_θ_discretization" => [0, π/4],
"U_ϕ_discretization" => [0, -π/2],
"U_λ_discretization" => [0, π/2],

"objective" => "minimize_depth", 
"decomposition_type" => "exact",
"optimizer" => "cplex"
                            
)

#------------------------------#
#      Optimization model      #
#------------------------------#
qcm_optimizer = get_solver(params)
result_qc = QCO.run_QCModel(params, qcm_optimizer, model_type = "compact_formulation")
