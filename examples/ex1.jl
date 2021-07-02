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
        return QCO.kron_single_gate(2, QCO.RXGate(π/4), "q1") * QCO.kron_single_gate(2, QCO.RYGate(π/4), "q2") * QCO.kron_single_gate(2, QCO.RZGate(π/4), "q1")
    elseif gate == 3
        return QCO.iSwapGate()
    elseif gate == 4
        return QCO.kron_single_gate(3, QCO.RXGate(π/4), "q3")
    elseif gate == 5
        return Array{Complex{Float64},2}([-1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1]) 
    end
end

params = Dict{String, Any}(
"num_qubits" => 2,
"depth" => 3,

"elementary_gates" => ["RZ", "cnot_12", "Identity"],
# "elementary_gates" => ["T1", "T2", "H1", "H2", "S1", "S2", "cnot_12", "cnot_21", "Identity"],
# "elementary_gates" => ["H1", "H2", "T1", "T2", "Tdagger1", "Tdagger2", "cnot_12", "Identity"],  

"target_gate" => target_gate(5),

"RX_discretization" => [0, π/4],
"RY_discretization" => [π/4],
"RZ_discretization" => [2*π],

"U_θ_discretization" => [0, π/2, π],
"U_ϕ_discretization" => [0, π/4, -π],
"U_λ_discretization" => [0, π/2, π],    

"objective" => "minimize_depth", 
"decomposition_type" => "exact",
"optimizer" => "cplex"
                           
)

#------------------------------#
#      Optimization model      #
#------------------------------#
qcm_optimizer = get_solver(params)
result_qc = QCO.run_QCModel(params, 
                            qcm_optimizer, 
                            model_type = "compact_formulation", 
                            eliminate_identical_gates = true)