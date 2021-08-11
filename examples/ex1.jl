import QuantumCircuitOpt as QCO
using JuMP
using CPLEX

# using Cbc

include("solver.jl")
include("2qubit_gates.jl")
include("3qubit_gates.jl")

#-------------------------------#
#      User-defined inputs      #
#-------------------------------#
function input_circuit()
         # [(depth, gate)]
    return [(1, "CNot_21"), 
            (2, "S_1"), 
            (3, "H_2"), 
            (4, "S_2")
            ]
end

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
        return Array{Complex{Float64},2}([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1]) 
    end
end

params = Dict{String, Any}(
"num_qubits" => 3,
"depth" => 7,

# "elementary_gates" => ["T_1", "T_2", "T_3", "H_3", "CNot_12", "CNot_13", "CNot_23", "Tdagger_2", "Tdagger_3", "Identity"],
"elementary_gates" => ["CV_12", "CV_23", "CV_13", "CVdagger_12", "CVdagger_23", "CVdagger_13", "CNot_12", "CNot_32", "CNot_23", "CNot_13", "Identity"],

"target_gate" => QCO.CSwapGate(),

# "input_circuit" => toffoli_circuit(),

# "RX_discretization" => [0, π/4],
# "RY_discretization" => [π/4],
# "RZ_discretization" => [2*π],

# "U_θ_discretization" => [-π/4, 0, π/4],
# "U_ϕ_discretization" => [-π/2, 0, π/2],
# "U_λ_discretization" => [-π/2, 0, π/4],

"objective" => "minimize_depth", 
"optimizer" => "cplex",
# "optimizer_presolve" => false                           
)

#------------------------------#
#      Optimization model      #
#------------------------------#
qcm_optimizer = get_solver(params)

# data = QCO.get_data(params)
result_qc = QCO.run_QCModel(params, 
                            qcm_optimizer, 
                            model_type = "compact_formulation")
