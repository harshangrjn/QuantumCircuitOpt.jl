# Testing file for creating useful circuit decompositions

import QuantumCircuitOpt as QCOpt
using JuMP
using CPLEX
using Gurobi

# using Cbc
include("optimizers.jl")
include("2qubit_gates.jl")
include("3qubit_gates.jl")

#-------------------------------#
#      User-defined inputs      #
#-------------------------------#
function input_circuit()
         # [(depth, gate)]
    return [(1, "CNot_2_1"), 
            (2, "S_1"), 
            (3, "H_2"), 
            (4, "S_2")
            ]
end

function target_gate(gate::Int)
    if gate == 1
        return QCOpt.kron_single_qubit_gate(2, QCOpt.U3Gate(0,0,π/4), "q1")
    elseif gate == 2
        return QCOpt.kron_single_qubit_gate(2, QCOpt.RXGate(π/4), "q1") * QCOpt.kron_single_qubit_gate(2, QCOpt.RYGate(π/4), "q2") * QCOpt.kron_single_qubit_gate(2, QCOpt.RZGate(π/4), "q1")
    elseif gate == 3
        return QCOpt.iSwapGate()
    elseif gate == 4
        return QCOpt.kron_single_qubit_gate(3, QCOpt.RXGate(π/4), "q3")
    elseif gate == 5
        return Array{Complex{Float64},2}([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1])
    end
end

params = Dict{String, Any}(
"num_qubits" => 3,
"maximum_depth" => 7,
"elementary_gates" => ["RY_3", "CNot_1_2", "CNot_2_3", "CNot_1_3", "Identity"],
# "elementary_gates" => ["CV_1_2", "CV_2_3", "CV_1_3", "CVdagger_1_2", "CVdagger_2_3", "CVdagger_1_3", "CNot_1_2", "CNot_3_2", "CNot_2_3", "CNot_1_3", "Identity"],
"target_gate" => target_gate(5),
"objective" => "minimize_depth", 

"RY_discretization" => -π:π/4:π,

# "input_circuit" => toffoli_circuit(),
)

#------------------------------#
#      Optimization model      #
#------------------------------#
qcopt_optimizer = get_gurobi()

# data = QCOpt.get_data(params)
result_qc = QCOpt.run_QCModel(params, qcopt_optimizer)
