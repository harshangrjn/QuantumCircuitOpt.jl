function RX_on_5qubits()

    println(">>>>> RX Gate on fourth qubit using U3Gate <<<<<")
 
    return Dict{String, Any}(
    
    "num_qubits" => 5, 
    "maximum_depth" => 3,
    "elementary_gates" => ["H_1xCNot_2_3xI_4xI_5", "RX_2", "CU3_3_4", "Identity"], 
    "target_gate" => QCOpt.kron_two_qubit_gate(5, QCOpt.CRXGate(π/4), "q3", "q4"),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal",

    "CU3_θ_discretization" => [0, π/4],
    "CU3_ϕ_discretization" => [0, -π/2],
    "CU3_λ_discretization" => [0, π/2],
    "RX_discretization" => [π/2],
    )
end

# function test_multi_controlled_gate()

#     println(">>>>> Testing Multi controlled gate <<<<<")

#     num_qubits = 5

#     function target_gate()
#         CV_1_4 = QCOpt.get_unitary("CV_1_4", num_qubits);
#         CVdagger_1_4 = QCOpt.get_unitary("CVdagger_1_4", num_qubits);
#         CV_4_5 = QCOpt.get_unitary("CV_4_5", num_qubits);
#         CVdagger_4_5 = QCOpt.get_unitary("CVdagger_4_5", num_qubits);
#         CV_2_4 = QCOpt.get_unitary("CV_2_4", num_qubits);
#         CVdagger_2_4 = QCOpt.get_unitary("CVdagger_2_4", num_qubits);
#         CX_1_2 = QCOpt.get_unitary("CX_1_2", num_qubits);
#         CX_3_4 = QCOpt.get_unitary("CX_3_4", num_qubits);

#         return CV_4_5 * CX_3_4 * CVdagger_4_5 * CV_2_4 * CX_1_2 * CVdagger_2_4 * 
#                CV_1_4 * CV_4_5 * CX_3_4 * CVdagger_4_5 * CVdagger_1_4 * 
#                CV_2_4 * CX_1_2 * CVdagger_2_4
#     end
 
#     return Dict{String, Any}(
    
#     "num_qubits" => num_qubits, 
#     "maximum_depth" => 4,
#     "elementary_gates" => ["MC_gates", "Identity"], 
#     "MC_gates" => Dict{String, Any}("MC_gate_1" => Dict{String, Any}("control" => [1,2], "target" => 4),
#                                     "MC_gate_2" => Dict{String, Any}("control" => [3,4], "target" => 5)
#                                    ),
#     "target_gate" => target_gate(),
#     "objective" => "minimize_depth", 
#     "decomposition_type" => "exact_optimal",
#     )

# end