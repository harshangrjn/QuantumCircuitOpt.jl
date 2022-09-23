function decompose_RX_on_5qubits()

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

function verification_of_cat_state()

    println(">>>>> verification_of_cat_state <<<<<")
 
    num_qubits = 7

    function target_gate()

        CNot_4_5 = QCOpt.get_full_sized_gate("CNot_4_5", num_qubits);
        CNot_3_5 = QCOpt.get_full_sized_gate("CNot_3_5", num_qubits);
        CNot_3_6 = QCOpt.get_full_sized_gate("CNot_3_6", num_qubits);
        CNot_2_6 = QCOpt.get_full_sized_gate("CNot_2_6", num_qubits);
        CNot_2_7 = QCOpt.get_full_sized_gate("CNot_2_7", num_qubits);
        CNot_1_7 = QCOpt.get_full_sized_gate("CNot_1_7", num_qubits);
        

        return CNot_4_5 * CNot_3_5 * CNot_3_6 * CNot_2_6 * CNot_2_7 * CNot_1_7
    end
    return Dict{String, Any}(
    
    "num_qubits" => 7, 
    "maximum_depth" => 6,
    "elementary_gates" => ["CNot_1_5", "CNot_2_5", "CNot_2_6", "CNot_3_6", "CNot_3_7", "CNot_4_7", "CNot_5_6", "CNot_5_7", "CNot_6_7" ], 
    "target_gate" => target_gate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal",

    )
end