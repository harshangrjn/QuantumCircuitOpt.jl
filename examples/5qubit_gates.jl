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

function decompose_exact_random_target1()

    num_qubits = 5
    elementary_gates = ["H_1", "H_2", "Y_3", "T_4", "T_5", "CH_1_2", "CNot_2_3", "CNot_3_4", "CH_4_5", "CZ_1_5", "Identity"]
    depth = 5

    function target_gate()

        H_1 = QCOpt.get_full_sized_gate("H_1", num_qubits);
        H_2 = QCOpt.get_full_sized_gate("H_2", num_qubits);
        H_3 = QCOpt.get_full_sized_gate("H_3", num_qubits);
        Y_3 = QCOpt.get_full_sized_gate("Y_3", num_qubits);
        T_1 = QCOpt.get_full_sized_gate("T_1", num_qubits);
        T_3 = QCOpt.get_full_sized_gate("T_3", num_qubits);
        T_4 = QCOpt.get_full_sized_gate("T_4", num_qubits);
        T_5 = QCOpt.get_full_sized_gate("T_5", num_qubits);
    
        CH_1_2 = QCOpt.get_full_sized_gate("CH_1_2", num_qubits);
        CH_4_5 = QCOpt.get_full_sized_gate("CH_4_5", num_qubits);
        CNot_2_3 = QCOpt.get_full_sized_gate("CNot_2_3", num_qubits);
        CNot_3_4 = QCOpt.get_full_sized_gate("CNot_3_4", num_qubits);
        CNot_1_5 = QCOpt.get_full_sized_gate("CNot_1_5", num_qubits);
        CZ_1_5 = QCOpt.get_full_sized_gate("CZ_1_5", num_qubits);
        CVdagger_1_2 = QCOpt.get_full_sized_gate("CVdagger_1_2", num_qubits);
        CVdagger_4_2 = QCOpt.get_full_sized_gate("CVdagger_4_2", num_qubits);

        return Y_3 * CNot_2_3 * Y_3 * CH_4_5 * H_1

    end

    return Dict{String, Any}(
    
    "num_qubits" => num_qubits, 
    "maximum_depth" => (depth + 1),
    "elementary_gates" => elementary_gates, 
    "target_gate" => target_gate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal"    
    )

end