function CNot_41()
    # Reference: https://doi.org/10.1109/DSD.2018.00005

    return Dict{String, Any}(
        
    "num_qubits" => 4,
    "maximum_depth" => 10,
    "elementary_gates" => ["H_1", "H_2", "H_3", "CNot_1_3", "CNot_4_3", "Identity"],
    "target_gate" => QCOpt.unitary("CNot_4_1", 4),
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal"
    )
end

function quantum_fulladder()
    #Reference-1: https://doi.org/10.1109/DATE.2005.249
    #Reference-2: https://doi.org/10.1109/TCAD.2005.858352 
    
    num_qubits = 4

    function target_gate_1()

        CV_1_2 = QCOpt.unitary("CV_1_2", num_qubits);
        CV_4_2 = QCOpt.unitary("CV_4_2", num_qubits);
        CVdagger_1_2 = QCOpt.unitary("CVdagger_1_2", num_qubits);
        CVdagger_3_2 = QCOpt.unitary("CVdagger_3_2", num_qubits);
        CNot_3_1 = QCOpt.unitary("CNot_3_1", num_qubits);
        CNot_4_3 = QCOpt.unitary("CNot_4_3", num_qubits);
        CNot_2_4 = QCOpt.unitary("CNot_2_4", num_qubits);

        return CNot_2_4 * CV_4_2 * CVdagger_1_2 * CVdagger_3_2 * CNot_4_3 * CNot_3_1 * CV_1_2
    end

    return Dict{String, Any}(

    "num_qubits" => num_qubits,
    "maximum_depth" => 7,
    "elementary_gates" => ["CV_1_2", "CV_4_2", "CV_3_2", "CVdagger_1_2", "CVdagger_4_2", "CVdagger_3_2", "CNot_3_1", "CNot_4_3", "CNot_2_4", "CNot_4_1", "Identity"],
    "target_gate" => target_gate_1(),
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal"
    )
end

function double_toffoli()

    # Reference: https://doi.org/10.1109/TCAD.2005.858352 
    
    num_qubits = 4

    function target_gate()
        CV_3_4 = QCOpt.unitary("CV_3_4", num_qubits);
        CV_2_4 = QCOpt.unitary("CV_2_4", num_qubits);
        CVdagger_2_4 = QCOpt.unitary("CVdagger_2_4", num_qubits);        
        CNot_1_3 = QCOpt.unitary("CNot_1_3", num_qubits);
        CNot_3_2 = QCOpt.unitary("CNot_3_2", num_qubits);

        return CV_2_4 * CNot_1_3 * CNot_3_2 * CV_3_4 * CVdagger_2_4 * CNot_3_2 * CNot_1_3
    end

    return Dict{String, Any}(

    "num_qubits" => 4,
    "maximum_depth" => 7,
    "elementary_gates" => ["CV_1_2", "CV_2_4", "CV_3_4", "CVdagger_1_2", "CVdagger_2_4", "CVdagger_3_4", "CNot_1_3", "CNot_3_2", "CNot_2_3", "Identity"],
    "target_gate" => target_gate(),
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal"
    )
end

function double_peres()

    # Reference: https://doi.org/10.1109/TCAD.2005.858352 
    
    num_qubits = 4

    function target_gate()
        CV_1_4 = QCOpt.unitary("CV_1_4", num_qubits);
        CV_2_4 = QCOpt.unitary("CV_2_4", num_qubits);
        CV_3_4 = QCOpt.unitary("CV_3_4", num_qubits);
        CVdagger_3_4 = QCOpt.unitary("CVdagger_3_4", num_qubits);        
        CNot_1_2 = QCOpt.unitary("CNot_1_2", num_qubits);
        CNot_2_3 = QCOpt.unitary("CNot_2_3", num_qubits);

        return CV_2_4 * CNot_1_2 * CV_3_4 * CV_1_4 * CNot_2_3 * CVdagger_3_4
    end

    return Dict{String, Any}(

    "num_qubits" => 4,
    "maximum_depth" => 7,
    "elementary_gates" => ["CV_1_4", "CV_2_4", "CV_3_4", "CVdagger_1_4", "CVdagger_2_4", "CVdagger_3_4", "CNot_1_2", "CNot_2_3", "Identity"],
    "target_gate" => target_gate(),
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal"
    )
end

function qubit_routing_circuit()

    # Reference: https://doi.org/10.1007/s10957-023-02229-w
    
    num_qubits = 4

    function target_gate()
        # Gates which do not satisfy qubit connectivity
        H_1 = QCOpt.unitary("H_1", num_qubits)
        S_1 = QCOpt.unitary("S_1", num_qubits)
        X_4 = QCOpt.unitary("X_4", num_qubits)
        CX_2_4 = QCOpt.unitary("CX_2_4", num_qubits);
        CY_2_3 = QCOpt.unitary("CY_2_3", num_qubits);
        CZ_1_2 = QCOpt.unitary("CZ_1_2", num_qubits);
        
        return H_1 * CX_2_4 * CZ_1_2 * X_4 * S_1 * CY_2_3
    end

    return Dict{String, Any}(

    "num_qubits" => 4,
    "maximum_depth" => 8,
    # gates which satisfy qubit connectivity :  |1⟩ -> |2⟩ -> |3⟩ -> |4⟩
    "elementary_gates" => ["H_1", "S_1", "X_3", "X_4", "CY_2_3", "CX_2_3", "CX_3_2", "CX_3_4", "CX_4_3", "CZ_1_2", "Identity"],
    "target_gate" => target_gate(),
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal"
    )
end