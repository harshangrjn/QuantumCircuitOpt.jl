function CNot_41()
    # Reference: https://doi.org/10.1109/DSD.2018.00005

    return Dict{String, Any}(
        
    "num_qubits" => 4,
    "maximum_depth" => 10,
    "elementary_gates" => ["H_1", "H_2", "H_3", "CNot_1_3", "CNot_4_3", "Identity"],
    "target_gate" => QCOpt.get_full_sized_gate("CNot_4_1", 4),
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal"
    )
end

function quantum_fulladder()
    #Reference-1: https://doi.org/10.1109/DATE.2005.249
    #Reference-2: https://doi.org/10.1109/TCAD.2005.858352 
    
    num_qubits = 4

    function target_gate_1()

        CV_1_2 = QCOpt.get_full_sized_gate("CV_1_2", num_qubits);
        CV_4_2 = QCOpt.get_full_sized_gate("CV_4_2", num_qubits);
        CVdagger_1_2 = QCOpt.get_full_sized_gate("CVdagger_1_2", num_qubits);
        CVdagger_3_2 = QCOpt.get_full_sized_gate("CVdagger_3_2", num_qubits);
        CNot_3_1 = QCOpt.get_full_sized_gate("CNot_3_1", num_qubits);
        CNot_4_3 = QCOpt.get_full_sized_gate("CNot_4_3", num_qubits);
        CNot_2_4 = QCOpt.get_full_sized_gate("CNot_2_4", num_qubits);

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
        CV_3_4 = QCOpt.get_full_sized_gate("CV_3_4", num_qubits);
        CV_2_4 = QCOpt.get_full_sized_gate("CV_2_4", num_qubits);
        CVdagger_2_4 = QCOpt.get_full_sized_gate("CVdagger_2_4", num_qubits);        
        CNot_1_3 = QCOpt.get_full_sized_gate("CNot_1_3", num_qubits);
        CNot_3_2 = QCOpt.get_full_sized_gate("CNot_3_2", num_qubits);

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
        CV_1_4 = QCOpt.get_full_sized_gate("CV_1_4", num_qubits);
        CV_2_4 = QCOpt.get_full_sized_gate("CV_2_4", num_qubits);
        CV_3_4 = QCOpt.get_full_sized_gate("CV_3_4", num_qubits);
        CVdagger_3_4 = QCOpt.get_full_sized_gate("CVdagger_3_4", num_qubits);        
        CNot_1_2 = QCOpt.get_full_sized_gate("CNot_1_2", num_qubits);
        CNot_2_3 = QCOpt.get_full_sized_gate("CNot_2_3", num_qubits);

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

function decompose_0_reflection_on_4qubits()

    println(">>>>> Reflection of the |0><0| state <<<<<")

    num_qubits = 4

    function target_gate()
        return Array{Complex{Float64},2}(Diagonal([-1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1]))
    end

    return Dict{String, Any}(
        "num_qubits" => num_qubits,
        "maximum_depth" => 10,
        "elementary_gates" => ["H_1", "H_2", "H_3", "H_4", "Z_1", "Z_2", "Z_3", "Z_4",  "T_1", "T_2", "T_3", "T_4",  "Tdagger_1", "Tdagger_2", "Tdagger_3", "Tdagger_4",  "CZ_1_2", "CZ_1_3", "CZ_1_4", "CZ_2_3", "CZ_2_4", "CZ_3_4",  "Identity"],
        "target_gate" => target_gate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "exact_feasible")
end