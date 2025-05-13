function RX_on_q3()

    println(">>>>> RX Gate on third qubit using U3Gate <<<<<")
 
    return Dict{String, Any}(
    "num_qubits" => 3, 
    "maximum_depth" => 3,
    "elementary_gates" => ["U3_3", "Identity"], 
    "target_gate" => QCOpt.kron_single_qubit_gate(3, QCOpt.RXGate(π/4), "q3"),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal",
       
    "U3_θ_discretization" => [0, π/4],
    "U3_ϕ_discretization" => [0, -π/2],
    "U3_λ_discretization" => [0, π/2],
    )
end

function toffoli()

    function toffoli_circuit()
        # [(depth, gate)]
        return [(1, "T_1"),              
                (2, "T_2"),                   
                (3, "H_3"),              
                (4, "CNot_2_3"),           
                (5, "Tdagger_3"),        
                (6, "CNot_1_3"),          
                (7, "T_3"),               
                (8, "CNot_2_3"),                 
                (9, "Tdagger_3"),          
                (10, "CNot_1_3"),           
                (11, "T_3"),             
                (12, "H_3"),             
                (13, "CNot_1_2"),         
                (14, "Tdagger_2"),       
                (15, "CNot_1_2")          
                ] 
    end

    println(">>>>> Toffoli gate <<<<<")
 
    return Dict{String, Any}(
    "num_qubits" => 3, 
    "maximum_depth" => 15,
    "elementary_gates" => ["T_1", "T_2", "T_3", "H_3", "CNot_1_2", "CNot_1_3", "CNot_2_3", "Tdagger_1", "Tdagger_2", "Tdagger_3", "Identity"], 
    "target_gate" => QCOpt.ToffoliGate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal",

    "input_circuit" => toffoli_circuit(),
    "set_cnot_lower_bound" => 6,
    "set_cnot_upper_bound" => 6,
    )
end

function toffoli_using_kronecker()

    println(">>>>> Toffoli gate using Kronecker <<<<<")
 
    return Dict{String, Any}(
    "num_qubits" => 3,
    "maximum_depth" => 12,
    "elementary_gates" => ["T_3", "H_3", "CNot_1_2", "CNot_1_3", "CNot_2_3", "Tdagger_3", "I_1xT_2xT_3", "CNot_1_2xH_3", "T_1xTdagger_2xI_3", "Identity"], 
    "target_gate" => QCOpt.ToffoliGate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal",

    "set_cnot_lower_bound" => 6,
    )
end

function toffoli_using_2qubit_gates()

    # Reference: https://doi.org/10.1109/TCAD.2005.858352
    println(">>>>> Toffoli gate with controlled gates <<<<<")

    return Dict{String, Any}(
    "num_qubits" => 3,
    "maximum_depth" => 5,
    "elementary_gates" => ["CV_1_3", "CV_2_3", "CV_1_2", "CVdagger_1_3", "CVdagger_2_3", "CVdagger_1_2", "CNot_2_1", "CNot_1_2", "Identity"],
    "target_gate" => QCOpt.ToffoliGate(),
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal",
    )
end

function CNot_1_3()

    return Dict{String, Any}(
    "num_qubits" => 3,
    "maximum_depth" => 8,
    # "elementary_gates" => ["CNot_1_2", "CNot_2_3", "Identity"], 
    "elementary_gates" => ["H_1", "H_3", "H_2", "CNot_2_1", "CNot_3_2", "Identity"],
    "target_gate" => QCOpt.unitary("CNot_1_3", 3),
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal",
    )
end

function fredkin()

    return Dict{String, Any}(
    "num_qubits" => 3,
    "maximum_depth" => 7,
    # Reference: https://doi.org/10.1103/PhysRevA.53.2855
    "elementary_gates" => ["CV_1_2", "CV_2_3", "CV_1_3", "CVdagger_1_2", "CVdagger_2_3", "CVdagger_1_3", "CNot_1_2", "CNot_3_2", "CNot_2_3", "CNot_1_3", "Identity"],
    "target_gate" => QCOpt.CSwapGate(), #also Fredkin
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal"
    )
end

function fredkin_using_SSwap()

    return Dict{String, Any}(
    "num_qubits" => 3,
    "maximum_depth" => 8,
    "elementary_gates" => ["SSwap_2_3", "SX_1", "SX_2", "SX_3", "H_1", "H_2", "H_3", "CNot_1_2", "CNot_3_2", "CNot_2_3", "CNot_1_3", "Identity"],
    "target_gate" => QCOpt.CSwapGate(), #also Fredkin
    "decomposition_type" => "exact_optimal",
    "objective" => "minimize_depth"
    )
end

function toffoli_left()

    # Reference: https://arxiv.org/pdf/0803.2316.pdf 

    function target_gate()
        H_3 = QCOpt.unitary("H_3", 3)
        Tdagger_3 = QCOpt.unitary("Tdagger_3", 3)
        T_3 = QCOpt.unitary("T_3", 3)
        cnot_23 = QCOpt.unitary("CNot_2_3", 3)
        CNot_1_3 = QCOpt.unitary("CNot_1_3", 3)

        return H_3 * cnot_23 * Tdagger_3 * CNot_1_3 * T_3 * cnot_23 * Tdagger_3
    end

    return Dict{String, Any}(
    "num_qubits" => 3,
    "maximum_depth" => 7,
    "elementary_gates" => ["H_3", "T_1", "T_2", "T_3", "Tdagger_1", "Tdagger_2", "Tdagger_3", "CNot_1_2", "CNot_2_3", "CNot_1_3", "Identity"],
    "target_gate" => target_gate(),
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal"
    )
end

function toffoli_right()

    # Reference: https://arxiv.org/pdf/0803.2316.pdf 
    
    function target_gate()
        H_3 = QCOpt.unitary("H_3", 3)
        Tdagger_2 = QCOpt.unitary("Tdagger_2", 3)
        T_1 = QCOpt.unitary("T_1", 3)
        T_2 = QCOpt.unitary("T_2", 3)
        T_3 = QCOpt.unitary("T_3", 3)
        CNot_1_2 = QCOpt.unitary("CNot_1_2", 3)
        CNot_1_3 = QCOpt.unitary("CNot_1_3", 3)

        return CNot_1_3 * T_2 * T_3 * CNot_1_2 * H_3 * T_1 * Tdagger_2 * CNot_1_2
    end

    return Dict{String, Any}(
    "num_qubits" => 3,
    "maximum_depth" => 8,
    "elementary_gates" => ["H_3", "T_1", "T_2", "T_3", "Tdagger_1", "Tdagger_2", "Tdagger_3", "CNot_1_2", "CNot_2_3", "CNot_1_3", "Identity"],
    "target_gate" => target_gate(), 
    "objective" => "minimize_depth"
    )
end

function miller()
    # Reference: https://doi.org/10.1109/TCAD.2005.858352

    function target_gate()
        CV_1_3 = QCOpt.unitary("CV_1_3", 3)
        CV_2_3 = QCOpt.unitary("CV_2_3", 3)
        CVdagger_2_3 = QCOpt.unitary("CVdagger_2_3", 3)
        CNot_1_2 = QCOpt.unitary("CNot_1_2", 3)
        CNot_3_1 = QCOpt.unitary("CNot_3_1", 3)
        CNot_3_2 = QCOpt.unitary("CNot_3_2", 3)

        return CNot_3_1 * CNot_3_2 * CV_2_3 * CNot_1_2 * CVdagger_2_3 * CV_1_3 * CNot_3_1 * CNot_1_2
    end

    return Dict{String, Any}(
    "num_qubits" => 3,
    "maximum_depth" => 10,
    "elementary_gates" => ["CV_1_3", "CV_2_3", "CVdagger_2_3", "CNot_1_2", "CNot_3_1", "CNot_3_2", "Identity"],
    "target_gate" => target_gate(),
    "decomposition_type" => "exact_optimal",
    "objective" => "minimize_depth"
    )
end

function relative_toffoli()
    #Reference: https://arxiv.org/pdf/1508.03273.pdf

    return Dict{String, Any}(
        "num_qubits" => 3,
        "maximum_depth" => 9,
        "elementary_gates" => ["H_3", "T_3", "Tdagger_3", "CNot_1_2", "CNot_2_3", "CNot_1_3", "Identity"],
        "target_gate" => QCOpt.RCCXGate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "exact_optimal",

        "set_cnot_upper_bound" => 3
    )
end

function margolus()
    #Reference: https://arxiv.org/pdf/quant-ph/0312225.pdf

    return Dict{String, Any}(
        "num_qubits" => 3,
        "maximum_depth" => 7,
        "elementary_gates" => ["RY_3", "CNot_1_2", "CNot_2_3", "CNot_1_3", "Identity"],
        "target_gate" => QCOpt.MargolusGate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "exact_optimal",

        "RY_discretization" => -π:π/4:π,
        # "set_cnot_lower_bound" => 3,
        # "set_cnot_upper_bound" => 3
        )
end

function CiSwap()
    #Reference: https://doi.org/10.1103/PhysRevResearch.2.033097

    return Dict{String, Any}(
        "num_qubits" => 3,
        "maximum_depth" => 12,
        "elementary_gates" => ["H_3", "T_3", "Tdagger_3", "CNot_3_2", "CNot_2_3", "CNot_1_3", "Identity"],
        "target_gate" => QCOpt.CiSwapGate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "exact_optimal",
        )
end

function QFT3()

    # QFT3 circuit using Controlled-Phase gates
    return Dict{String, Any}(
        "num_qubits" => 3,
        "maximum_depth" => 7,   
        # "elementary_gates" => ["H_1", "H_2", "H_3", "CS_2_1", "CS_3_2", "CT_3_1", "Swap_1_3", "Identity"],
        "elementary_gates" => ["H_1", "H_2", "H_3", "CU3_2_1", "CU3_3_2", "CU3_3_1", "Swap_1_3", "Identity"],
        "CU3_θ_discretization" => [0],
        "CU3_ϕ_discretization" => [0],
        "CU3_λ_discretization" => 0:π/4:π,
        "target_gate" => QCOpt.QFT3Gate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "exact_optimal",
        )
end