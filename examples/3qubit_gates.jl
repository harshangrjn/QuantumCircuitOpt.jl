function decompose_RX_on_q3()

    println(">>>>> RX Gate on third qubit using U3Gate <<<<<")
 
    params = Dict{String, Any}(
    
    "num_qubits" => 3, 
    "depth" => 3,    

    "elementary_gates" => ["U3_3", "Identity"], 
    "target_gate" => QCO.kron_single_qubit_gate(3, QCO.RXGate(π/4), "q3"),
       
    "U3_θ_discretization" => [0, π/4],
    "U3_ϕ_discretization" => [0, -π/2],
    "U3_λ_discretization" => [0, π/2],    
 
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cplex",
    "optimizer_presolve" => false, #turning this true will give infeasiblity in cplex - most probably a bug in cplex's presolve
    
    )

    return params
    
end

function decompose_toffoli()

    println(">>>>> Toffoli gate <<<<<")
 
    params = Dict{String, Any}(
    
    "num_qubits" => 3, 
    "depth" => 15,    

    "elementary_gates" => ["T_1", "T_2", "T_3", "H_3", "CNot_1_2", "CNot_1_3", "CNot_2_3", "Tdagger_1", "Tdagger_2", "Tdagger_3", "Identity"], 
    "target_gate" => QCO.ToffoliGate(),
    "input_circuit" => toffoli_circuit(),
    
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    "relax_integrality" => false,
    
    "optimizer" => "cplex",
    
    )

    return params
    
end

function decompose_toffoli_using_kronecker()

    println(">>>>> Toffoli gate using Kronecker <<<<<")
 
    params = Dict{String, Any}(
    
    "num_qubits" => 3, 
    "depth" => 12,    
  
    "elementary_gates" => ["T_3", "H_3", "CNot_1_2", "CNot_1_3", "CNot_2_3", "Tdagger_3", "I_1xT_2xT_3", "CNot_1_2xH_3", "T_1xTdagger_2xI_3", "Identity"], 
    "target_gate" => QCO.ToffoliGate(),
    
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    "relax_integrality" => false,
    
    "optimizer" => "cplex",
    
    )

    return params
    
end

function decompose_toffoli_using_Rotations()

    println(">>>>> Toffoli gate using rotations <<<<<")
 
    params = Dict{String, Any}(
    
    "num_qubits" => 3,
    "depth" => 10,

    "elementary_gates" => ["RZ_3", "CNot_1_3", "CNot_2_3", "Identity"],
    "RZ_discretization" => [-π/2, π/2, π/4],

    "target_gate" => QCO.ToffoliGate(),
    
    "objective" => "minimize_depth",
    "decomposition_type" => "exact",
    "relax_integrality" => false,
    
    "optimizer" => "cplex",
    )

    return params
    
end

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

function decompose_CNot_13()

    params = Dict{String, Any}(
    "num_qubits" => 3,
    "depth" => 8,

    # "elementary_gates" => ["CNot_1_2", "CNot_2_3", "Identity"], 
    "elementary_gates" => ["H_1", "H_3", "H_2", "CNot_2_1", "CNot_3_2", "Identity"],
    "target_gate" => QCO.get_full_sized_gate("CNot_1_3", 3),

    "objective" => "minimize_depth", 
    "optimizer" => "cplex"   
    )

    return params
end

function decompose_FredkinGate()

    params = Dict{String, Any}(
    "num_qubits" => 3,
    "depth" => 7,

    # Reference: https://doi.org/10.1103/PhysRevA.53.2855
    "elementary_gates" => ["CV_1_2", "CV_2_3", "CV_1_3", "CVdagger_1_2", "CVdagger_2_3", "CVdagger_1_3", "CNot_1_2", "CNot_3_2", "CNot_2_3", "CNot_1_3", "Identity"],
    "target_gate" => QCO.CSwapGate(), #also Fredkin

    "objective" => "minimize_depth", 
    "optimizer" => "cplex"
    )
    
    return params
end

function decompose_toffoli_left()

    # Reference: https://arxiv.org/pdf/0803.2316.pdf 

    function target_gate()
        H_3 = QCO.get_full_sized_gate("H_3", 3)
        Tdagger_3 = QCO.get_full_sized_gate("Tdagger_3", 3)
        T_3 = QCO.get_full_sized_gate("T_3", 3)
        cnot_23 = QCO.get_full_sized_gate("CNot_2_3", 3)
        cnot_13 = QCO.get_full_sized_gate("CNot_1_3", 3)

        return H_3 * cnot_23 * Tdagger_3 * cnot_13 * T_3 * cnot_23 * Tdagger_3
    end

    params = Dict{String, Any}(
    "num_qubits" => 3,
    "depth" => 7,

    "elementary_gates" => ["H_3", "T_1", "T_2", "T_3", "Tdagger_1", "Tdagger_2", "Tdagger_3", "CNot_1_2", "CNot_2_3", "CNot_1_3", "Identity"],
    "target_gate" => target_gate(), 

    "objective" => "minimize_depth", 
    "optimizer" => "cplex"
    )
    
    return params
end

function decompose_toffoli_right()

    # Reference: https://arxiv.org/pdf/0803.2316.pdf 
    
    function target_gate()
        H_3 = QCO.get_full_sized_gate("H_3", 3)
        Tdagger_2 = QCO.get_full_sized_gate("Tdagger_2", 3)
        T_1 = QCO.get_full_sized_gate("T_1", 3)
        T_2 = QCO.get_full_sized_gate("T_2", 3)
        T_3 = QCO.get_full_sized_gate("T_3", 3)
        cnot_12 = QCO.get_full_sized_gate("CNot_1_2", 3)
        cnot_13 = QCO.get_full_sized_gate("CNot_1_3", 3)

        return cnot_13 * T_2 * T_3 * cnot_12 * H_3 * T_1 * Tdagger_2 * cnot_12
    end

    params = Dict{String, Any}(
    "num_qubits" => 3,
    "depth" => 8,

    "elementary_gates" => ["H_3", "T_1", "T_2", "T_3", "Tdagger_1", "Tdagger_2", "Tdagger_3", "CNot_1_2", "CNot_2_3", "CNot_1_3", "Identity"],
    "target_gate" => target_gate(), 

    "objective" => "minimize_depth", 
    "optimizer" => "cplex"
    )
    
    return params
end