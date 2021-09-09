function decompose_RX_on_q3()

    println(">>>>> RX Gate on third qubit using U3Gate <<<<<")
 
    params = Dict{String, Any}(
    
    "num_qubits" => 3, 
    "depth" => 3,    

    "elementary_gates" => ["U3_3", "Identity"], 
    "target_gate" => QCO.kron_single_qubit_gate(3, QCO.RXGate(π/4), "q3"),
       
    "U_θ_discretization" => [0, π/4],
    "U_ϕ_discretization" => [0, -π/2],
    "U_λ_discretization" => [0, π/2],    
 
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

    "elementary_gates" => ["T_1", "T_2", "T_3", "H_3", "CNot_12", "CNot_13", "CNot_23", "Tdagger_2", "Tdagger_3", "Identity"], 
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
  
    "elementary_gates" => ["T_3", "H_3", "CNot_12", "CNot_13", "CNot_23", "Tdagger_3", "I_1xT_2xT_3", "CNot_12xH_3", "T_1xTdagger_2xI_3", "Identity"], 
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

    "elementary_gates" => ["RZ_3", "CNot_13", "CNot_23", "Identity"],
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
            (4, "CNot_23"),           
            (5, "Tdagger_3"),        
            (6, "CNot_13"),          
            (7, "T_3"),               
            (8, "CNot_23"),                 
            (9, "Tdagger_3"),          
            (10, "CNot_13"),           
            (11, "T_3"),             
            (12, "H_3"),             
            (13, "CNot_12"),         
            (14, "Tdagger_2"),       
            (15, "CNot_12")          
            ] 
end

function decompose_CNot_13()

    params = Dict{String, Any}(
    "num_qubits" => 3,
    "depth" => 8,

    # "elementary_gates" => ["CNot_12", "CNot_23", "Identity"],
    "elementary_gates" => ["H_1", "H_3", "H_2", "CNot_21", "CNot_32", "Identity"],
    "target_gate" => QCO.get_full_sized_gate("CNot_13", 3),

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
    "elementary_gates" => ["CV_12", "CV_23", "CV_13", "CVdagger_12", "CVdagger_23", "CVdagger_13", "CNot_12", "CNot_32", "CNot_23", "CNot_13", "Identity"],
    "target_gate" => QCO.CSwapGate(), #also Fredkin

    "objective" => "minimize_depth", 
    "optimizer" => "cplex"
    )
    
    return params
end