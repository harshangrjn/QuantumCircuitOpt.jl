function decompose_RX_on_q3()

    println(">>>>> RX Gate on third qubit using U3Gate <<<<<")
 
    params = Dict{String, Any}(
    
    "num_qubits" => 3, 
    "depth" => 3,    

    "elementary_gates" => ["U3_3", "Identity"], 
    "target_gate" => QCO.kron_single_gate(3, QCO.RXGate(π/4), "q3"),
       
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
    "relax_integrality" => true,
    
    "optimizer" => "cplex",
    
    )

    return params
    
end

function toffoli_circuit()
    # [(depth, gate)]
    return [(1, "H_3"), 
            (2, "CNot_23"), 
            (3, "Tdagger_3"), 
            (4, "CNot_13"), 
            (5, "T_3"), 
            (6, "CNot_23"), 
            (7, "Tdagger_3"), 
            (8, "CNot_13"), 
            (9, "T_2"), 
            (10, "T_3"), 
            (11, "CNot_12"), 
            (12, "H_3"), 
            (13, "T_1"), 
            (14, "Tdagger_2"),
            (15, "CNot_12")
            ] 
end
