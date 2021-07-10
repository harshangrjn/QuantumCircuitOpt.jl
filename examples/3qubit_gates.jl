function test_RX_on_q3()

    println(">>>>> RX Gate on third qubit using U3Gate <<<<<")
 
    params = Dict{String, Any}(
    
    "num_qubits" => 3, 
    "depth" => 3,    

    "elementary_gates" => ["U3", "Identity"], 
    "target_gate" => QCO.kron_single_gate(3, QCO.RXGate(π/4), "q3"),
       
    "U_θ_discretization" => [0, π/4],
    "U_ϕ_discretization" => [0, -π/2],
    "U_λ_discretization" => [0, π/2],    
 
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cplex",
    
    )

    return params
    
end

function test_toffoli()

    println(">>>>> Toffoli gate <<<<<")
 
    params = Dict{String, Any}(
    
    "num_qubits" => 3, 
    "depth" => 15,    

    "elementary_gates" => ["T1", "T2", "T3", "H3", "cnot_12", "cnot_13", "cnot_23", "Tdagger2", "Tdagger3", "Identity"], 
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
    return [(1, "H3"), 
        (2, "cnot_23"), 
        (3, "Tdagger3"), 
        (4, "cnot_13"), 
        (5, "T3"), 
        (6, "cnot_23"), 
        (7, "Tdagger3"), 
        (8, "cnot_13"), 
        (9, "T2"), 
        (10, "T3"), 
        (11, "cnot_12"), 
        (12, "H3"), 
        (13, "T1"), 
        (14, "Tdagger2"),
        (15, "cnot_12")
        ] 
end