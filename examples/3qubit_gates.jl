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
    
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cplex",
    
    )

    return params
    
end