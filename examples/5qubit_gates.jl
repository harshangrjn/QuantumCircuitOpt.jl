function decompose_RX_on_5qubits()

    println(">>>>> RX Gate on fourth qubit using U3Gate <<<<<")
 
    params = Dict{String, Any}(
    
    "num_qubits" => 5, 
    "maximum_depth" => 3,    

    "elementary_gates" => ["H_1xCNot_2_3xI_4xI_5", "RX_2", "CU3_3_4", "Identity"], 
    "target_gate" => QCO.kron_two_qubit_gate(5, QCO.CRXGate(π/4), "q3", "q4"),
       
    "CU3_θ_discretization" => [0, π/4],
    "CU3_ϕ_discretization" => [0, -π/2],
    "CU3_λ_discretization" => [0, π/2],    

    "RX_discretization" => [π/2],
 
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "gurobi"
    )

    return params
    
end