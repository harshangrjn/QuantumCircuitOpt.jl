function decompose_hadamard()

    println(">>>>> Hadamard Gate <<<<<")
 
    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 3,    

    "elementary_gates" => ["U3", "CNot_12", "Identity"], 
    "target_gate" => QCO.kron_single_gate(2, QCO.HGate(), "q1"),
       
    "U_θ_discretization" => [0, π/2],
    "U_ϕ_discretization" => [0, π/2],
    "U_λ_discretization" => [0, π],
 
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cplex",
    "MIP_feasiblity_emphasis" => true
    
    )

    return params
    
end

function decompose_controlled_Z()

    println(">>>>> Controlled-Z Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 4,    

    "elementary_gates" => ["U3", "CNot_12", "Identity"], 
    "target_gate" => QCO.CZGate(),
       
    "U_θ_discretization" => [-π/2, 0, π/2],
    "U_ϕ_discretization" => [-π/2, 0, π/2],
    "U_λ_discretization" => [-π/2, 0, π/2],

    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cplex",
    "optimizer_presolve" => false
                                
    )

    return params
    
end

function decompose_controlled_V()

    println(">>>>> Controlled-V Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 7,    

    "elementary_gates" => ["H_1", "H_2", "T_1", "T_2", "Tdagger_1", "CNot_12", "CNot_21"],
    "target_gate" => QCO.CVGate(),
    
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cplex"
                                
    )

    return params
    
end

function decompose_controlled_H()

    println(">>>>> Controlled-H Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 5,    

    "elementary_gates" => ["U3", "CNot_12", "Identity"], 
    "target_gate" => QCO.CHGate(),

    "U_θ_discretization" => [-π/4, 0, π/4],
    "U_ϕ_discretization" => [0],
    "U_λ_discretization" => [0],
    
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cplex",
    "MIP_feasiblity_emphasis" => true
                                
    )

    return params
    
end

function decompose_controlled_H_with_R()

    println(">>>>> Controlled-H with R Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 5,    

    "elementary_gates" => ["RY", "CNot_12", "Identity"], 
    "target_gate" => QCO.CHGate(),
       
    "RX_discretization" => [], 
    "RY_discretization" => [-π/4, π/4, π/2, -π/2], 
    "RZ_discretization" => [], 
  
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cplex"
                                
    )

    return params
    
end

function decompose_magic_M()
    
    println(">>>>> M Gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 5,    
    
        "elementary_gates" => ["U3", "CNot_12", "CNot_21", "Identity"], 
        "target_gate" => QCO.MGate(),   
           
        "U_θ_discretization" => [0, π/2],
        "U_ϕ_discretization" => [-π/2, π/2],
        "U_λ_discretization" => [-π/2, π],
          
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",
        
        "optimizer" => "cplex"
                                    
        )
    
        return params
    
end

function decompose_magic_M_using_SHCnot()
    
    println(">>>>> M gate using S, H and CNOT Gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 5,    
    
        "elementary_gates" => ["S_1", "S_2", "H_1", "H_2", "CNot_12", "CNot_21", "Identity"], 
        "target_gate" => QCO.MGate(),
           
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",
           
        "optimizer" => "cplex"
                                    
        )
    
        return params
    
end

function decompose_S()

    println(">>>>> S Gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 3,    
    
        "elementary_gates" => ["U3", "CNot_12", "Identity"], 
        "target_gate" => QCO.kron_single_gate(2, QCO.SGate(), "q1"),
           
        "U_θ_discretization" => [0, π/4, π/2],
        "U_ϕ_discretization" => [-π/2, 0, π/2],
        "U_λ_discretization" => [0, π],
         
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",
    
        "optimizer" => "cplex",
                                    
        )
    
        return params
    
end

function decompose_cnot_21()

    println(">>>>> CNOT_21 Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 5,    

    "elementary_gates" => ["H_1", "H_2", "Identity", "CNot_12"],  
    "target_gate" => QCO.CNotRevGate(),
 
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cplex",
                                
    )

    return params
    
end

function decompose_cnot_21_with_U()
    
    println(">>>>> CNOT_21 using U3 and CNOT Gates <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 6,    
    
        "elementary_gates" => ["U3", "CNot_12", "Identity"], 
        "target_gate" => QCO.CNotRevGate(),   
           
        "U_θ_discretization" => [-π/2, π/2],
        "U_ϕ_discretization" => [0, π/2],
        "U_λ_discretization" => [0],
          
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact", 
        
        "optimizer" => "cplex",
                              
        )
    
        return params
    
end

function decompose_swap()

    println(">>>>> SWAP Gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 5,    
    
        "elementary_gates" => ["CNot_21", "CNot_12", "Identity"], 
        "target_gate" => QCO.SwapGate(),   
           
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact", 
        
        "optimizer" => "cplex",
                                    
        )
    
        return params
    
end

function decompose_W()

    println(">>>>> W Gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 5,    
    
        "elementary_gates" => ["U3", "CNot_21", "CNot_12", "Identity"], 
        "target_gate" => QCO.WGate(),   

        "U_θ_discretization" => [-π/4, π/4],
        "U_ϕ_discretization" => [0],
        "U_λ_discretization" => [0],
                   
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact", 
        
        "optimizer" => "cplex",
        
        )
    
        return params
    
end

function decompose_W_using_HCnot()

    println(">>>>> W Gate using H and CNOT gates <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 6,    
    
        "elementary_gates" => ["CH_12", "CNot_21", "CNot_12", "Identity"], 
        "target_gate" => QCO.WGate(),   
                  
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact", 
        
        "optimizer" => "cplex",
                                    
        )
    
        return params
    
end

function decompose_HCoinGate()

    println(">>>>> Hadamard Coin gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 14,    
    
        "elementary_gates" => ["Y_1", "Y_2", "Z_1", "Z_2", "T_2", "Tdagger_1", "Sdagger_1", "SX_1", "SXdagger_2", "CNot_21", "CNot_12", "Identity"], 
        "target_gate" => -QCO.HCoinGate(),   
                  
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact", 
        
        "optimizer" => "cplex",
        # "MIP_feasiblity_emphasis" => true
                                    
        )
    
        return params
    
end

function decompose_GroverDiffusionGate()

    println(">>>>> Grover's Diffusion Operator <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 10,    
    
        "elementary_gates" => ["X_1", "X_2", "H_1", "H_2", "CNot_12", "Identity"], 
        "target_gate" => QCO.GroverDiffusionGate(),   
                  
        "objective" => "minimize_depth",
        "decomposition_type" => "exact", 
        
        "optimizer" => "cplex"
                                    
        )
    
        return params
    
end
