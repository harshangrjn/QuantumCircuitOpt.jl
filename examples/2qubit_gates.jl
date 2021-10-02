function decompose_hadamard()

    println(">>>>> Hadamard Gate <<<<<")
 
    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 3,    

    "elementary_gates" => ["U3_1", "CNot_1_2", "Identity"], 
    "target_gate" => QCO.kron_single_qubit_gate(2, QCO.HGate(), "q1"),
       
    "U3_θ_discretization" => [0, π/2],
    "U3_ϕ_discretization" => [0, π/2],
    "U3_λ_discretization" => [0, π],
 
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

    "elementary_gates" => ["U3_2", "CNot_1_2", "Identity"], 
    "target_gate" => QCO.CZGate(),

    "U3_θ_discretization" => -π/2:π/2:π/2,
    "U3_ϕ_discretization" => -π/2:π/2:π/2,
    "U3_λ_discretization" => -π/2:π/2:π/2,

    "set_cnot_lower_bound" => 1,

    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",

    "optimizer" => "gurobi",
    "optimizer_presolve" => false
                                
    )

    return params
    
end

function decompose_controlled_V()

    println(">>>>> Controlled-V Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 7,

    "elementary_gates" => ["H_1", "H_2", "T_1", "T_2", "Tdagger_1", "CNot_1_2", "CNot_2_1", "Identity"],
    "target_gate" => QCO.CVGate(),
    
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",

    "set_cnot_lower_bound" => 2,
    
    "optimizer" => "gurobi"
                                
    )

    return params
    
end

function decompose_controlled_H()

    println(">>>>> Controlled-H Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 5,    

    "elementary_gates" => ["U3_2", "CNot_1_2", "Identity"], 
    "target_gate" => QCO.CHGate(),

    "U3_θ_discretization" => [-π/4, 0, π/4],
    "U3_ϕ_discretization" => [0],
    "U3_λ_discretization" => [0],
    
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

    "elementary_gates" => ["RY_1", "RY_2", "CNot_1_2", "Identity"], 
    "target_gate" => QCO.CHGate(),
       
    "RX_discretization" => [], 
    "RY_discretization" => [-π/4, π/4, π/2, -π/2], 
    "RZ_discretization" => [], 

    "set_cnot_lower_bound" => 1,
    "set_cnot_upper_bound" => 2,
  
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    "identify_real_gates" => true,
    
    "optimizer" => "gurobi"
                                
    )

    return params
    
end

function decompose_magic_M()
    
    println(">>>>> M Gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 6,    
    
        "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "CNot_2_1", "Identity"], 
        "target_gate" => QCO.MGate(),   
           
        "U3_θ_discretization" => [0, π/2],
        "U3_ϕ_discretization" => [-π/2, π/2],
        "U3_λ_discretization" => [-π/2, π],
          
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",

        "set_cnot_lower_bound" => 1,
        
        "optimizer" => "gurobi"
                                    
        )
    
        return params
    
end

function decompose_magic_M_using_SHCnot()
    
    println(">>>>> M gate using S, H and CNOT Gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 5,    
    
        "elementary_gates" => ["S_1", "S_2", "H_1", "H_2", "CNot_1_2", "CNot_2_1", "Identity"], 
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
    
        "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "Identity"], 
        "target_gate" => QCO.kron_single_qubit_gate(2, QCO.SGate(), "q1"),
           
        "U3_θ_discretization" => [0, π/4, π/2],
        "U3_ϕ_discretization" => [-π/2, 0, π/2],
        "U3_λ_discretization" => [0, π],
         
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",
    
        "optimizer" => "cplex",
                                    
        )
    
        return params
    
end

function decompose_revcnot()

    println(">>>>> Reverse CNOT Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 5,    

    "elementary_gates" => ["H_1", "H_2", "Identity", "CNot_1_2"],  
    "target_gate" => QCO.CNotRevGate(),
 
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cplex",
                                
    )

    return params
    
end

function decompose_revcnot_with_U()
    
    println(">>>>> Reverse CNOT using U3 and CNOT Gates <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 6,    
    
        "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "Identity"], 
        "target_gate" => QCO.CNotRevGate(),   
           
        "U3_θ_discretization" => [-π/2, π/2],
        "U3_ϕ_discretization" => [0, π/2],
        "U3_λ_discretization" => [0],
          
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
    
        "elementary_gates" => ["CNot_2_1", "CNot_1_2", "Identity"], 
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
    
        "elementary_gates" => ["U3_2", "CNot_2_1", "CNot_1_2", "Identity"], 
        "target_gate" => QCO.WGate(),

        "U3_θ_discretization" => [-π/4, 0, π/4],
        "U3_ϕ_discretization" => [0, π/4],
        "U3_λ_discretization" => [0, π/2],
                   
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
    
        "elementary_gates" => ["CH_1_2", "CNot_2_1", "CNot_1_2", "Identity"], 
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
    
        "elementary_gates" => ["Y_1", "Y_2", "Z_1", "Z_2", "T_2", "Tdagger_1", "Sdagger_1", "SX_1", "SXdagger_2", "CNot_2_1", "CNot_1_2", "Identity"], 
        "target_gate" => QCO.HCoinGate(),   
                  
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact", 
        "identify_real_gates" => true,
        
        "optimizer" => "cplex",
        
                                    
        )
    
        return params
    
end

function decompose_GroverDiffusionGate()

    println(">>>>> Grover's Diffusion Operator <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 10,    
    
        "elementary_gates" => ["X_1", "X_1xX_2", "H_1xH_2", "X_2", "H_1", "H_2", "CNot_1_2", "Identity"], 
        "target_gate" => QCO.GroverDiffusionGate(),   
                  
        "objective" => "minimize_depth",
        "decomposition_type" => "exact", 
        "identify_real_gates" => true,
        
        "optimizer" => "cplex"
     
        )
    
        return params
    
end

function decompose_iSwapGate()

    println(">>>>> iSwap Gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 10,    
    
        "elementary_gates" => ["T_1", "T_2", "H_1", "H_2", "CNot_1_2", "CNot_2_1", "Identity"],
        "target_gate" => QCO.iSwapGate(),   
                  
        "objective" => "minimize_depth",
        "decomposition_type" => "exact", 
        
        "optimizer" => "cplex"
                                    
        )
    
        return params
    
end

