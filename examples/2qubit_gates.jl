function test_hadamard()
    println(">>>>> Decomposing Hadamard Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 3,    

    "elementary_gates" => ["U3", "cnot_12", "Identity"], 
    "target_gate" => "H1",
       
    "U_θ_discretization" => [0, π/2],
    "U_ϕ_discretization" => [0, π/2],
    "U_λ_discretization" => [0, π],
 
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cplex",
    
    )

    return params
    
end

function test_controlled_Z()
    println(">>>>> Decomposing Controlled-Z Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 4,    

    "elementary_gates" => ["U3", "cnot_12", "Identity"], 
    "target_gate" => "controlled_Z",
       
    "U_θ_discretization" => [-π/2, 0, π/2],
    "U_ϕ_discretization" => [0, π/2],
    "U_λ_discretization" => [0, π/2],

    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cplex",
                                
    )

    return params
    
end

function test_controlled_V()

    println(">>>>> Decomposing Controlled-V Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 7,    

    "elementary_gates" => ["H1", "H2", "T1", "T2", "T1_conjugate", "cnot_12", "cnot_21"],
    "target_gate" => "controlled_V",
    
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cplex",
                                
    )

    return params
    
end

function test_controlled_H()

    println(">>>>> Decomposing Controlled-H Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 5,    

    "elementary_gates" => ["U3", "cnot_12", "Identity"], 
    "target_gate" => "controlled_H_12",

    "U_θ_discretization" => [-π/4, 0, π/4],
    "U_ϕ_discretization" => [0],
    "U_λ_discretization" => [0],
    
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cplex",
                                
    )

    return params
    
end

function test_controlled_H_with_R()

    println(">>>>> Decomposing Controlled-H with R Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 5,    

    "elementary_gates" => ["R_y", "cnot_12", "Identity"], 
    "target_gate" => "controlled_H_12",
       
    "R_x_discretization" => [], 
    "R_y_discretization" => [-π/4, π/4, π/2, -π/2], 
    "R_z_discretization" => [], 
  
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cplex",
                                
    )

    return params
    
end

function test_controlled_R2()

    println(">>>>> Decomposing Controlled-R2 Gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 5,    
    
        "elementary_gates" => ["U3", "cnot_12", "Identity"], 
        "target_gate" => "controlled_R2",
           
        "U_θ_discretization" => [0],
        "U_ϕ_discretization" => [0],
        "U_λ_discretization" => [-π/4, π/4],
           
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",  
        
        "optimizer" => "cplex",
                                    
        )
    
        return params
    
end

function test_magic_M()
    
    println(">>>>> Decomposing M Gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 4,    
    
        "elementary_gates" => ["U3", "cnot_12", "cnot_21", "Identity"], 
        "target_gate" => "magic_M",   
           
        "U_θ_discretization" => [0, π/2],
        "U_ϕ_discretization" => [-π/2, π/2],
        "U_λ_discretization" => [-π/2, π],
          
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",
        
        "optimizer" => "cplex",
                                    
        )
    
        return params
    
end

function test_magic_M_using_SHCnot()
    
    println(">>>>> Decomposing M gate using S, H and CNOT Gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 5,    
    
        "elementary_gates" => ["S1", "S2", "H1", "H2", "cnot_12", "cnot_21", "Identity"], 
        "target_gate" => "magic_M",
           
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",
           
        "optimizer" => "cplex",
                                    
        )
    
        return params
    
end

function test_S()

    println(">>>>> Decomposing S Gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 3,    
    
        "elementary_gates" => ["U3", "cnot_12", "Identity"], 
        "target_gate" => "S1",
           
        "U_θ_discretization" => [0, π/2],
        "U_ϕ_discretization" => [0, π/2],
        "U_λ_discretization" => [0, π],
         
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",
    
        "optimizer" => "cplex",
                                    
        )
    
        return params
    
end

function test_cnot_21()

    println(">>>>> Decomposing CNOT_21 Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 5,    

    "elementary_gates" => ["H1", "H2", "Identity", "cnot_12"],  
    "target_gate" => "cnot_21",
 
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cplex",
                                
    )

    return params
    
end

function test_cnot_21_with_U()

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 5,    
    
        "elementary_gates" => ["U3", "cnot_12", "Identity"], 
        "target_gate" => "cnot_21",   
           
        "U_θ_discretization" => [-π/2, π/2],
        "U_ϕ_discretization" => [0, π/2],
        "U_λ_discretization" => [0],
          
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact", 
        
        "optimizer" => "cplex",
                              
        )
    
        return params
    
end

function test_swap()

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 5,    
    
        "elementary_gates" => ["cnot_21", "cnot_12", "Identity"], 
        "target_gate" => "swap",   
           
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact", 
        
        "optimizer" => "cplex",
                                    
        )
    
        return params
    
end

function test_W()

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 5,    
    
        "elementary_gates" => ["U3", "cnot_21", "cnot_12", "Identity"], 
        "target_gate" => "W_hermitian",   

        "U_θ_discretization" => [-π/4, π/4],
        "U_ϕ_discretization" => [0],
        "U_λ_discretization" => [0],
                   
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact", 
        
        "optimizer" => "cplex",
        
        )
    
        return params
    
end

function test_W_using_HCnot()

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 6,    
    
        "elementary_gates" => ["controlled_H_12", "cnot_21", "cnot_12", "Identity"], 
        "target_gate" => "W_hermitian",   
                  
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact", 
        
        "optimizer" => "cplex",
                                    
        )
    
        return params
    
end