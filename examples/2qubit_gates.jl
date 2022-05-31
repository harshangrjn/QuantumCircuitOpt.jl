function decompose_hadamard()

    println(">>>>> Hadamard Gate <<<<<")
 
    return Dict{String, Any}(
    
    "num_qubits" => 2, 
    "maximum_depth" => 3,    
    "elementary_gates" => ["U3_1", "CNot_1_2", "Identity"], 
    "target_gate" => QCOpt.get_full_sized_gate("H_1", 2),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
       
    "U3_θ_discretization" => [0, π/2],
    "U3_ϕ_discretization" => [0, π/2],
    "U3_λ_discretization" => [0, π],
    )
 
end

function decompose_controlled_Z()

    println(">>>>> Controlled-Z Gate <<<<<")

    return Dict{String, Any}(
    
    "num_qubits" => 2, 
    "maximum_depth" => 4,
    "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "Identity"], 
    "target_gate" => QCOpt.CZGate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",

    "U3_θ_discretization" => -π:π/2:π,
    "U3_ϕ_discretization" => -π:π/2:π,
    "U3_λ_discretization" => -π:π/2:π,
    )

end

function decompose_controlled_V()

    println(">>>>> Controlled-V Gate <<<<<")

    return Dict{String, Any}(
    
    "num_qubits" => 2, 
    "maximum_depth" => 7,
    "elementary_gates" => ["H_1", "H_2", "T_1", "T_2", "Tdagger_1", "Tdagger_2", "CNot_1_2", "CNot_2_1", "Identity"],
    "target_gate" => QCOpt.CVGate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    )

end

function decompose_controlled_H()

    println(">>>>> Controlled-H Gate <<<<<")

    return Dict{String, Any}(
    
    "num_qubits" => 2, 
    "maximum_depth" => 5,
    "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "CNot_2_1", "Identity"], 
    "target_gate" => QCOpt.CHGate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",

    "U3_θ_discretization" => -2*π:π/4:2*π,
    "U3_ϕ_discretization" => [0],
    "U3_λ_discretization" => [0],
    )

end

function decompose_controlled_H_with_R()

    println(">>>>> Controlled-H with R Gate <<<<<")

    return Dict{String, Any}(
    
    "num_qubits" => 2, 
    "maximum_depth" => 5,
    "elementary_gates" => ["RY_1", "RY_2", "CNot_1_2", "Identity"], 
    "target_gate" => QCOpt.CHGate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
       
    "RY_discretization" => [-π/4, π/4, π/2, -π/2],
    "set_cnot_lower_bound" => 1,
    "set_cnot_upper_bound" => 2,           
    )

end

function decompose_magic()
    
    println(">>>>> M Gate <<<<<")

    return Dict{String, Any}(
    
        "num_qubits" => 2, 
        "maximum_depth" => 4,
        "elementary_gates" => ["U3_1", "U3_2", "CNot_2_1", "CNot_1_2", "Identity"], 
        "target_gate" => QCOpt.MGate(),   
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",

        "U3_θ_discretization" => -π:π/2:π,
        "U3_ϕ_discretization" => -π:π/2:π,
        "U3_λ_discretization" => -π:π/2:π,                         
        )
     
end

function decompose_magic_using_CNOT_1_2()
    
    println(">>>>> M Gate <<<<<")

    return Dict{String, Any}(
    
        "num_qubits" => 2, 
        "maximum_depth" => 5,
        "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "Identity"], 
        "target_gate" => QCOpt.MGate(),
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",
           
        "U3_θ_discretization" => -π:π/2:π,
        "U3_ϕ_discretization" => -π:π:π,
        "U3_λ_discretization" => -π:π/2:π,            
        )
    
end

function decompose_magic_using_SHCnot()
    
    println(">>>>> M gate using S, H and CNOT Gate <<<<<")

    return Dict{String, Any}(
    
        "num_qubits" => 2, 
        "maximum_depth" => 5,
        "elementary_gates" => ["S_1", "S_2", "H_1", "H_2", "CNot_1_2", "CNot_2_1", "Identity"], 
        "target_gate" => QCOpt.MGate(),
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",
        )
    
end

function decompose_S()

    println(">>>>> S Gate <<<<<")

    return Dict{String, Any}(
    
        "num_qubits" => 2, 
        "maximum_depth" => 3,
        "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "Identity"], 
        "target_gate" => QCOpt.get_full_sized_gate("S_1", 2),
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",
           
        "U3_θ_discretization" => [0, π/4, π/2],
        "U3_ϕ_discretization" => [-π/2, 0, π/2],
        "U3_λ_discretization" => [0, π],           
        )
end

function decompose_revcnot()

    println(">>>>> Reverse CNOT Gate <<<<<")

    return Dict{String, Any}(
    
    "num_qubits" => 2, 
    "maximum_depth" => 5,
    "elementary_gates" => ["H_1", "H_2", "CNot_1_2", "Identity"],  
    "target_gate" => QCOpt.CNotRevGate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",       
    )

end

function decompose_revcnot_with_U()
    
    println(">>>>> Reverse CNOT using U3 and CNOT Gates <<<<<")

    return Dict{String, Any}(
    
        "num_qubits" => 2, 
        "maximum_depth" => 5,    
        "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "Identity"], 
        "target_gate" => QCOpt.CNotRevGate(),
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",

        "U3_θ_discretization" => -π:π/4:π,
        "U3_ϕ_discretization" => [0],
        "U3_λ_discretization" => [0],
        )
    
end

function decompose_swap()

    println(">>>>> SWAP Gate <<<<<")

    return Dict{String, Any}(
    
        "num_qubits" => 2, 
        "maximum_depth" => 5,    
        "elementary_gates" => ["CNot_2_1", "CNot_1_2", "Identity"], 
        "target_gate" => QCOpt.SwapGate(),
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",           
        )
    
end

function decompose_W()

    println(">>>>> W Gate <<<<<")

    return Dict{String, Any}(
    
        "num_qubits" => 2, 
        "maximum_depth" => 5,
        "elementary_gates" => ["U3_2", "CNot_2_1", "CNot_1_2", "Identity"], 
        "target_gate" => QCOpt.WGate(),
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",

        "U3_θ_discretization" => [-π/4, 0, π/4],
        "U3_ϕ_discretization" => [0, π/4],
        "U3_λ_discretization" => [0, π/2],
        )
    
end

function decompose_W_using_HCnot()

    println(">>>>> W Gate using H and CNOT gates <<<<<")

    return Dict{String, Any}(
    
        "num_qubits" => 2, 
        "maximum_depth" => 6,
        "elementary_gates" => ["CH_1_2", "CNot_2_1", "CNot_1_2", "Identity"], 
        "target_gate" => QCOpt.WGate(),
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",
        )
    
end

function decompose_HCoinGate()

    println(">>>>> Hadamard Coin gate <<<<<")

    return Dict{String, Any}(
    
        "num_qubits" => 2, 
        "maximum_depth" => 10,
        "elementary_gates" => ["Y_1", "Y_2", "Z_1", "Z_2", "T_2", "Tdagger_1", "Sdagger_1", "SX_1", "SXdagger_2", "CNot_2_1", "CNot_1_2", "Identity"], 
        "target_gate" => QCOpt.HCoinGate(),  
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",
        )
    
end

function decompose_GroverDiffusion_using_HX()

    println(">>>>> Grover's Diffusion Operator <<<<<")

    return Dict{String, Any}(
        
        #Reference: https://arxiv.org/pdf/1804.03719.pdf (Fig 6)

        "num_qubits" => 2, 
        "maximum_depth" => 10,
        "elementary_gates" => ["X_1", "X_1xX_2", "H_1xH_2", "X_2", "H_1", "H_2", "CNot_1_2", "Identity"],
        "target_gate" => QCOpt.GroverDiffusionGate(),          
        "objective" => "minimize_depth",
        "decomposition_type" => "exact",
        )
    
end

function decompose_GroverDiffusion_using_Clifford()

    println(">>>>> Grover's Diffusion Operator <<<<<")

    return Dict{String, Any}(
    
        "num_qubits" => 2,
        "maximum_depth" => 6,
        "elementary_gates" => ["X_1", "X_2", "H_1", "H_2", "S_1", "S_2", "T_1", "T_2", "Y_1", "Y_2", "CNot_1_2", "Identity"], 
        "target_gate" => QCOpt.GroverDiffusionGate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "exact",
        )
    
end


function decompose_GroverDiffusion_using_U3()

    println(">>>>> Grover's Diffusion Operator using U3 gate <<<<<")

    return Dict{String, Any}(
    
        "num_qubits" => 2, 
        "maximum_depth" => 10,
        "elementary_gates" => ["U3_1", "CNot_1_2", "Identity"], 
        "target_gate" => QCOpt.GroverDiffusionGate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "exact",

        "U3_θ_discretization" => -π:π/2:π,
        "U3_ϕ_discretization" => -π:π/2:π,
        "U3_λ_discretization" => [0],
        )

end

function decompose_iSwap()

    println(">>>>> iSwap Gate <<<<<")

    return Dict{String, Any}(
    
        "num_qubits" => 2, 
        "maximum_depth" => 10,
        "elementary_gates" => ["T_1", "T_2", "Tdagger_1", "Tdagger_2", "H_1", "H_2", "CNot_1_2", "CNot_2_1", "Identity"],
        # "elementary_gates" => ["S_1", "S_2", "Sdagger_1", "Sdagger_2", "H_1", "H_2", "CNot_1_2", "CNot_2_1", "Identity"],
        "target_gate" => QCOpt.iSwapGate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "exact",                           
        )

end

function decompose_qft2_using_R()
    println(">>>>> QFT2 Gate using R and CNOT gates <<<<<")

    return Dict{String, Any}(
    
        "num_qubits" => 2, 
        "maximum_depth" => 10,
        "elementary_gates" => ["RX_1", "RZ_1", "RZ_2", "CNot_1_2", "CNot_2_1", "Identity"],
        "target_gate" => QCOpt.QFT2Gate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "approximate",
        
        "RX_discretization" => [π/2],
        "RZ_discretization" => [-π/4, π/2, 3*π/4, 7*π/4],                      
        )
    
end

function decompose_qft2_using_HT()
    println(">>>>> QFT2 Gate using H, T, CNOT gates <<<<<")

    return Dict{String, Any}(
    
        "num_qubits" => 2, 
        "maximum_depth" => 8,
        "elementary_gates" => ["H_1", "H_2", "T_1", "T_2", "Tdagger_1", "Tdagger_2", "CNot_1_2", "CNot_2_1", "Identity"],
        "target_gate" => QCOpt.QFT2Gate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "exact", 
        )
    
end

function decompose_RGate()
    println(">>>>> RGate testing <<<<<")
    
    R1 = qco.get_full_sized_gate("R_1", 2; angle = [π/6, π/3])
    R2 = qco.get_full_sized_gate("R_2", 2; angle = [π/3, π/6])
    CNot_1_2 = QCOpt.get_full_sized_gate("CNot_1_2", 2)
    T = QCOpt.round_complex_values(R2 * CNot_1_2 * R1)

    return Dict{String, Any}(
    
        "num_qubits" => 2, 
        "maximum_depth" => 3,
        "elementary_gates" => ["R_1", "R_2", "CNot_1_2", "Identity"],
        "R_θ_discretization" => [π/3, π/6],
        "R_ϕ_discretization" => [π/3, π/6],
        "target_gate" => T,
        "objective" => "minimize_depth",
        "decomposition_type" => "exact", 
        )
    
end