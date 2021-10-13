#=
Settings for results in the paper "QuantumCircuitOpt: An Open-source Frameworkfor Provably Optimal 
Quantum Circuit Design"

Version of QCOpt: v0.3.0

Last updated: Oct 12, 2021
=#

function decompose_controlled_Z()

    println(">>>>> Controlled-Z Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 4,

    "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "Identity"], 
    "target_gate" => QCOpt.CZGate(),

    "U3_θ_discretization" => -π:π/2:π,
    "U3_ϕ_discretization" => -π:π/2:π,
    "U3_λ_discretization" => -π:π/2:π,

    "objective" => "minimize_depth", 
    "decomposition_type" => "exact")

    return params
    
end

function decompose_controlled_V()

    println(">>>>> Controlled-V Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 7,

    "elementary_gates" => ["H_1", "H_2", "T_1", "T_2", "Tdagger_1", "Tdagger_2", "CNot_1_2", "CNot_2_1", "Identity"],
    "target_gate" => QCOpt.CVGate(),
    
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact")

    return params
    
end

function decompose_controlled_H()

    println(">>>>> Controlled-H Gate <<<<<")

    params = Dict{String, Any}(
    
    "num_qubits" => 2, 
    "depth" => 5,    

    "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "CNot_2_1", "Identity"], 
    "target_gate" => QCOpt.CHGate(),

    "U3_θ_discretization" => -2*π:π/4:2*π,
    "U3_ϕ_discretization" => [0],
    "U3_λ_discretization" => [0],
    
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact")

    return params
    
end

function decompose_magic_using_U3_CNot_2_1()
    
    println(">>>>> Magic basis using U3 and CNot_2_1 <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 5,    
    
        "elementary_gates" => ["U3_1", "U3_2", "CNot_2_1", "Identity"], 
        "target_gate" => QCOpt.MGate(),   
        
        "U3_θ_discretization" => -π:π/2:π,
        "U3_ϕ_discretization" => -π:π/2:π,
        "U3_λ_discretization" => -π:π/2:π,
          
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact")
    
        return params 
end

function decompose_iSwapGate()

    println(">>>>> iSwap Gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 10,    
    
        "elementary_gates" => ["T_1", "T_2", "Tdagger_1", "Tdagger_2", "H_1", "H_2", "CNot_1_2", "CNot_2_1", "Identity"],
        "target_gate" => QCOpt.iSwapGate(),   
                  
        "objective" => "minimize_depth",
        "decomposition_type" => "exact")
    
        return params
    
end

function decompose_GroverDiffusion_using_Clifford()

    println(">>>>> Grover's Diffusion Operator <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2,
        "depth" => 6,    
    
        "elementary_gates" => ["X_1", "X_2", "H_1", "H_2", "S_1", "S_2", "T_1", "T_2", "Y_1", "Y_2", "CNot_1_2", "Identity"], 
        
        "target_gate" => QCOpt.GroverDiffusionGate(),
                  
        "objective" => "minimize_depth",
        "decomposition_type" => "exact")
    
        return params
    
end


function decompose_GroverDiffusion_using_U3()

    println(">>>>> Grover's Diffusion Operator using U3 gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 10,    
    
        "elementary_gates" => ["U3_1", "CNot_1_2", "Identity"], 

        "U3_θ_discretization" => -π:π/2:π,
        "U3_ϕ_discretization" => -π:π/2:π,
        "U3_λ_discretization" => [0],
        
        "target_gate" => QCOpt.GroverDiffusionGate(),
                  
        "objective" => "minimize_depth",
        "decomposition_type" => "exact", 
        "identify_real_gates" => true)
    
        return params
    
end

function decompose_magic_using_SH_CNot_1_2()
    
    println(">>>>> M gate using S, H and CNOT_1_2 Gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 10,    
    
        "elementary_gates" => ["S_1", "S_2", "H_1", "H_2", "CNot_1_2", "Identity"], 
        "target_gate" => QCOpt.MGate(),
           
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact")
    
        return params 
end


function decompose_magic_using_SH_CNot_2_1()
    
    println(">>>>> M gate using S, H and CNOT_2_1 Gate <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 10,    
    
        "elementary_gates" => ["S_1", "S_2", "H_1", "H_2", "CNot_2_1", "Identity"], 
        "target_gate" => QCOpt.MGate(),
           
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact")
    
        return params 
end

function decompose_magic_using_U3_CNot_1_2()
    
    println(">>>>> Magic basis using U3 and CNot_1_2 <<<<<")

    params = Dict{String, Any}(
    
        "num_qubits" => 2, 
        "depth" => 5,    
    
        "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "Identity"], 
        "target_gate" => QCOpt.MGate(),   
        
        "U3_θ_discretization" => -π:π/2:π,
        "U3_ϕ_discretization" => -π:π:π,
        "U3_λ_discretization" => -π:π/2:π,
          
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact")
    
        return params 
end

function decompose_toffoli_with_controlled_gates()

    # Reference: https://doi.org/10.1109/TCAD.2005.858352
    println(">>>>> Toffoli gate using controlled gates <<<<<")
 
    params = Dict{String, Any}(
    
    "num_qubits" => 3,
    "depth" => 5,

    "elementary_gates" => ["CV_1_3", "CV_2_3", "CV_1_2", "CVdagger_1_3", "CVdagger_2_3", "CVdagger_1_2", "CNot_2_1", "CNot_1_2", "Identity"],

    "target_gate" => QCOpt.ToffoliGate(),
    
    "objective" => "minimize_depth",
    "decomposition_type" => "exact"
    )

    return params
    
end

function decompose_Fredkin()

    println(">>>>> Fredkin gate using controlled gates <<<<<")

    params = Dict{String, Any}(
    "num_qubits" => 3,
    "depth" => 7,

    # Reference: https://doi.org/10.1103/PhysRevA.53.2855
    "elementary_gates" => ["CV_1_2", "CV_2_3", "CV_1_3", "CVdagger_1_2", "CVdagger_2_3", "CVdagger_1_3", "CNot_1_2", "CNot_3_2", "CNot_2_3", "CNot_1_3", "Identity"],
    "target_gate" => QCOpt.CSwapGate(), #also Fredkin

    "objective" => "minimize_depth", 
    "decomposition_type" => "exact"
    )
    
    return params
end

function decompose_double_toffoli()

    println(">>>>> Double Toffoli using controlled gates <<<<<")

    # Reference: https://doi.org/10.1109/TCAD.2005.858352 
    
    num_qubits = 4

    function target_gate()
        CV_3_4 = QCOpt.get_full_sized_gate("CV_3_4", num_qubits);
        CV_2_4 = QCOpt.get_full_sized_gate("CV_2_4", num_qubits);
        CVdagger_2_4 = QCOpt.get_full_sized_gate("CVdagger_2_4", num_qubits);        
        CNot_1_3 = QCOpt.get_full_sized_gate("CNot_1_3", num_qubits);
        CNot_3_2 = QCOpt.get_full_sized_gate("CNot_3_2", num_qubits);

        return CV_2_4 * CNot_1_3 * CNot_3_2 * CV_3_4 * CVdagger_2_4 * CNot_3_2 * CNot_1_3
    end

    params = Dict{String, Any}(
    "num_qubits" => 4,
    "depth" => 7,

    "elementary_gates" => ["CV_1_2", "CV_2_4", "CV_3_4", "CVdagger_1_2", "CVdagger_2_4", "CVdagger_3_4", "CNot_1_3", "CNot_3_2", "CNot_2_3", "Identity"],
    "target_gate" => target_gate(),

    "objective" => "minimize_depth",
    "decomposition_type" => "exact"
    )
    
    return params
end

function decompose_quantum_fulladder()
    #Reference-1: https://doi.org/10.1109/DATE.2005.249
    #Reference-2: https://doi.org/10.1109/TCAD.2005.858352 
    
    println(">>>>> Quantum Full adder using controlled gates <<<<<")

    num_qubits = 4

    function target_gate()

        CV_1_2 = QCOpt.get_full_sized_gate("CV_1_2", num_qubits);
        CV_4_2 = QCOpt.get_full_sized_gate("CV_4_2", num_qubits);
        CVdagger_1_2 = QCOpt.get_full_sized_gate("CVdagger_1_2", num_qubits);
        CVdagger_3_2 = QCOpt.get_full_sized_gate("CVdagger_3_2", num_qubits);
        CNot_3_1 = QCOpt.get_full_sized_gate("CNot_3_1", num_qubits);
        CNot_4_3 = QCOpt.get_full_sized_gate("CNot_4_3", num_qubits);
        CNot_2_4 = QCOpt.get_full_sized_gate("CNot_2_4", num_qubits);

        return CNot_2_4 * CV_4_2 * CVdagger_1_2 * CVdagger_3_2 * CNot_4_3 * CNot_3_1 * CV_1_2
    end

    params = Dict{String, Any}(
        "num_qubits" => num_qubits,
        "depth" => 7,
    
        "elementary_gates" => ["CV_1_2", "CV_4_2", "CV_3_2", "CVdagger_1_2", "CVdagger_4_2", "CVdagger_3_2", "CNot_3_1", "CNot_4_3", "CNot_2_4", "CNot_4_1", "Identity"],
        
        "target_gate" => target_gate(),
        "identify_real_gates" => true,
    
        "objective" => "minimize_depth",  
        "optimizer" => "gurobi"
        
        )
    
        return params
end