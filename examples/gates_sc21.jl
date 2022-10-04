#=
Settings for results in the paper "QuantumCircuitOpt: An Open-source Frameworkfor Provably Optimal 
Quantum Circuit Design"

Version of QCOpt: v0.3.0

Last updated: Oct 12, 2021
=#

function controlled_Z()

    println(">>>>> Controlled-Z Gate <<<<<")

    return Dict{String, Any}(
    "num_qubits" => 2, 
    "maximum_depth" => 4,
    "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "Identity"], 
    "target_gate" => QCOpt.CZGate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal",

    "U3_θ_discretization" => -π:π/2:π,
    "U3_ϕ_discretization" => -π:π/2:π,
    "U3_λ_discretization" => -π:π/2:π,)    
end

function controlled_V()

    println(">>>>> Controlled-V Gate <<<<<")

    return Dict{String, Any}(
    "num_qubits" => 2, 
    "maximum_depth" => 7,
    "elementary_gates" => ["H_1", "H_2", "T_1", "T_2", "Tdagger_1", "Tdagger_2", "CNot_1_2", "CNot_2_1", "Identity"],
    "target_gate" => QCOpt.CVGate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal")
end

function controlled_H()

    println(">>>>> Controlled-H Gate <<<<<")

    return Dict{String, Any}(
    "num_qubits" => 2, 
    "maximum_depth" => 5,    
    "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "CNot_2_1", "Identity"], 
    "target_gate" => QCOpt.CHGate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal",

    "U3_θ_discretization" => -2*π:π/4:2*π,
    "U3_ϕ_discretization" => [0],
    "U3_λ_discretization" => [0],
    )
end

function magic_using_U3_CNot_2_1()
    
    println(">>>>> Magic basis using U3 and CNot_2_1 <<<<<")

    return Dict{String, Any}(
        "num_qubits" => 2, 
        "maximum_depth" => 5,
        "elementary_gates" => ["U3_1", "U3_2", "CNot_2_1", "Identity"], 
        "target_gate" => QCOpt.MGate(),
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact_optimal",
        
        "U3_θ_discretization" => -π:π/2:π,
        "U3_ϕ_discretization" => -π:π/2:π,
        "U3_λ_discretization" => -π:π/2:π,
    )     
end

function iSwapGate()

    println(">>>>> iSwap Gate <<<<<")

    return Dict{String, Any}(
        "num_qubits" => 2, 
        "maximum_depth" => 10,
        "elementary_gates" => ["T_1", "T_2", "Tdagger_1", "Tdagger_2", "H_1", "H_2", "CNot_1_2", "CNot_2_1", "Identity"],
        "target_gate" => QCOpt.iSwapGate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "exact_optimal")
    
        
    
end

function GroverDiffusion_using_Clifford()

    println(">>>>> Grover's Diffusion Operator <<<<<")

    return Dict{String, Any}(
        "num_qubits" => 2,
        "maximum_depth" => 6,
        "elementary_gates" => ["X_1", "X_2", "H_1", "H_2", "S_1", "S_2", "T_1", "T_2", "Y_1", "Y_2", "CNot_1_2", "Identity"], 
        "target_gate" => QCOpt.GroverDiffusionGate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "exact_optimal")
end


function GroverDiffusion_using_U3()

    println(">>>>> Grover's Diffusion Operator using U3 gate <<<<<")

    return Dict{String, Any}(
        "num_qubits" => 2, 
        "maximum_depth" => 10,
        "elementary_gates" => ["U3_1", "CNot_1_2", "Identity"], 
        "target_gate" => QCOpt.GroverDiffusionGate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "exact_optimal",

        "U3_θ_discretization" => -π:π/2:π,
        "U3_ϕ_discretization" => -π:π/2:π,
        "U3_λ_discretization" => [0])
end

function magic_using_SH_CNot_1_2()
    
    println(">>>>> M gate using S, H and CNOT_1_2 Gate <<<<<")

    return Dict{String, Any}(
        "num_qubits" => 2, 
        "maximum_depth" => 10,
        "elementary_gates" => ["S_1", "S_2", "H_1", "H_2", "CNot_1_2", "Identity"], 
        "target_gate" => QCOpt.MGate(),
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact_optimal") 
end


function magic_using_SH_CNot_2_1()
    
    println(">>>>> M gate using S, H and CNOT_2_1 Gate <<<<<")

    return Dict{String, Any}(
        "num_qubits" => 2, 
        "maximum_depth" => 10,
        "elementary_gates" => ["S_1", "S_2", "H_1", "H_2", "CNot_2_1", "Identity"], 
        "target_gate" => QCOpt.MGate(),
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact_optimal")
end

function magic_using_U3_CNot_1_2()
    
    println(">>>>> Magic basis using U3 and CNot_1_2 <<<<<")

    return Dict{String, Any}(
        "num_qubits" => 2, 
        "maximum_depth" => 5,
        "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "Identity"], 
        "target_gate" => QCOpt.MGate(),
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact_optimal",
        
        "U3_θ_discretization" => -π:π/2:π,
        "U3_ϕ_discretization" => -π:π:π,
        "U3_λ_discretization" => -π:π/2:π,)
end

function toffoli_with_controlled_gates()

    # Reference: https://doi.org/10.1109/TCAD.2005.858352
    println(">>>>> Toffoli gate using controlled gates <<<<<")
 
    return Dict{String, Any}(
    "num_qubits" => 3,
    "maximum_depth" => 5,
    "elementary_gates" => ["CV_1_3", "CV_2_3", "CV_1_2", "CVdagger_1_3", "CVdagger_2_3", "CVdagger_1_2", "CNot_2_1", "CNot_1_2", "Identity"],
    "target_gate" => QCOpt.ToffoliGate(),
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal")
end

function Fredkin()

    println(">>>>> Fredkin gate using controlled gates <<<<<")

    return Dict{String, Any}(
    "num_qubits" => 3,
    "maximum_depth" => 7,
    # Reference: https://doi.org/10.1103/PhysRevA.53.2855
    "elementary_gates" => ["CV_1_2", "CV_2_3", "CV_1_3", "CVdagger_1_2", "CVdagger_2_3", "CVdagger_1_3", "CNot_1_2", "CNot_3_2", "CNot_2_3", "CNot_1_3", "Identity"],
    "target_gate" => QCOpt.CSwapGate(), #also Fredkin
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal")
end

function double_toffoli()

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

    return Dict{String, Any}(
    "num_qubits" => 4,
    "maximum_depth" => 7,
    "elementary_gates" => ["CV_1_2", "CV_2_4", "CV_3_4", "CVdagger_1_2", "CVdagger_2_4", "CVdagger_3_4", "CNot_1_3", "CNot_3_2", "CNot_2_3", "Identity"],
    "target_gate" => target_gate(),
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal")
end

function quantum_fulladder()
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

    return Dict{String, Any}(
        "num_qubits" => num_qubits,
        "maximum_depth" => 7,
        "elementary_gates" => ["CV_1_2", "CV_4_2", "CV_3_2", "CVdagger_1_2", "CVdagger_4_2", "CVdagger_3_2", "CNot_3_1", "CNot_4_3", "CNot_2_4", "CNot_4_1", "Identity"],
        "target_gate" => target_gate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "exact_optimal")
end