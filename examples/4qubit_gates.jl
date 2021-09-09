function decompose_CNot_41()
    # Reference: https://doi.org/10.1109/DSD.2018.00005

    params = Dict{String, Any}(
    "num_qubits" => 4,
    "depth" => 10,

    "elementary_gates" => ["H_1", "H_2", "H_3", "CNot_13", "CNot_43", "Identity"],
    "target_gate" => QCO.get_full_sized_gate("CNot_41", 4),
    "identify_real_gates" => true,

    "objective" => "minimize_depth", 
    "optimizer" => "cplex"
    
    )

    return params
end
