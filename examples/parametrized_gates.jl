function parametrized_hermitian_gates()

    function target_gate(num_qubits, i, j, t)
        A = Array{Complex{Float64},2}(zeros(2^num_qubits, 2^num_qubits))

        A[i,j] = 1.0 + 0.0im
        A[j,i] = 1.0 + 0.0im 

        return exp(im*A*t)
    end

    num_qubits = 3
    A_i = 1 # 1 <= i,j <= num_qubits
    A_j = 3 # 1 <= i,j <= num_qubits
    t = 2.5

    return Dict{String, Any}(
    "num_qubits" => num_qubits, 
    "maximum_depth" => 8,
    # "elementary_gates" => ["RX_1", "RX_2", "CRX_1_2", "CRX_2_1", "CNot_1_2", "CNot_2_1", "Identity"],
    "elementary_gates" => ["RX_1", "RX_2", "RX_3", "CRX_1_2", "CRX_2_1", "CRX_2_3", "CRX_3_2", "CNot_1_2", "CNot_2_1", "CNot_2_3", "CNot_3_2", "Identity"],
    "target_gate" => target_gate(num_qubits, A_i, A_j, t),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal",

    "RX_discretization" => [t*2,-t*2],
    "CRX_discretization" => [t*2,-t*2],
    )
end