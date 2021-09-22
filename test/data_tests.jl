# Unit tests for functions in data.jl

@testset "Tests: building elementary universal gate" begin
    test_angle = π/3
    pauli_Y = QCO.YGate()
    H = QCO.HGate()
    U3_1 = [(1/sqrt(2))+0.0im     0.0-(1/sqrt(2))im
             0.0+(1/sqrt(2))im   -(1/sqrt(2))+0.0im]

    test_U3_1 = QCO.U3Gate(π/2,π/2,π/2)
    @test isapprox(U3_1, test_U3_1)
    # @test ceil.(real(U3_1), digits=6) == ceil.(real(test_U3_1), digits=6)
    # @test ceil.(imag(U3_1), digits=6) == ceil.(imag(test_U3_1), digits=6)

    test_U3_2 = QCO.U3Gate(π,π/2,π/2)
    @test isapprox(test_U3_2, pauli_Y)

    test_U2_2 = QCO.U2Gate(0,π)
    @test isapprox(test_U2_2, H)

    test_U3_3 = QCO.U3Gate(test_angle, -π/2, π/2)
    @test isapprox(QCO.RXGate(test_angle), test_U3_3)

    test_U3_4 = QCO.U3Gate(test_angle, 0, 0)
    @test isapprox(QCO.RYGate(test_angle), test_U3_4)

    test_U1_1 = exp(-((test_angle)/2)im)*QCO.U1Gate(test_angle)
    @test isapprox(QCO.RZGate(test_angle), test_U1_1)
end 

@testset "Tests: get_full_sized_gate" begin

    # 2-qubit gates
    params = Dict{String, Any}(
    "num_qubits" => 2,
    "depth" => 2,
    "elementary_gates" => ["T_1", "T_2", "Tdagger_1", "Tdagger_2", "S_1", "S_2", "Sdagger_1", "Sdagger_2", "SX_1", "SX_2", "SXdagger_1", "SXdagger_2", "X_1", "X_2", "Y_1", "Y_2", "Z_1", "Z_2", "CZ_1_2", "CH_1_2", "CV_1_2", "Phase_1", "Phase_2", "CSX_1_2", "DCX_1_2", "Sycamore_1_2"],
    "Phase_discretization" => [π],
    "target_gate" => QCO.IGate(2),               
    )

    data = QCO.get_data(params, eliminate_identical_gates = false)
    @test length(keys(data["gates_dict"])) == 26
    
    # kron gates
    params = Dict{String, Any}(
    "num_qubits" => 3,
    "depth" => 2,
    "elementary_gates" => ["I_1xH_2xT_3", "Tdagger_1xS_2xSdagger_3", "SX_1xSXdagger_2xX_3", "Y_1xCNot_2_3", "CNot_2_1xZ_3", "CV_2_1xI_3", "I_1xCV_2_3", "CVdagger_2_1xI_3", "I_1xCVdagger_2_3", "CX_2_1xI_3", "I_1xCX_2_3", "CY_2_1xI_3", "I_1xCY_2_3", "CZ_2_1xI_3", "I_1xCZ_2_3", "CH_2_1xI_3", "I_1xCH_2_3", "CSX_2_1xI_3", "I_1xCSX_2_3", "Swap_2_1xI_3", "I_1xiSwap_2_3", "DCX_2_1xI_3", "I_1xSycamore_2_3"],
    "target_gate" => QCO.IGate(3),         
    )

    data = QCO.get_data(params, eliminate_identical_gates = false)
    @test length(keys(data["gates_dict"])) == 23

    # 3-qubit gates
    params = Dict{String, Any}(
    "num_qubits" => 3,
    "depth" => 2,
    "elementary_gates" => ["H_3", "T_3", "Tdagger_3", "Sdagger_3", "SX_3", "SXdagger_3", "X_3", "Y_3", "Z_3", "CNot_1_2", "CNot_2_3", "CNot_2_1", "CNot_3_2", "CNot_1_3", "CNot_3_1"],
    "target_gate" => QCO.IGate(3)
    )

    # 4-qubit gates
    params = Dict{String, Any}(
        "num_qubits" => 4,
        "depth" => 2,
        "elementary_gates" => ["CU3_1_3", "CRX_3_4", "CNot_2_4", "CH_1_4"],
        "CRX_discretization" => [π],
        "CU3_θ_discretization" => [π/2],
        "CU3_ϕ_discretization" => [-π/2],
        "CU3_λ_discretization" => [π/4],
        "target_gate" => QCO.IGate(4)
        )
    data = QCO.get_data(params)
    @test length(keys(data["gates_dict"])) == 4

    #5-qubit gates
    params = Dict{String, Any}(
        "num_qubits" => 5,
        "depth" => 2,
        "elementary_gates" => ["H_1", "T_2", "Tdagger_3", "Sdagger_4", "SX_5", "CX_1_2", "CX_2_1", "CY_1_2", "CY_2_1", "CRX_2_3", "CNot_3_4", "CH_5_4", "CRY_1_3", "CRZ_2_4", "CV_3_5", "CZ_4_1", "CSX_5_1"],
        "CRX_discretization" => [π],
        "CRY_discretization" => [π],
        "CRZ_discretization" => [π],
        "target_gate" => QCO.IGate(5)
        )

    data = QCO.get_data(params)
    @test length(keys(data["gates_dict"])) == 17
end

@testset "Tests: get_input_circuit_dict" begin
    
    function input_circuit_1()
        # [(depth, gate)]
        return [(1, "CNot_2_1"), 
                (2, "S_1"), 
                (3, "H_2"), 
                (4, "S_2")
                ]
    end

    params = Dict{String, Any}(
    "num_qubits" => 2,
    "depth" => 5,
    "elementary_gates" => ["S_1", "S_2", "H_1", "H_2", "CNot_1_2", "CNot_2_1", "Identity"], 
    "target_gate" => QCO.MGate(),
    "input_circuit" => input_circuit_1(),
    )

    data = QCO.get_data(params)

    @test "input_circuit" in keys(data)
    @test data["input_circuit"]["1"]["depth"] == 1
    @test data["input_circuit"]["1"]["gate"] == "CNot_2_1"
    @test data["input_circuit"]["2"]["depth"] == 2
    @test data["input_circuit"]["2"]["gate"] == "S_1"
    @test data["input_circuit"]["3"]["depth"] == 3
    @test data["input_circuit"]["3"]["gate"] == "H_2"
    @test data["input_circuit"]["4"]["depth"] == 4
    @test data["input_circuit"]["4"]["gate"] == "S_2"

    function input_circuit_2()
        # [(depth, gate)]
        return [(1, "CNot_2_1"), 
                (2, "S_1"), 
                (3, "H_2"), 
                (4, "T_1")
                ]
    end

    params["input_circuit"] = input_circuit_2()
    data = QCO.get_data(params)
    @test !("input_circuit" in keys(data))

    function input_circuit_3()
        # [(depth, gate)]
        return [(1, "CNot_2_1"), 
                (2, "S_1"), 
                (2, "H_2"), 
                (4, "S_2")
                ]
    end

    params["input_circuit"] = input_circuit_3()
    data = QCO.get_data(params)
    @test !("input_circuit" in keys(data))

    function input_circuit_4()
        # [(depth, gate)]
        return [(1, "CNot_2_1"), 
                (2, "S_1"), 
                (3, "H_2"), 
                (4, "S_2"),
                (5, "Identity"),
                (6, "S_1"),
                ]
    end

    params["input_circuit"] = input_circuit_4()
    data = QCO.get_data(params)
    @test !("input_circuit" in keys(data))

end
