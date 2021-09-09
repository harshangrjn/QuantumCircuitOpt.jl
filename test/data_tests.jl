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
    "elementary_gates" => ["T_1", "T_2", "Tdagger_1", "Tdagger_2", "S_1", "S_2", "Sdagger_1", "Sdagger_2", "SX_1", "SX_2", "SXdagger_1", "SXdagger_2", "X_1", "X_2", "Y_1", "Y_2", "Z_1", "Z_2", "CZ_12", "CH_12", "CV_12", "Phase_1", "Phase_2", "DCX_12"],
    "Ph_discretization" => [π],
    "target_gate" => QCO.IGate(2),               
    )

    data = QCO.get_data(params, eliminate_identical_gates = false)
    @test length(keys(data["gates_dict"])) == 24
    
    # kron gates
    params = Dict{String, Any}(
    "num_qubits" => 3,
    "depth" => 2,
    "elementary_gates" => ["I_1xH_2xT_3", "Tdagger_1xS_2xSdagger_3", "SX_1xSXdagger_2xX_3", "Y_1xCNot_23", "CNot_21xZ_3", "CV_21xI_3", "I_1xCV_23", "CVdagger_21xI_3", "I_1xCVdagger_23", "CX_21xI_3", "I_1xCX_23", "CY_21xI_3", "I_1xCY_23", "CZ_21xI_3", "I_1xCZ_23", "CH_21xI_3", "I_1xCH_23", "CSX_21xI_3", "I_1xCSX_23", "Swap_21xI_3", "I_1xiSwap_23", "DCX_21xI_3"],
    "target_gate" => QCO.IGate(3),         
    )

    data = QCO.get_data(params, eliminate_identical_gates = false)
    @test length(keys(data["gates_dict"])) == 22

    # 3-qubit gates
    params = Dict{String, Any}(
    "num_qubits" => 3,
    "depth" => 2,
    "elementary_gates" => ["H_3", "T_3", "Tdagger_3", "Sdagger_3", "SX_3", "SXdagger_3", "X_3", "Y_3", "Z_3", "CNot_12", "CNot_23", "CNot_21", "CNot_32", "CNot_13", "CNot_31"],
    "target_gate" => QCO.IGate(3)
    )

    # 4-qubit gates
    params = Dict{String, Any}(
        "num_qubits" => 4,
        "depth" => 2,
        "elementary_gates" => ["CU3_13", "CRX_34", "CNot_24", "CH_14"],
        "CRX_discretization" => [π],
        "CU_θ_discretization" => [π/2],
        "CU_ϕ_discretization" => [-π/2],
        "CU_λ_discretization" => [π/4],
        "target_gate" => QCO.IGate(4)
        )
    data = QCO.get_data(params)
    @test length(keys(data["gates_dict"])) == 4

    #5-qubit gates
    params = Dict{String, Any}(
        "num_qubits" => 5,
        "depth" => 2,
        "elementary_gates" => ["H_1", "T_2", "Tdagger_3", "Sdagger_4", "SX_5", "CX_12", "CX_21", "CY_12", "CY_21", "CRX_23", "CNot_34", "CH_54", "CRY_13", "CRZ_24", "CV_35", "CZ_41", "CSX_51"],
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
        return [(1, "CNot_21"), 
                (2, "S_1"), 
                (3, "H_2"), 
                (4, "S_2")
                ]
    end

    params = Dict{String, Any}(
    "num_qubits" => 2,
    "depth" => 5,
    "elementary_gates" => ["S_1", "S_2", "H_1", "H_2", "CNot_12", "CNot_21", "Identity"], 
    "target_gate" => QCO.MGate(),
    "input_circuit" => input_circuit_1(),
    )

    data = QCO.get_data(params)

    @test "input_circuit" in keys(data)
    @test data["input_circuit"]["1"]["depth"] == 1
    @test data["input_circuit"]["1"]["gate"] == "CNot_21"
    @test data["input_circuit"]["2"]["depth"] == 2
    @test data["input_circuit"]["2"]["gate"] == "S_1"
    @test data["input_circuit"]["3"]["depth"] == 3
    @test data["input_circuit"]["3"]["gate"] == "H_2"
    @test data["input_circuit"]["4"]["depth"] == 4
    @test data["input_circuit"]["4"]["gate"] == "S_2"

    function input_circuit_2()
        # [(depth, gate)]
        return [(1, "CNot_21"), 
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
        return [(1, "CNot_21"), 
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
        return [(1, "CNot_21"), 
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
