# Unit tests for functions in data.jl

@testset "building elementary universal gate tests" begin
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

@testset "get_full_sized_gate tests" begin

    # 2-qubit gates
    params = Dict{String, Any}(
    "num_qubits" => 2,
    "depth" => 2,
    "elementary_gates" => ["T_1", "T_2", "Tdagger_1", "Tdagger_2", "S_1", "S_2", "Sdagger_1", "Sdagger_2", "SX_1", "SX_2", "SXdagger_1", "SXdagger_2", "X_1", "X_2", "Y_1", "Y_2", "Z_1", "Z_2", "CNotSwap", "H_1⊗H_2", "CZ_12", "CH_12", "CV_12", "Swap", "M_12", "QFT_12", "CSX_12", "W_12", "HCoin"],
    "target_gate" => QCO.IGate(2),               
    )

    data = QCO.get_data(params, eliminate_identical_gates = false)
    @test length(keys(data["gates_dict"])) == 29 

    # 3-qubit gates
    params = Dict{String, Any}(
    "num_qubits" => 3,
    "depth" => 2,
    "elementary_gates" => ["H_3", "T_3", "Tdagger_3", "Sdagger_3", "SX_3", "SXdagger_3", "X_3", "Y_3", "Z_3", "Toffoli", "CSwap", "CCZ", "Peres", "CNot_12", "CNot_23", "CNot_21", "CNot_32", "CNot_13", "CNot_31"],
    "target_gate" => QCO.IGate(3)
    )

    data = QCO.get_data(params)
    @test length(keys(data["gates_dict"])) == 19
end

@testset "get_input_circuit_dict tests" begin
    
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
