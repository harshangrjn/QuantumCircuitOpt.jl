# Unit tests for functions in data.jl

@testset "building elementary universal gate tests" begin
    test_angle = π/3
    pauli_Y = QCO.get_elementary_gates(2)["pauli_Y"]
    H = QCO.get_elementary_gates(2)["hadamard_H"]
    R_gates = QCO.get_pauli_rotationum_gates(test_angle)
    U3_1 = [(1/sqrt(2))+0.0im     0.0-(1/sqrt(2))im
             0.0+(1/sqrt(2))im   -(1/sqrt(2))+0.0im]

    test_U3_1 = QCO.get_u3_gate(π/2,π/2,π/2)
    @test isapprox(U3_1, test_U3_1)
    # @test ceil.(real(U3_1), digits=6) == ceil.(real(test_U3_1), digits=6)
    # @test ceil.(imag(U3_1), digits=6) == ceil.(imag(test_U3_1), digits=6)

    test_U3_2 = QCO.get_u3_gate(π,π/2,π/2)
    @test isapprox(test_U3_2, pauli_Y)

    test_U2_2 = QCO.get_u2_gate(0,π)
    @test isapprox(test_U2_2, H)

    test_U3_3 = QCO.get_u3_gate(test_angle, -π/2, π/2)
    @test isapprox(R_gates["R_x"], test_U3_3)

    test_U3_4 = QCO.get_u3_gate(test_angle, 0, 0)
    @test isapprox(R_gates["R_y"], test_U3_4)

    test_U1_1 = exp(-((test_angle)/2)im)*QCO.get_u1_gate(test_angle)
    @test isapprox(R_gates["R_z"], test_U1_1)
end 