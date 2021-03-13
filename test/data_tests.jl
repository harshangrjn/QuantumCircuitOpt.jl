# Unit tests for functions in data.jl

@testset "building elementary universal gate tests" begin
    test_U3_1 = QCO.get_universal_gate(π/2,π/2,π/2)
    U3_1 = [(1/sqrt(2))+0.0im     0.0-(1/sqrt(2))im
             0.0+(1/sqrt(2))im   -(1/sqrt(2))+0.0im]
    @test ceil.(real(U3_1), digits=6) == ceil.(real(test_U3_1), digits=6)
    @test ceil.(imag(U3_1), digits=6) == ceil.(imag(test_U3_1), digits=6)

    test_U3_2 = QCO.get_universal_gate(π,π/2,π/2)
    G = QCO.get_elementary_gates(2)
    @test test_U3_2 == G["pauli_Y"]
end 