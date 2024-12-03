@testset "Tests: types" begin
    gate = QCO.GateData("H_1", 2)
    @test gate.isreal
    gate = QCO.GateData("S_3", 3)
    @test !gate.isreal

    # Add more tests for U3 and R gates.
end
