@testset "types testing" begin
    gate = QCO.GateData("H1", 2)
    @test gate.isreal
    gate = QCO.GateData("S3", 3)
    @test !gate.isreal

    # Add more tests for U3 and R gates.
end