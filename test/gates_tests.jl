# Unit tests for functions in gates.jl

@testset "building elementary universal gate tests" begin
    @test isapprox(QCO.U2Gate(0,-π/4), QCO.U3Gate(π/2,0,-π/4), atol=1E-6)
    @test isapprox(QCO.U1Gate(-π/4), QCO.U3Gate(0,0,-π/4), atol=1E-6)

    ctrl_qbit_1 = Array{Complex{Float64},2}([1 0; 0 0])
    ctrl_qbit_2 = Array{Complex{Float64},2}([0 0; 0 1])
    
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.XGate()), QCO.CXGate(), atol = 1E-6)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.YGate()), QCO.CYGate(), atol = 1E-6)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.ZGate()), QCO.CZGate(), atol = 1E-6)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.HGate()), QCO.CHGate(), atol = 1E-6)
    @test isapprox(QCO.CXGate(), QCO.CNotGate(), atol = 1E-6)

    ϴ1 = π/3
    ϴ2 = -2*π/3
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RXGate(ϴ1)), QCO.CRXGate(ϴ1), atol = 1E-6)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RXGate(ϴ2)), QCO.CRXGate(ϴ2), atol = 1E-6)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RYGate(ϴ1)), QCO.CRYGate(ϴ1), atol = 1E-6)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RYGate(ϴ2)), QCO.CRYGate(ϴ2), atol = 1E-6)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RZGate(ϴ1)), QCO.CRZGate(ϴ1), atol = 1E-6)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RZGate(ϴ2)), QCO.CRZGate(ϴ2), atol = 1E-6)

    θ = π/3
    ϕ = -2*π/3
    λ = π/6
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.U3Gate(θ,ϕ,λ)), QCO.CU3Gate(θ,ϕ,λ), atol = 1E-6)

    @test isapprox(QCO.SwapGate(), QCO.CNotGate() * QCO.CNotRevGate() * QCO.CNotGate(), atol = 1E-6)

    @test isapprox(QCO.iSwapGate(), QCO.CNotRevGate() * QCO.CNotGate() * kron(QCO.SGate(), QCO.IGate(1)) * QCO.CNotRevGate(), atol = 1E-6)

    @test isapprox(QCO.PhaseGate(λ), QCO.U3Gate(0,0,λ), atol = 1E-6)

    @test isapprox(QCO.DCXGate(), QCO.CNotGate() * QCO.CNotRevGate(), atol = 1E-6)

    S1_S2 = kron(QCO.SGate(), QCO.SGate())
    H2 = kron(QCO.IGate(1), QCO.HGate())
    @test isapprox(QCO.MGate(), QCO.CNotRevGate() * H2 * S1_S2)

    @test isapprox(QCO.NegIGate(), kron(QCO.RZGate(2*π), QCO.IGate(1)), atol = 1E-6)

    Z1 = QCO.get_full_sized_gate("Z1", 2);
    Z2 = QCO.get_full_sized_gate("Z2", 2);
    T2 = QCO.get_full_sized_gate("T2", 2);
    Y2 = QCO.get_full_sized_gate("Y2", 2);
    cnot_12 = QCO.get_full_sized_gate("cnot_12", 2);
    cnot_21 = QCO.get_full_sized_gate("cnot_21", 2);
    Sdagger1 = QCO.get_full_sized_gate("Sdagger1", 2);
    Tdagger1 = QCO.get_full_sized_gate("Tdagger1", 2);
    SX1 = QCO.get_full_sized_gate("SX1", 2);
    SXdagger2 = QCO.get_full_sized_gate("SXdagger2", 2);
    @test isapprox(-QCO.HCoinGate(), Z2 * cnot_21 * SXdagger2 * Y2 * cnot_21 * Tdagger1 * Z1 * T2 * Sdagger1 * cnot_12 * Sdagger1 * SX1 * cnot_21 * cnot_12, atol = 1E-6)

end

