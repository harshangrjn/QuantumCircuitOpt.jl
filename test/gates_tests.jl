# Unit tests for functions in gates.jl
@testset "Gates tests: elementary universal gates in gates.jl" begin
    @test isapprox(QCO.U2Gate(0, -π/4), QCO.U3Gate(π/2, 0, -π/4), atol=tol_0)
    @test isapprox(QCO.U1Gate(-π/4),   QCO.U3Gate(0,   0, -π/4),   atol=tol_0)

    ctrl_qbit_1 = Matrix{ComplexF64}([1 0; 0 0])
    ctrl_qbit_2 = Matrix{ComplexF64}([0 0; 0 1])
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.XGate()), QCO.CXGate(), atol=tol_0)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.YGate()), QCO.CYGate(), atol=tol_0)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.ZGate()), QCO.CZGate(), atol=tol_0)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.HGate()), QCO.CHGate(), atol=tol_0)
    @test isapprox(QCO.CXGate(), QCO.CNotGate(), atol=tol_0)

    ϴ1 = π/3
    ϴ2 = -2π/3
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RXGate(ϴ1)), QCO.CRXGate(ϴ1), atol=tol_0)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RXGate(ϴ2)), QCO.CRXGate(ϴ2), atol=tol_0)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RYGate(ϴ1)), QCO.CRYGate(ϴ1), atol=tol_0)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RYGate(ϴ2)), QCO.CRYGate(ϴ2), atol=tol_0)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RZGate(ϴ1)), QCO.CRZGate(ϴ1), atol=tol_0)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RZGate(ϴ2)), QCO.CRZGate(ϴ2), atol=tol_0)

    θ = π/3
    ϕ = -2π/3
    λ = π/6
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.U3Gate(θ, ϕ, λ)), QCO.CU3Gate(θ, ϕ, λ), atol=tol_0)
    @test isapprox(QCO.SwapGate(), QCO.CNotGate() * QCO.CNotRevGate() * QCO.CNotGate(), atol=tol_0)
    @test isapprox(QCO.iSwapGate(), QCO.CNotRevGate() * QCO.CNotGate() * kron(QCO.SGate(), QCO.IGate(1)) * QCO.CNotRevGate(), atol=tol_0)
    @test isapprox(QCO.PhaseGate(λ), QCO.U3Gate(0, 0, λ), atol=tol_0)
    @test isapprox(QCO.DCXGate(), QCO.CNotGate() * QCO.CNotRevGate(), atol=tol_0)

    S1_S2 = kron(QCO.SGate(), QCO.SGate())
    H2    = kron(QCO.IGate(1), QCO.HGate())
    @test isapprox(QCO.MGate(), QCO.CNotRevGate() * H2 * S1_S2)

    Z1       = QCO.unitary("Z_1", 2)
    Z2       = QCO.unitary("Z_2", 2)
    T2       = QCO.unitary("T_2", 2)
    Y_2      = QCO.unitary("Y_2", 2)
    CNot_1_2 = QCO.unitary("CNot_1_2", 2)
    CNot_2_1 = QCO.unitary("CNot_2_1", 2)
    Sdagger1 = QCO.unitary("Sdagger_1", 2)
    Tdagger1 = QCO.unitary("Tdagger_1", 2)
    SX1      = QCO.unitary("SX_1", 2)
    SXdagger2 = QCO.unitary("SXdagger_2", 2)
    @test isapprox(-QCO.HCoinGate(), Z2 * CNot_2_1 * SXdagger2 * Y_2 * CNot_2_1 * Tdagger1 * Z1 *
                                     T2 * Sdagger1 * CNot_1_2 * Sdagger1 * SX1 * CNot_2_1 * CNot_1_2,
                   atol=tol_0)

    CV_23       = QCO.unitary("CV_2_3", 3)
    CNot_1_2    = QCO.unitary("CNot_1_2", 3)
    CVdagger_23 = QCO.unitary("CVdagger_2_3", 3)
    CV_13       = QCO.unitary("CV_1_3", 3)
    @test isapprox(QCO.ToffoliGate(), CV_23 * CNot_1_2 * CVdagger_23 * CNot_1_2 * CV_13, atol=tol_0)

    CVdagger_12 = QCO.unitary("CVdagger_1_2", 3)
    CVdagger_31 = QCO.unitary("CVdagger_3_1", 3)
    @test isapprox(QCO.IGate(3), CVdagger_12^4 * CVdagger_31^4, atol=tol_0)

    H_2      = QCO.unitary("H_2", 2)
    T_1      = QCO.unitary("T_1", 2)
    Tdagger2 = QCO.unitary("Tdagger_2", 2)
    @test isapprox(QCO.CVdaggerGate(), H_2 * Tdagger1 * CNot_2_1 * T_1 * Tdagger2 * CNot_2_1 * H_2, atol=tol_0)
    @test isapprox(QCO.CVGate(), QCO.CNotGate() * QCO.CVdaggerGate(), atol=tol_0)

    # Next 3 tests from: https://american-cse.org/csci2015/data/9795a059.pdf
    CV_31       = QCO.unitary("CV_3_1", 3)
    CNot_2_3    = QCO.unitary("CNot_2_3", 3)
    CVdagger_31 = QCO.unitary("CVdagger_3_1", 3)
    CV_21       = QCO.unitary("CV_2_1", 3)
    CVdagger_21 = QCO.unitary("CVdagger_2_1", 3)
    @test isapprox(CV_31 * CNot_2_3 * CVdagger_31 * CNot_2_3 * CV_21,
                   CVdagger_21 * CNot_2_3 * CV_31 * CNot_2_3 * CVdagger_31, atol=tol_0)

    CV_12       = QCO.unitary("CV_1_2", 3)
    CV_32       = QCO.unitary("CV_3_2", 3)
    CVdagger_32 = QCO.unitary("CVdagger_3_2", 3)
    CVdagger_12 = QCO.unitary("CVdagger_1_2", 3)
    CNot_1_3    = QCO.unitary("CNot_1_3", 3)
    @test isapprox(CV_32 * CNot_1_3 * CVdagger_32 * CNot_1_3 * CV_12,
                   CVdagger_32 * CNot_1_3 * CV_32 * CNot_1_3 * CVdagger_12, atol=tol_0)

    CVdagger_13 = QCO.unitary("CVdagger_1_3", 3)
    CNot_2_1     = QCO.unitary("CNot_2_1", 3);
    @test isapprox(CV_13 * CNot_2_1 * CVdagger_13 * CNot_2_1 * CV_23,
                   CVdagger_13 * CNot_2_1 * CV_13 * CNot_2_1 * CVdagger_23, atol=tol_0)

    # 2-qubit Grover's diffusion operator
    H_1      = QCO.unitary("H_1", 2)
    X_2      = QCO.unitary("X_2", 2)
    CNot_1_2 = QCO.unitary("CNot_1_2", 2)
    @test isapprox(QCO.GroverDiffusionGate(), Y_2 * H_1 * X_2 * CNot_1_2 * H_1 * Y_2, atol=tol_0)

    # 2-qubit QFT
    SWAP = QCO.unitary("Swap_1_2", 2)
    CU   = QCO.unitary("CU3_2_1", 2, angle=[0, π/4, π/4])
    @test isapprox(QCO.QFT2Gate(), H_1 * CU * H_2 * SWAP)

    # 3-qubit QFT
    CS_2_1   = QCO.unitary("CS_2_1", 3)
    CS_3_2   = QCO.unitary("CS_3_2", 3)
    CT_3_1   = QCO.unitary("CT_3_1", 3)
    SWAP_1_3 = QCO.unitary("Swap_1_3", 3)
    H_1      = QCO.unitary("H_1", 3)
    H_2      = QCO.unitary("H_2", 3)
    H_3      = QCO.unitary("H_3", 3)
    @test isapprox(QCO.QFT3Gate(), H_1 * CS_2_1 * CT_3_1 * H_2 * CS_3_2 * H_3 * SWAP_1_3)

    # CRY Decomp
    CNOT12   = QCO.unitary("CNot_1_2", 2)
    U3_negpi = QCO.unitary("U3_2", 2, angle=[-π/2, 0, 0])
    U3_pi    = QCO.unitary("U3_2", 2, angle=[ π/2,  0, 0])
    @test isapprox(QCO.CRYGate(π), CNOT12 * U3_negpi * CNOT12 * U3_pi, atol=tol_0)

    # Identity tests
    Id     = QCO.IGate(2)
    CRXRev = QCO.unitary("CRX_2_1", 2, angle=π)
    @test isapprox(Id, CRXRev^4, atol=tol_0)
    CRYRev = QCO.unitary("CRY_2_1", 2, angle=π)
    @test isapprox(Id, CRYRev^4, atol=tol_0)
    CRZRev = QCO.unitary("CRZ_2_1", 2, angle=π)
    @test isapprox(Id, CRZRev^4, atol=tol_0)

    # Rev gate tests for involution
    @test isapprox(QCO.CXRevGate(), QCO.CNotRevGate(), atol=tol_0)
    @test isapprox(QCO.CXGate() * QCO.SwapGate() * QCO.CXRevGate() * QCO.SwapGate(), QCO.IGate(2), atol=tol_0)
    @test isapprox(QCO.CYGate() * QCO.SwapGate() * QCO.CYRevGate() * QCO.SwapGate(), QCO.IGate(2), atol=tol_0)
    @test isapprox(QCO.CZGate() * QCO.SwapGate() * QCO.CZGate()    * QCO.SwapGate(), QCO.IGate(2), atol=tol_0)
    @test isapprox(QCO.CHGate() * QCO.SwapGate() * QCO.CHRevGate() * QCO.SwapGate(), QCO.IGate(2), atol=tol_0)

    # Rev gate tests for non-involution
    @test isapprox(QCO.CVGate() * QCO.SwapGate() * QCO.CVRevGate() * QCO.SwapGate(), QCO.CVGate()^2, atol=tol_0)
    @test isapprox(QCO.CSXGate() * QCO.SwapGate() * QCO.CSXRevGate() * QCO.SwapGate(), QCO.CSXGate()^2, atol=tol_0)
    @test isapprox(QCO.WGate(), QCO.HCoinGate(), atol=tol_0)

    # Fredkin = CSwapGate
    CNot_3_2    = QCO.unitary("CNot_3_2", 3)
    CV_23       = QCO.unitary("CV_2_3", 3)
    CV_13       = QCO.unitary("CV_1_3", 3)
    CNot_1_2    = QCO.unitary("CNot_1_2", 3)
    CVdagger_23 = QCO.unitary("CVdagger_2_3", 3)
    @test isapprox(QCO.CSwapGate(), CNot_3_2 * CV_23 * CV_13 * CNot_1_2 * CVdagger_23 * CNot_1_2 * CNot_3_2, atol=tol_0)

    # CCZGate test (Toffoli with target qubit conjugated by Hadamard)
    H_3 = QCO.unitary("H_3", 3)
    @test isapprox(QCO.ToffoliGate(), H_3 * QCO.CCZGate() * H_3, atol=tol_0)

    # Peres test
    CVdagger_13 = QCO.unitary("CVdagger_1_3", 3)
    @test isapprox(QCO.PeresGate(), QCO.ToffoliGate() * CNot_1_2, atol=tol_0)
    @test isapprox(QCO.PeresGate(), CVdagger_13 * CVdagger_23 * CNot_1_2 * CV_23, atol=tol_0)

    # iSwap test
    @test isapprox(QCO.unitary("iSwap_2_1", 3), QCO.unitary("iSwap_1_2", 3), atol=tol_0)
    @test isapprox(QCO.unitary("iSwap_2_3", 4), QCO.unitary("iSwap_3_2", 4), atol=tol_0)
    @test isapprox(QCO.unitary("iSwap_3_5", 5), QCO.unitary("iSwap_5_3", 5), atol=tol_0)

    # Sycamore gate test
    @test isapprox(QCO.SycamoreGate()^12, QCO.IGate(2), atol=tol_0)

    # RCCX test
    T_3       = QCO.unitary("T_3", 3)
    Tdagger_3 = QCO.unitary("Tdagger_3", 3)
    @test isapprox(QCO.RCCXGate(), H_3 * Tdagger_3 * CNot_2_3 * T_3 * CNot_1_3 * Tdagger_3 * CNot_2_3 * T_3 * H_3, atol=tol_0)

    # Margolus test
    RY_a = QCO.unitary("RY_3", 3, angle=-π/4)
    RY_b = QCO.unitary("RY_3", 3, angle= π/4)
    @test isapprox(QCO.MargolusGate(), RY_a * CNot_2_3 * RY_a * CNot_1_3 * RY_b * CNot_2_3 * RY_b, atol=tol_0)

    # CiSwap test
    @test isapprox(QCO.CiSwapGate(), CNot_3_2 * H_3 * CNot_2_3 * T_3 * CNot_1_3 * Tdagger_3 * CNot_2_3 * T_3 *
                                     CNot_1_3 * Tdagger_3 * H_3 * CNot_3_2, atol=tol_0)

    # 2-angle RGate test
    θ = π/3
    ϕ = π/6
    R1 = QCO.RGate(θ, ϕ)
    U3 = QCO.U3Gate(θ, ϕ, -ϕ)
    U3[1,2] *= im
    U3[2,1] *= -im
    @test isapprox(R1, U3, atol=tol_0)

    # CS, CSdagger gate test
    H1    = QCO.unitary("H_1", 2)
    H2    = QCO.unitary("H_2", 2)
    CS_12 = QCO.unitary("CS_1_2", 2)
    @test isapprox(H1 * CS_12 * H2 * QCO.SwapGate(), QCO.QFT2Gate(), atol=tol_0)
    @test isapprox(QCO.CSGate() * QCO.CSdaggerGate(), QCO.IGate(2), atol=tol_0)

    # SSwapGate test
    @test isapprox(QCO.SSwapGate()^2, QCO.SwapGate(), atol=1e-6)

    # CT, CTdagger gate test
    @test isapprox(QCO.CTGate() * QCO.CTdaggerGate(), QCO.IGate(2), atol=tol_0)

    # Tests for multi-qubit symbol expansions
    M1 = QCO.unitary("I_1xCZ_2_4xH_5", 5)
    M2 = kron(kron(QCO.IGate(1), QCO.unitary("CZ_1_3", 3)), QCO.HGate())
    @test isapprox(M1, M2, atol=tol_0)

    M1 = QCO.unitary("CV_4_1xH_5xZ_6", 6)
    M2 = kron(kron(QCO.unitary("CV_4_1", 4), QCO.HGate()), QCO.ZGate())
    @test isapprox(M1, M2, atol=tol_0)

    # Test multi_controlled_gate with an ancilla qubit
    M1 = QCO.multi_controlled_gate(QCO.HGate(), [1], 2, 3)
    M2 = kron(QCO.CHGate(), QCO.IGate(1))
    @test isapprox(M1, M2, atol=tol_0)

    # Test Fermionic-Swap gate
    @test isapprox(QCO.fSwapGate(), QCO.SwapGate()*QCO.CZGate(), atol=tol_0)
end
