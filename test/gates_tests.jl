# Unit tests for functions in gates.jl

@testset "Tests: elementary universal gates in gates.jl" begin
    @test isapprox(QCO.U2Gate(0,-π/4), QCO.U3Gate(π/2,0,-π/4), atol=tol_0)
    @test isapprox(QCO.U1Gate(-π/4), QCO.U3Gate(0,0,-π/4), atol=tol_0)

    ctrl_qbit_1 = Array{Complex{Float64},2}([1 0; 0 0])
    ctrl_qbit_2 = Array{Complex{Float64},2}([0 0; 0 1])
    
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.XGate()), QCO.CXGate(), atol = tol_0)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.YGate()), QCO.CYGate(), atol = tol_0)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.ZGate()), QCO.CZGate(), atol = tol_0)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.HGate()), QCO.CHGate(), atol = tol_0)
    @test isapprox(QCO.CXGate(), QCO.CNotGate(), atol = tol_0)

    ϴ1 = π/3
    ϴ2 = -2*π/3
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RXGate(ϴ1)), QCO.CRXGate(ϴ1), atol = tol_0)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RXGate(ϴ2)), QCO.CRXGate(ϴ2), atol = tol_0)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RYGate(ϴ1)), QCO.CRYGate(ϴ1), atol = tol_0)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RYGate(ϴ2)), QCO.CRYGate(ϴ2), atol = tol_0)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RZGate(ϴ1)), QCO.CRZGate(ϴ1), atol = tol_0)
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.RZGate(ϴ2)), QCO.CRZGate(ϴ2), atol = tol_0)

    θ = π/3
    ϕ = -2*π/3
    λ = π/6
    @test isapprox(kron(ctrl_qbit_1, QCO.IGate(1)) + kron(ctrl_qbit_2, QCO.U3Gate(θ,ϕ,λ)), QCO.CU3Gate(θ,ϕ,λ), atol = tol_0)

    @test isapprox(QCO.SwapGate(), QCO.CNotGate() * QCO.CNotRevGate() * QCO.CNotGate(), atol = tol_0)

    @test isapprox(QCO.iSwapGate(), QCO.CNotRevGate() * QCO.CNotGate() * kron(QCO.SGate(), QCO.IGate(1)) * QCO.CNotRevGate(), atol = tol_0)

    @test isapprox(QCO.PhaseGate(λ), QCO.U3Gate(0,0,λ), atol = tol_0)

    @test isapprox(QCO.DCXGate(), QCO.CNotGate() * QCO.CNotRevGate(), atol = tol_0)

    S1_S2 = kron(QCO.SGate(), QCO.SGate())
    H2 = kron(QCO.IGate(1), QCO.HGate())
    @test isapprox(QCO.MGate(), QCO.CNotRevGate() * H2 * S1_S2)

    Z1 = QCO.get_full_sized_gate("Z_1", 2);
    Z2 = QCO.get_full_sized_gate("Z_2", 2);
    T2 = QCO.get_full_sized_gate("T_2", 2);
    Y_2 = QCO.get_full_sized_gate("Y_2", 2);
    cnot_12 = QCO.get_full_sized_gate("CNot_1_2", 2);
    cnot_21 = QCO.get_full_sized_gate("CNot_2_1", 2);
    Sdagger1 = QCO.get_full_sized_gate("Sdagger_1", 2);
    Tdagger1 = QCO.get_full_sized_gate("Tdagger_1", 2);
    SX1 = QCO.get_full_sized_gate("SX_1", 2);
    SXdagger2 = QCO.get_full_sized_gate("SXdagger_2", 2);
    @test isapprox(-QCO.HCoinGate(), Z2 * cnot_21 * SXdagger2 * Y_2 * cnot_21 * Tdagger1 * Z1 * T2 * Sdagger1 * cnot_12 * Sdagger1 * SX1 * cnot_21 * cnot_12, atol = tol_0)

    CV_23 = QCO.get_full_sized_gate("CV_2_3", 3);
    cnot_12 = QCO.get_full_sized_gate("CNot_1_2", 3);
    CVdagger_23 = QCO.get_full_sized_gate("CVdagger_2_3", 3);
    CV_13 = QCO.get_full_sized_gate("CV_1_3", 3);
    @test isapprox(QCO.ToffoliGate(), CV_23 * cnot_12 * CVdagger_23 * cnot_12 * CV_13, atol = tol_0)

    CVdagger_12 = QCO.get_full_sized_gate("CVdagger_1_2", 3);
    CVdagger_31 = QCO.get_full_sized_gate("CVdagger_3_1", 3);
    @test isapprox(QCO.IGate(3), CVdagger_12 * CVdagger_12 * CVdagger_12 * CVdagger_12 * CVdagger_31 * CVdagger_31 * CVdagger_31 * CVdagger_31, atol = tol_0)

    H_2 = QCO.get_full_sized_gate("H_2", 2);
    T_1 = QCO.get_full_sized_gate("T_1", 2);
    Tdagger2 = QCO.get_full_sized_gate("Tdagger_2", 2);
    @test isapprox(QCO.CVdaggerGate(), H_2 * Tdagger1 * cnot_21 * T_1 * Tdagger2 * cnot_21 * H_2, atol = tol_0)

    cnot_12 = QCO.get_full_sized_gate("CNot_1_2", 2);
    @test isapprox(QCO.CVGate(), cnot_12 * QCO.CVdaggerGate(), atol = tol_0)

    # Next 3 tests from: https://american-cse.org/csci2015/data/9795a059.pdf

    CV_31 = QCO.get_full_sized_gate("CV_3_1", 3);
    CNot_23 = QCO.get_full_sized_gate("CNot_2_3", 3);
    CVdagger_31 = QCO.get_full_sized_gate("CVdagger_3_1", 3);
    CV_21 = QCO.get_full_sized_gate("CV_2_1", 3);

    CVdagger_21 = QCO.get_full_sized_gate("CVdagger_2_1", 3);

    @test isapprox(CV_31 * CNot_23 * CVdagger_31 * CNot_23 * CV_21, CVdagger_21 * CNot_23 * CV_31 * CNot_23 * CVdagger_31, atol = tol_0)

    CV_12 = QCO.get_full_sized_gate("CV_1_2", 3);
    CV_32 = QCO.get_full_sized_gate("CV_3_2", 3);
    CVdagger_32 = QCO.get_full_sized_gate("CVdagger_3_2", 3);
    CVdagger_12 = QCO.get_full_sized_gate("CVdagger_1_2", 3);
    CNot_13 = QCO.get_full_sized_gate("CNot_1_3", 3);

    @test isapprox(CV_32 * CNot_13 * CVdagger_32 * CNot_13 * CV_12, CVdagger_32 * CNot_13 * CV_32 * CNot_13 * CVdagger_12, atol = tol_0)

    CVdagger_13 = QCO.get_full_sized_gate("CVdagger_1_3", 3);
    cnot_21 = QCO.get_full_sized_gate("CNot_2_1", 3);
    @test isapprox(CV_13 * cnot_21 * CVdagger_13 * cnot_21 * CV_23, CVdagger_13 * cnot_21 * CV_13 * cnot_21 * CVdagger_23, atol = tol_0)

    # 2-qubit Grover's diffusion operator (Ref: https://arxiv.org/pdf/1804.03719.pdf)
    H_1 = QCO.get_full_sized_gate("H_1", 2);
    X_2 = QCO.get_full_sized_gate("X_2", 2);
    CNot_1_2 = QCO.get_full_sized_gate("CNot_1_2", 2);

    @test isapprox(QCO.GroverDiffusionGate(), Y_2 * H_1 * X_2 * CNot_1_2 * H_1 * Y_2, atol=tol_0)

    # 2 Qubit QFT (Ref: https://www.cs.bham.ac.uk/internal/courses/intro-mqc/current/lecture06_handout.pdf)

    SWAP = QCO.get_full_sized_gate("Swap_1_2", 2);
    CU = QCO.get_full_sized_gate("CU3_2_1", 2, angle = [0, π/4, π/4]);

    @test isapprox(QCO.QFT2Gate(), H_1 * CU * H_2 * SWAP)

    # CRY Decomp (Ref: https://quantumcomputing.stackexchange.com/questions/2143/how-can-a-controlled-ry-be-made-from-cnots-and-rotations)

    CNOT12 = QCO.get_full_sized_gate("CNot_1_2", 2);
    U3_negpi = QCO.get_full_sized_gate("U3_2", 2, angle=[-pi/2, 0, 0]);
    U3_pi = QCO.get_full_sized_gate("U3_2", 2, angle=[pi/2, 0, 0]);

    @test isapprox(QCO.CRYGate(pi), CNOT12 * U3_negpi * CNOT12 * U3_pi, atol = tol_0)

    # Identity tests

    Id = QCO.IGate(2);
    CRXRev = QCO.get_full_sized_gate("CRX_2_1", 2, angle=pi);

    @test isapprox(Id, CRXRev * CRXRev * CRXRev * CRXRev, atol = tol_0)

    CRYRev = QCO.get_full_sized_gate("CRY_2_1", 2, angle=pi);
    @test isapprox(Id, CRYRev * CRYRev * CRYRev * CRYRev, atol = tol_0)

    CRZRev = QCO.get_full_sized_gate("CRZ_2_1", 2, angle=pi);
    @test isapprox(Id, CRZRev * CRZRev * CRZRev * CRZRev, atol = tol_0)

    # Rev gate tests for involution
    @test isapprox(QCO.CXRevGate(), QCO.CNotRevGate(), atol = tol_0)
    @test isapprox(QCO.CXGate() * QCO.SwapGate() * QCO.CXRevGate() * QCO.SwapGate(), QCO.IGate(2), atol = tol_0)
    @test isapprox(QCO.CYGate() * QCO.SwapGate() * QCO.CYRevGate() * QCO.SwapGate(), QCO.IGate(2), atol = tol_0)
    @test isapprox(QCO.CZGate() * QCO.SwapGate() * QCO.CZRevGate() * QCO.SwapGate(), QCO.IGate(2), atol = tol_0)
    @test isapprox(QCO.CHGate() * QCO.SwapGate() * QCO.CHRevGate() * QCO.SwapGate(), QCO.IGate(2), atol = tol_0)
    
    # Rev gate tests for non-involution
    @test isapprox(QCO.CVGate() * QCO.SwapGate() * QCO.CVRevGate() * QCO.SwapGate(), QCO.CVGate()^2, atol = tol_0)
    @test isapprox(QCO.CSXGate() * QCO.SwapGate() * QCO.CSXRevGate() * QCO.SwapGate(), QCO.CSXGate()^2, atol = tol_0)

    @test isapprox(QCO.WGate(), QCO.HCoinGate(), atol = tol_0)

    # Fredkin test = CSwapGate
    CNot_32 = QCO.get_full_sized_gate("CNot_3_2", 3)
    CV_23 = QCO.get_full_sized_gate("CV_2_3", 3)
    CV_13 = QCO.get_full_sized_gate("CV_1_3", 3)
    CNot_12 = QCO.get_full_sized_gate("CNot_1_2", 3)
    CVdagger_23 = QCO.get_full_sized_gate("CVdagger_2_3", 3)
    @test isapprox(QCO.CSwapGate(), CNot_32 * CV_23 * CV_13 * CNot_12 * CVdagger_23 * CNot_12 * CNot_32, atol = tol_0)

    # CCZGate test: CCZ is equivalent to Toffoli when the target qubit is conjugated by Hadamard gates
    H_3 = QCO.get_full_sized_gate("H_3", 3)
    @test isapprox(QCO.ToffoliGate(), H_3 * QCO.CCZGate() * H_3, atol = tol_0)

    # Peres test 
    CVdagger_13 = QCO.get_full_sized_gate("CVdagger_1_3", 3)
    @test isapprox(QCO.PeresGate(), QCO.ToffoliGate() * CNot_12, atol = tol_0)
    @test isapprox(QCO.PeresGate(), CVdagger_13 * CVdagger_23 * CNot_12 * CV_23, atol = tol_0)

    # iSwap test 
    @test isapprox(QCO.get_full_sized_gate("iSwap_2_1", 3), QCO.get_full_sized_gate("iSwap_1_2", 3), atol = tol_0)
    @test isapprox(QCO.get_full_sized_gate("iSwap_2_3", 4), QCO.get_full_sized_gate("iSwap_3_2", 4), atol = tol_0)
    @test isapprox(QCO.get_full_sized_gate("iSwap_3_5", 5), QCO.get_full_sized_gate("iSwap_5_3", 5), atol = tol_0)

    # Sycamore test
    @test isapprox(QCO.SycamoreGate()^12, QCO.IGate(2), atol = tol_0)

    # RCCX test
    T_3 = QCO.get_full_sized_gate("T_3", 3)
    Tdagger_3 = QCO.get_full_sized_gate("Tdagger_3", 3)
    @test isapprox(QCO.RCCXGate(),  H_3 * Tdagger_3 * CNot_23 * T_3 * CNot_13 * Tdagger_3 * CNot_23 * T_3 * H_3, atol = tol_0)

    # Margolus test 
    RY_a = QCO.get_full_sized_gate("RY_3", 3, angle = -π/4)
    RY_b = QCO.get_full_sized_gate("RY_3", 3, angle = π/4)
    @test isapprox(QCO.MargolusGate(),  RY_a * CNot_23 * RY_a * CNot_13 * RY_b * CNot_23 * RY_b, atol = tol_0)
end

