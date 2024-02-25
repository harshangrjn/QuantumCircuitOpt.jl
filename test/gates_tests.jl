# Unit tests for functions in gates.jl

@testset "Gates tests: elementary universal gates in gates.jl" begin
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
    H2    = kron(QCO.IGate(1), QCO.HGate())
    @test isapprox(QCO.MGate(), QCO.CNotRevGate() * H2 * S1_S2)

    Z1        = QCO.get_unitary("Z_1", 2);
    Z2        = QCO.get_unitary("Z_2", 2);
    T2        = QCO.get_unitary("T_2", 2);
    Y_2       = QCO.get_unitary("Y_2", 2);
    CNot_1_2  = QCO.get_unitary("CNot_1_2", 2);
    CNot_2_1  = QCO.get_unitary("CNot_2_1", 2);
    Sdagger1  = QCO.get_unitary("Sdagger_1", 2);
    Tdagger1  = QCO.get_unitary("Tdagger_1", 2);
    SX1       = QCO.get_unitary("SX_1", 2);
    SXdagger2 = QCO.get_unitary("SXdagger_2", 2);
    @test isapprox(-QCO.HCoinGate(), Z2 * CNot_2_1 * SXdagger2 * Y_2 * CNot_2_1 * Tdagger1 * Z1 * T2 * Sdagger1 * CNot_1_2 * Sdagger1 * SX1 * CNot_2_1 * CNot_1_2, atol = tol_0)

    CV_23       = QCO.get_unitary("CV_2_3", 3);
    CNot_1_2    = QCO.get_unitary("CNot_1_2", 3);
    CVdagger_23 = QCO.get_unitary("CVdagger_2_3", 3);
    CV_13       = QCO.get_unitary("CV_1_3", 3);
    @test isapprox(QCO.ToffoliGate(), CV_23 * CNot_1_2 * CVdagger_23 * CNot_1_2 * CV_13, atol = tol_0)

    CVdagger_12 = QCO.get_unitary("CVdagger_1_2", 3);
    CVdagger_31 = QCO.get_unitary("CVdagger_3_1", 3);
    @test isapprox(QCO.IGate(3), CVdagger_12 * CVdagger_12 * CVdagger_12 * CVdagger_12 * CVdagger_31 * CVdagger_31 * CVdagger_31 * CVdagger_31, atol = tol_0)

    H_2      = QCO.get_unitary("H_2", 2);
    T_1      = QCO.get_unitary("T_1", 2);
    Tdagger2 = QCO.get_unitary("Tdagger_2", 2);
    @test isapprox(QCO.CVdaggerGate(), H_2 * Tdagger1 * CNot_2_1 * T_1 * Tdagger2 * CNot_2_1 * H_2, atol = tol_0)

    @test isapprox(QCO.CVGate(), QCO.CNotGate() * QCO.CVdaggerGate(), atol = tol_0)

    # Next 3 tests from: https://american-cse.org/csci2015/data/9795a059.pdf

    CV_31       = QCO.get_unitary("CV_3_1", 3);
    CNot_2_3    = QCO.get_unitary("CNot_2_3", 3);
    CVdagger_31 = QCO.get_unitary("CVdagger_3_1", 3);
    CV_21       = QCO.get_unitary("CV_2_1", 3);

    CVdagger_21 = QCO.get_unitary("CVdagger_2_1", 3);

    @test isapprox(CV_31 * CNot_2_3 * CVdagger_31 * CNot_2_3 * CV_21, CVdagger_21 * CNot_2_3 * CV_31 * CNot_2_3 * CVdagger_31, atol = tol_0)

    CV_12       = QCO.get_unitary("CV_1_2", 3);
    CV_32       = QCO.get_unitary("CV_3_2", 3);
    CVdagger_32 = QCO.get_unitary("CVdagger_3_2", 3);
    CVdagger_12 = QCO.get_unitary("CVdagger_1_2", 3);
    CNot_1_3    = QCO.get_unitary("CNot_1_3", 3);

    @test isapprox(CV_32 * CNot_1_3 * CVdagger_32 * CNot_1_3 * CV_12, CVdagger_32 * CNot_1_3 * CV_32 * CNot_1_3 * CVdagger_12, atol = tol_0)

    CVdagger_13  = QCO.get_unitary("CVdagger_1_3", 3);
    CNot_2_1     = QCO.get_unitary("CNot_2_1", 3);
    @test isapprox(CV_13 * CNot_2_1 * CVdagger_13 * CNot_2_1 * CV_23, CVdagger_13 * CNot_2_1 * CV_13 * CNot_2_1 * CVdagger_23, atol = tol_0)

    # 2-qubit Grover's diffusion operator (Ref: https://arxiv.org/pdf/1804.03719.pdf)
    H_1      = QCO.get_unitary("H_1", 2);
    X_2      = QCO.get_unitary("X_2", 2);
    CNot_1_2 = QCO.get_unitary("CNot_1_2", 2);

    @test isapprox(QCO.GroverDiffusionGate(), Y_2 * H_1 * X_2 * CNot_1_2 * H_1 * Y_2, atol=tol_0)

    # 2 Qubit QFT (Ref: https://www.cs.bham.ac.uk/internal/courses/intro-mqc/current/lecture06_handout.pdf)

    SWAP = QCO.get_unitary("Swap_1_2", 2);
    CU   = QCO.get_unitary("CU3_2_1", 2, angle = [0, π/4, π/4]);

    @test isapprox(QCO.QFT2Gate(), H_1 * CU * H_2 * SWAP)

    # 3 Qubit QFT (Ref: Nielsen and Chuang, Quantum Computation and Quantum Information, Ch. 5.1)
    CS_2_1 = QCO.get_unitary("CS_2_1", 3)
    CS_3_2 = QCO.get_unitary("CS_3_2", 3)
    CT_3_1 = QCO.get_unitary("CT_3_1", 3)
    SWAP_1_3 = QCO.get_unitary("Swap_1_3", 3)
    H_1 = QCO.get_unitary("H_1", 3)
    H_2 = QCO.get_unitary("H_2", 3)
    H_3 = QCO.get_unitary("H_3", 3)
    @test isapprox(QCO.QFT3Gate(), H_1 * CS_2_1 * CT_3_1 * H_2 * CS_3_2 * H_3 * SWAP_1_3)

    # CRY Decomp (Ref: https://quantumcomputing.stackexchange.com/questions/2143/how-can-a-controlled-ry-be-made-from-cnots-and-rotations)

    CNOT12   = QCO.get_unitary("CNot_1_2", 2);
    U3_negpi = QCO.get_unitary("U3_2", 2, angle=[-pi/2, 0, 0]);
    U3_pi    = QCO.get_unitary("U3_2", 2, angle=[pi/2, 0, 0]);

    @test isapprox(QCO.CRYGate(pi), CNOT12 * U3_negpi * CNOT12 * U3_pi, atol = tol_0)

    # Identity tests

    Id     = QCO.IGate(2);
    CRXRev = QCO.get_unitary("CRX_2_1", 2, angle=pi);

    @test isapprox(Id, CRXRev * CRXRev * CRXRev * CRXRev, atol = tol_0)

    CRYRev = QCO.get_unitary("CRY_2_1", 2, angle=pi);
    @test isapprox(Id, CRYRev * CRYRev * CRYRev * CRYRev, atol = tol_0)

    CRZRev = QCO.get_unitary("CRZ_2_1", 2, angle=pi);
    @test isapprox(Id, CRZRev * CRZRev * CRZRev * CRZRev, atol = tol_0)

    # Rev gate tests for involution
    @test isapprox(QCO.CXRevGate(), QCO.CNotRevGate(), atol = tol_0)
    @test isapprox(QCO.CXGate() * QCO.SwapGate() * QCO.CXRevGate() * QCO.SwapGate(), QCO.IGate(2), atol = tol_0)
    @test isapprox(QCO.CYGate() * QCO.SwapGate() * QCO.CYRevGate() * QCO.SwapGate(), QCO.IGate(2), atol = tol_0)
    @test isapprox(QCO.CZGate() * QCO.SwapGate() * QCO.CZGate()    * QCO.SwapGate(), QCO.IGate(2), atol = tol_0)
    @test isapprox(QCO.CHGate() * QCO.SwapGate() * QCO.CHRevGate() * QCO.SwapGate(), QCO.IGate(2), atol = tol_0)
    
    # Rev gate tests for non-involution
    @test isapprox(QCO.CVGate() * QCO.SwapGate() * QCO.CVRevGate() * QCO.SwapGate(), QCO.CVGate()^2, atol = tol_0)
    @test isapprox(QCO.CSXGate() * QCO.SwapGate() * QCO.CSXRevGate() * QCO.SwapGate(), QCO.CSXGate()^2, atol = tol_0)

    @test isapprox(QCO.WGate(), QCO.HCoinGate(), atol = tol_0)

    # Fredkin test = CSwapGate
    CNot_3_2    = QCO.get_unitary("CNot_3_2", 3)
    CV_23       = QCO.get_unitary("CV_2_3", 3)
    CV_13       = QCO.get_unitary("CV_1_3", 3)
    CNot_1_2    = QCO.get_unitary("CNot_1_2", 3)
    CVdagger_23 = QCO.get_unitary("CVdagger_2_3", 3)
    @test isapprox(QCO.CSwapGate(), CNot_3_2 * CV_23 * CV_13 * CNot_1_2 * CVdagger_23 * CNot_1_2 * CNot_3_2, atol = tol_0)

    # CCZGate test: CCZ is equivalent to Toffoli when the target qubit is conjugated by Hadamard gates
    H_3 = QCO.get_unitary("H_3", 3)
    @test isapprox(QCO.ToffoliGate(), H_3 * QCO.CCZGate() * H_3, atol = tol_0)

    # Peres test 
    CVdagger_13 = QCO.get_unitary("CVdagger_1_3", 3)
    @test isapprox(QCO.PeresGate(), QCO.ToffoliGate() * CNot_1_2, atol = tol_0)
    @test isapprox(QCO.PeresGate(), CVdagger_13 * CVdagger_23 * CNot_1_2 * CV_23, atol = tol_0)

    # iSwap test 
    @test isapprox(QCO.get_unitary("iSwap_2_1", 3), QCO.get_unitary("iSwap_1_2", 3), atol = tol_0)
    @test isapprox(QCO.get_unitary("iSwap_2_3", 4), QCO.get_unitary("iSwap_3_2", 4), atol = tol_0)
    @test isapprox(QCO.get_unitary("iSwap_3_5", 5), QCO.get_unitary("iSwap_5_3", 5), atol = tol_0)

    # Sycamore test
    @test isapprox(QCO.SycamoreGate()^12, QCO.IGate(2), atol = tol_0)

    # RCCX test
    T_3       = QCO.get_unitary("T_3", 3)
    Tdagger_3 = QCO.get_unitary("Tdagger_3", 3)
    @test isapprox(QCO.RCCXGate(),  H_3 * Tdagger_3 * CNot_2_3 * T_3 * CNot_1_3 * Tdagger_3 * CNot_2_3 * T_3 * H_3, atol = tol_0)

    # Margolus test 
    RY_a = QCO.get_unitary("RY_3", 3, angle = -π/4)
    RY_b = QCO.get_unitary("RY_3", 3, angle = π/4)
    @test isapprox(QCO.MargolusGate(),  RY_a * CNot_2_3 * RY_a * CNot_1_3 * RY_b * CNot_2_3 * RY_b, atol = tol_0)
    
    # CiSwap test 
    @test isapprox(QCO.CiSwapGate(),  CNot_3_2 * H_3 * CNot_2_3 * T_3 * CNot_1_3 * Tdagger_3 * CNot_2_3 * T_3 * CNot_1_3 * Tdagger_3 * H_3 * CNot_3_2, atol = tol_0)

    # 2-angle RGate test
    # Ref: https://qiskit.org/documentation/stubs/qiskit.circuit.library.RGate.html#qiskit.circuit.library.RGate
    θ = π/3; ϕ = π/6;
    R1 = QCO.RGate(θ, ϕ)
    U3 = QCO.U3Gate(θ, ϕ, -ϕ)
    U3[1,2] = U3[1,2] * im; U3[2,1] = U3[2,1] * -im;
    @test isapprox(R1, U3; atol = tol_0)

    # CS, CSdagger gate test
    H1    = QCO.get_unitary("H_1", 2)
    H2    = QCO.get_unitary("H_2", 2)
    CS_12 = QCO.get_unitary("CS_1_2", 2)
    @test isapprox(H1 * CS_12 * H2 * QCO.SwapGate(), QCO.QFT2Gate(), atol = tol_0)
    @test isapprox(QCO.CSGate() *  QCO.CSdaggerGate(), QCO.IGate(2), atol = tol_0)

    #SSwapGate test
    @test isapprox(QCO.SSwapGate()^2, QCO.SwapGate(), atol = 1E-6)

    # CT, CTdagger gate test
    @test isapprox(QCO.CTGate() *  QCO.CTdaggerGate(), QCO.IGate(2), atol = tol_0)

    # Test for elementary gates with kron symbols
    M1 = QCO.get_unitary("I_1xCZ_2_4xH_5", 5)
    M2 = kron(kron(QCO.IGate(1), QCO.get_unitary("CZ_1_3", 3)), QCO.HGate())
    @test isapprox(M1, M2, atol = tol_0)

    M1 = QCO.get_unitary("CV_4_1xH_5xZ_6", 6)
    M2 = kron(kron(QCO.get_unitary("CV_4_1", 4), QCO.HGate()), QCO.ZGate())
    @test isapprox(M1, M2, atol = tol_0)

end