# Unit tests for functions in utility.jl

M_c = [complex(1,0)       complex(0,-1)
       complex(-0.5,0.5)  complex(0,0)]

M_r = [1.0   0.0   0.0  -1.0
       0.0   1.0   1.0  0.0
      -0.5  0.5    0.0  0.0
      -0.5  -0.5   0.0  0.0]

# Variable bounds to test auxiliary_variable_bounds function
l = [-2, -1, 0, 2]
u = [-1,  2, 3, 2]

@testset "Tests: complex to real matrix function" begin
    test_M_r = QCO.complex_to_real_gate(M_c)
    @test test_M_r == M_r
end

@testset "Tests: real to complex matrix function" begin
    test_M_c = QCO.real_to_complex_gate(M_r)
    @test test_M_c == M_c
end

@testset "Tests: auxiliary variable product bounds" begin
    m = Model()
    @variable(m, l[i] <= x[i=1:length(l)] <= u[i])
    # x1*x2 
    @test QCO.auxiliary_variable_bounds([x[1], x[2]]) == (-4, 2)
    # x1*x3 
    @test QCO.auxiliary_variable_bounds([x[1], x[3]]) == (-6, 0)
    # x1*x4
    @test QCO.auxiliary_variable_bounds([x[1], x[4]]) == (-4, -2)
    # x4*x4
    @test QCO.auxiliary_variable_bounds([x[4], x[4]]) == (4, 4)
    # x1*x2*x3 
    @test QCO.auxiliary_variable_bounds([x[1], x[2], x[3]]) == (-12, 6)
    @test QCO.auxiliary_variable_bounds([x[3], x[1], x[2]]) == (-12, 6)
    @test QCO.auxiliary_variable_bounds([x[2], x[1], x[3]]) == (-12, 6)
    # x1*x2*x3*x4
    @test QCO.auxiliary_variable_bounds([x[1], x[2], x[3], x[4]]) == (-24, 12)
    @test QCO.auxiliary_variable_bounds([x[4], x[1], x[3], x[2]]) == (-24, 12)
    @test QCO.auxiliary_variable_bounds([x[4], x[3], x[2], x[1]]) == (-24, 12)
end

@testset "Tests: unique_matrices and unique_idx" begin
    
    D = zeros(3,3,4)
    D[:,:,1] = Matrix(LA.I,3,3)
    D[:,:,2] = rand(3,3)
    D[:,:,3] = Matrix(LA.I,3,3)
    D[:,:,4] = rand(3,3)

    D[1,1,3] = 1 + 2*tol_0
    D[3,3,3] = 1 - tol_0
    D[1,2,3] = 0 - tol_0
    D[1,3,3] = 0 + tol_0

    D[isapprox.(D, 0, atol=tol_0)] .= 0

    D_unique, D_unique_idx = QCO.unique_matrices(D)

    @test D_unique_idx == [1,2,4]
    @test D_unique[:,:,1] == D[:,:,1]
    @test D_unique[:,:,2] == D[:,:,2]
    @test D_unique[:,:,3] == D[:,:,4]

end

@testset "Tests: commuting matrices" begin
    params = Dict{String, Any}(
       "num_qubits" => 2, 
       "depth" => 5,    

       "elementary_gates" => ["H_1", "H_2", "CNot_1_2", "Identity"],  
       "target_gate" => QCO.CNotRevGate()
       )
    
    data = QCO.get_data(params)
    C2,C2_Identity = QCO.get_commutative_gate_pairs(data["gates_dict"])
    @test length(C2) == 4
    @test length(C2_Identity) == 0
    @test isapprox(data["gates_real"][:,:,C2[1][1]] * data["gates_real"][:,:,C2[1][2]], data["gates_real"][:,:,C2[1][2]] * data["gates_real"][:,:,C2[1][1]], atol=tol_0)
end

@testset "Tests: kron_single_qubit_gate" begin
    I1 = QCO.kron_single_qubit_gate(4, QCO.IGate(1), "q1")
    I2 = QCO.kron_single_qubit_gate(4, QCO.IGate(1), "q2")
    I3 = QCO.kron_single_qubit_gate(4, QCO.IGate(1), "q3")
    I4 = QCO.kron_single_qubit_gate(4, QCO.IGate(1), "q4")
    @test I1 == I2 == I3 == I4
end

@testset "Tests: get_redundant_gate_product_pairs" begin 
     
    params = Dict{String, Any}(
        "num_qubits" => 3, 
        "depth" => 15,    
        "elementary_gates" => ["T_1", "T_2", "T_3", "H_3", "CNot_1_2", "CNot_1_3", "CNot_2_3", "Tdagger_2", "Tdagger_3", "Identity"], 
        "target_gate" => QCO.ToffoliGate()
        )
    data = QCO.get_data(params)
    redundant_pairs = QCO.get_redundant_gate_product_pairs(data["gates_dict"])
    @test length(redundant_pairs) == 0

end

@testset "Tests: get_idempotent_gates" begin
    
    params = Dict{String, Any}(
        "num_qubits" => 2, 
        "depth" => 14,    
        "elementary_gates" => ["Y_1", "Y_2", "Z_1", "Z_2", "T_2", "Tdagger_1", "Sdagger_1", "SX_1", "SXdagger_2", "CNot_2_1", "CNot_1_2", "Identity"], 
        "target_gate" => -QCO.HCoinGate())

    data = QCO.get_data(params)
    idempotent_gates = QCO.get_idempotent_gates(data["gates_dict"])
    @test idempotent_gates[1] == 6
    @test idempotent_gates[2] == 7
    @test "Tdagger_1" in data["gates_dict"]["$(idempotent_gates[1])"]["type"]
    @test "Sdagger_1" in data["gates_dict"]["$(idempotent_gates[2])"]["type"]
    @test "Z_1" in data["gates_dict"]["3"]["type"]
    
    M1 = data["gates_dict"]["$(idempotent_gates[1])"]["matrix"]
    M2 = data["gates_dict"]["$(idempotent_gates[2])"]["matrix"]
    
    @test isapprox(M1^2, data["gates_dict"]["7"]["matrix"], atol = tol_0)
    @test isapprox(M2^2, data["gates_dict"]["3"]["matrix"], atol = tol_0)
end

@testset "Tests: _parse_gate_string" begin
    v1 = QCO._parse_gate_string("RX_1", qubits=true)    
    @test v1[1] == 1
    
    v2 = QCO._parse_gate_string("CU3_51", qubits=true)    
    @test v2[1] == 51

    v3 = QCO._parse_gate_string("CRZ_9_1", qubits=true)    
    @test length(v3) == 2
    @test v3[1] == 9
    @test v3[2] == 1
end

@testset "Tests: is_gate_real" begin
    M_test = [0 1E-7im -1E-7im; 1E-10im 0 0; -1E-8im -1E-10im 0]
    @test QCO.is_gate_real(M_test)
    
    M_test = [0 1E-7im -1E-7im; 1E-10im 0 0; -1E-8im -1E-5im 0]
    @test !QCO.is_gate_real(M_test)
end

@testset "Tests: _get_constraint_slope_intercept" begin
    v1 = [0.0, 0.0]
    v2 = [1.0, 1.0]
    m,c = QCO._get_constraint_slope_intercept(v1, v2)
    @test m == 1 
    @test c == 0 

    m,c = QCO._get_constraint_slope_intercept(v2, v1)
    @test m == 1 
    @test c == 0 

    v2 = [0, 0.5]
    m,c = QCO._get_constraint_slope_intercept(v1, v2)
    @test !isfinite(m)
    @test !isfinite(c)

    v1 = [1.0, 2.0]
    v2 = [-3.0, 4.0]
    m,c = QCO._get_constraint_slope_intercept(v1, v2)
    @test m == -0.5 
    @test c == 2.5

    QCO._get_constraint_slope_intercept(v1, v1)
end