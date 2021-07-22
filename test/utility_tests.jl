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

@testset "complex to real matrix function tests" begin
    test_M_r = QCO.complex_to_real_matrix(M_c)
    @test test_M_r == M_r
end

@testset "real to complex matrix function tests" begin
    test_M_c = QCO.real_to_complex_matrix(M_r)
    @test test_M_c == M_c
end

@testset "auxiliary variable product bounds tests" begin
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

@testset "unique_matrices and unique_idx tests" begin
    
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

@testset "commuting matrices tests" begin
    params = Dict{String, Any}(
       "num_qubits" => 2, 
       "depth" => 5,    

       "elementary_gates" => ["H1", "H2", "cnot_12", "Identity"],  
       "target_gate" => QCO.CNotRevGate()
       )
    
    data = QCO.get_data(params)
    C2,C2_Identity = QCO.get_commutative_gate_pairs(data["gates_dict"])
    @test length(C2) == 4
    @test length(C2_Identity) == 0
    @test isapprox(data["gates_real"][:,:,C2[1][1]] * data["gates_real"][:,:,C2[1][2]], data["gates_real"][:,:,C2[1][2]] * data["gates_real"][:,:,C2[1][1]], atol=tol_0)
end

@testset "kron_single_gate tests" begin
    I1 = QCO.kron_single_gate(4, QCO.IGate(1), "q1")
    I2 = QCO.kron_single_gate(4, QCO.IGate(1), "q2")
    I3 = QCO.kron_single_gate(4, QCO.IGate(1), "q3")
    I4 = QCO.kron_single_gate(4, QCO.IGate(1), "q4")
    @test I1 == I2 == I3 == I4
end

@testset "get_redundant_gate_product_pairs tests" begin 
     
    params = Dict{String, Any}(
        "num_qubits" => 3, 
        "depth" => 15,    
        "elementary_gates" => ["T1", "T2", "T3", "H3", "cnot_12", "cnot_13", "cnot_23", "Tdagger2", "Tdagger3", "Identity"], 
        "target_gate" => QCO.ToffoliGate()
        )
    data = QCO.get_data(params)
    redundant_pairs = QCO.get_redundant_gate_product_pairs(data["gates_dict"])
    @test length(redundant_pairs) == 0

end

@testset "get_idempotent_gates tests" begin
    params = Dict{String, Any}(
        "num_qubits" => 2, 
        "depth" => 14,    
        "elementary_gates" => ["Y1", "Y2", "Z1", "Z2", "T2", "Tdagger1", "Sdagger1", "SX1", "SXdagger2", "cnot_21", "cnot_12", "Identity"], 
        "target_gate" => -QCO.HCoinGate())
    data = QCO.get_data(params)
    idempotent_gates = QCO.get_idempotent_gates(data["gates_dict"])
    @test idempotent_gates[1] == 6
    @test idempotent_gates[2] == 7
    @test "Tdagger1" in data["gates_dict"]["$(idempotent_gates[1])"]["type"]
    @test "Sdagger1" in data["gates_dict"]["$(idempotent_gates[2])"]["type"]
    @test "Z1" in data["gates_dict"]["3"]["type"]
    
    M1 = data["gates_dict"]["$(idempotent_gates[1])"]["matrix"]
    M2 = data["gates_dict"]["$(idempotent_gates[2])"]["matrix"]
    
    @test isapprox(M1^2, data["gates_dict"]["7"]["matrix"], atol = tol_0)
    @test isapprox(M2^2, data["gates_dict"]["3"]["matrix"], atol = tol_0)
end