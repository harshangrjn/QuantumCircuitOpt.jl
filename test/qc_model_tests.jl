@testset "Tests: Minimize depth for controlled-NOT gate" begin

    params = Dict{String, Any}(
    "num_qubits" => 2, 
    "depth" => 5,
    "elementary_gates" => ["H_1", "H_2", "CNot_1_2", "Identity"],  
    "initial_gate" => "Identity",
    "identify_real_gates" => true,
    "target_gate" => QCO.CNotRevGate(),
    "set_cnot_lower_bound" => 1,
    "set_cnot_upper_bound" => 1,
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact"
    )

    result_qc = QCO.run_QCModel(params, CBC, model_type = "compact_formulation_1", all_valid_constraints = 2)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 5, atol=tol_0)
    @test isapprox(result_qc["solution"]["z_onoff_var"][1,1], 1, atol=tol_0) || isapprox(result_qc["solution"]["z_onoff_var"][2,1], 1, atol=tol_0)
    @test isapprox(result_qc["solution"]["z_onoff_var"][1,2], 1, atol=tol_0) || isapprox(result_qc["solution"]["z_onoff_var"][2,2], 1, atol=tol_0)
    @test isapprox(result_qc["solution"]["z_onoff_var"][3,3], 1, atol=tol_0)
    @test isapprox(result_qc["solution"]["z_onoff_var"][1,4], 1, atol=tol_0) || isapprox(result_qc["solution"]["z_onoff_var"][2,4], 1, atol=tol_0)
    @test isapprox(result_qc["solution"]["z_onoff_var"][1,5], 1, atol=tol_0) || isapprox(result_qc["solution"]["z_onoff_var"][2,5], 1, atol=tol_0)
    
end

@testset "Tests: Minimum CNOT swap gate decomposition" begin

    params = Dict{String, Any}(
        "num_qubits" => 2,
        "depth" => 4,        
        "elementary_gates" => ["CNot_1_2", "CNot_2_1", "Identity"],
        "identify_real_gates" => true,
        "target_gate" => QCO.SwapGate(),  
        "objective" => "minimize_cnot", 
        "decomposition_type" => "exact"                      
        )

    result_qc = QCO.run_QCModel(params, CBC, model_type = "balas_formulation")

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 3, atol = tol_0)
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][1:2,:]), 3, atol=tol_0)
    
end

@testset "Tests: Minimum depth U3(0,0,π/4) gate" begin

    params = Dict{String, Any}(
    "num_qubits" => 2, 
    "depth" => 2,    
    "elementary_gates" => ["U3_1", "U3_2", "Identity"],  
    "target_gate" => QCO.kron_single_qubit_gate(2, QCO.U3Gate(0,0,π/4), "q1"),
    "U3_θ_discretization" => [0, π/2],
    "U3_ϕ_discretization" => [0],
    "U3_λ_discretization" => [0, π/4],
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact"                  
    )

    result_qc = QCO.run_QCModel(params, CBC, model_type = "compact_formulation")

    data = QCO.get_data(params)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 1, atol = tol_0)
    @test "Identity" in data["gates_dict"]["3"]["type"]
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][3,:]), 1, atol = tol_0)
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][4,:]), 1, atol = tol_0) 
    @test data["gates_dict"]["4"]["qubit_loc"] == "qubit_1"
end

@testset "Tests: Minimum depth CU3(0,0,π/4) gate" begin

    params = Dict{String, Any}(
    "num_qubits" => 2, 
    "depth" => 2,    
    "elementary_gates" => ["CU3_1_2", "CU3_2_1", "Identity"],  
    "target_gate" => QCO.CU3Gate(0, 0, π/4),
    "CU3_θ_discretization" => [0, π/2],
    "CU3_ϕ_discretization" => [0],
    "CU3_λ_discretization" => [0, π/4],
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact"                  
    )

    result_qc = QCO.run_QCModel(params, CBC, model_type = "compact_formulation")

    data = QCO.get_data(params)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 1, atol = tol_0)
    @test "Identity" in data["gates_dict"]["3"]["type"]
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][3,:]), 1, atol = tol_0)
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][4,:]), 1, atol = tol_0) 
    @test data["gates_dict"]["4"]["qubit_loc"] == "qubit_1_2"
end

@testset "Tests: Minimum depth Rev CU3(0,0,π/4) gate" begin

    params = Dict{String, Any}(
    "num_qubits" => 3, 
    "depth" => 2,    
    "elementary_gates" => ["CU3_3_1", "CU3_1_3", "Identity"],  
    "target_gate" => QCO.get_full_sized_gate("CU3_3_1", 3, angle = [0, 0, pi/4]),
    "CU3_θ_discretization" => [0, π/2],
    "CU3_ϕ_discretization" => [0],
    "CU3_λ_discretization" => [0, π/4],
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact"                  
    )

    result_qc = QCO.run_QCModel(params, CBC, model_type = "compact_formulation")

    data = QCO.get_data(params)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 1, atol = tol_0)
    @test "Identity" in data["gates_dict"]["3"]["type"]
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][3,:]), 1, atol = tol_0)
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][4,:]), 1, atol = tol_0) 
    @test data["gates_dict"]["4"]["qubit_loc"] == "qubit_3_1"
end

@testset "Tests: Minimum depth RX, RY, RZ gate decomposition" begin

    params = Dict{String, Any}(
    "num_qubits" => 2, 
    "depth" => 3,
    "elementary_gates" => ["RX_1", "RY_2", "RZ_1", "Identity"],  
    "target_gate" => QCO.kron_single_qubit_gate(2, QCO.RXGate(π/4), "q1") * QCO.kron_single_qubit_gate(2, QCO.RYGate(π/4), "q2") * QCO.kron_single_qubit_gate(2, QCO.RZGate(π/4), "q1"),
    "RX_discretization" => [0, π/4],
    "RY_discretization" => [π/4],
    "RZ_discretization" => [π/2, π/4],
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact"                            
    )

    result_qc = QCO.run_QCModel(params, CBC, model_type = "balas_formulation", commute_gate_constraints = true)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 3, atol = tol_0)

end

@testset "Tests: Minimum depth CRX, CRY, CRZ gate decomposition" begin

    params = Dict{String, Any}(
    "num_qubits" => 3, 
    "depth" => 3,
    "elementary_gates" => ["CRX_1_2", "CRY_2_3", "CRZ_3_1", "Identity"],  
    "target_gate" => QCO.get_full_sized_gate("CRX_1_2", 3, angle = pi/4) * QCO.get_full_sized_gate("CRY_2_3", 3, angle=pi/4) * QCO.get_full_sized_gate("CRZ_3_1", 3, angle=pi/4),
    "CRX_discretization" => [0, π/4],
    "CRY_discretization" => [π/4],
    "CRZ_discretization" => [π/2, π/4],
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact"                            
    )

    result_qc = QCO.run_QCModel(params, CBC, model_type = "balas_formulation", commute_gate_constraints = true)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 3, atol = tol_0)

end

@testset "Tests: Minimum depth Rev CRX, CRY, CRZ gate decomposition" begin

    params = Dict{String, Any}(
    "num_qubits" => 3, 
    "depth" => 3,
    "elementary_gates" => ["CRX_3_1", "CRY_3_1", "CRZ_1_3", "Identity"],  
    "target_gate" => QCO.get_full_sized_gate("CRX_3_1", 3, angle=pi/4) * QCO.get_full_sized_gate("CRY_3_1", 3, angle=pi/4) * QCO.get_full_sized_gate("CRZ_1_3", 3, angle = pi/2),
    "CRX_discretization" => [0, π/4],
    "CRY_discretization" => [π/4],
    "CRZ_discretization" => [π/2, π/4],
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact"                            
    )

    result_qc = QCO.run_QCModel(params, CBC, model_type = "balas_formulation", commute_gate_constraints = true)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 3, atol = tol_0)

end


@testset "Tests: 3-qubit RX gate decomposition" begin

    params = Dict{String, Any}(
    "num_qubits" => 3, 
    "depth" => 2,    
    "elementary_gates" => ["U3_1", "U3_2", "U3_3", "Identity"],  
    "target_gate" => QCO.kron_single_qubit_gate(3, QCO.RXGate(π/4), "q3"),
    "U3_θ_discretization" => [0, π/4],
    "U3_ϕ_discretization" => [0, -π/2],
    "U3_λ_discretization" => [0, π/2],    
    "objective" => "minimize_cnot", 
    "decomposition_type" => "exact"                   
    )

    result_qc = QCO.run_QCModel(params, CBC, model_type = "balas_formulation")
    
    data = QCO.get_data(params)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 1, atol = tol_0)
    if isapprox(sum(result_qc["solution"]["z_onoff_var"][14,:]), 1, atol=tol_0)
        @test data["gates_dict"]["14"]["qubit_loc"] == "qubit_3"
        @test isapprox(rad2deg(data["gates_dict"]["14"]["angle"]["θ"]),  45, atol=tol_0)
        @test isapprox(rad2deg(data["gates_dict"]["14"]["angle"]["ϕ"]), -90, atol=tol_0)
        @test isapprox(rad2deg(data["gates_dict"]["14"]["angle"]["λ"]),  90, atol=tol_0)
    end

end

@testset "Tests: 3-qubit CRX gate decomposition" begin

    params = Dict{String, Any}(
        "num_qubits" => 3, 
        "depth" => 2,    
        "elementary_gates" => ["CU3_1_2", "CU3_2_3", "CU3_1_3", "Identity"],  
        "target_gate" => QCO.get_full_sized_gate("CRX_1_3", 3, angle = pi/4),
        "CU3_θ_discretization" => [0, π/4],
        "CU3_ϕ_discretization" => [0, -π/2],
        "CU3_λ_discretization" => [0, π/2],    
        "objective" => "minimize_cnot", 
        "decomposition_type" => "exact"                  
    )

    result_qc = QCO.run_QCModel(params, CBC, model_type = "balas_formulation")
    
    data = QCO.get_data(params)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 1, atol = tol_0)
    if isapprox(sum(result_qc["solution"]["z_onoff_var"][18,:]), 1, atol=tol_0)
        @test data["gates_dict"]["18"]["qubit_loc"] == "qubit_1_3"
        @test isapprox(rad2deg(data["gates_dict"]["18"]["angle"]["θ"]),  45, atol=tol_0)
        @test isapprox(rad2deg(data["gates_dict"]["18"]["angle"]["ϕ"]), -90, atol=tol_0)
        @test isapprox(rad2deg(data["gates_dict"]["18"]["angle"]["λ"]),  90, atol=tol_0)
    end

end

@testset "Tests: feasibility problem" begin
    params = Dict{String, Any}(
        "num_qubits" => 2, 
        "depth" => 3,    
        "elementary_gates" => ["H_1", "H_2"],  
        "target_gate" => QCO.kron_single_qubit_gate(2, QCO.HGate(), "q1"),
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact"
        )

    result_qc = QCO.run_QCModel(params, CBC, all_valid_constraints = -1)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    data = QCO.get_data(params)

    for i in keys(data["gates_dict"])
        if data["gates_dict"][i]["type"] == "H_1"
            @test isapprox(sum(result_qc["solution"]["z_onoff_var"][parse(Int64, i),:]), 3 , atol = tol_0)
        end
    end
  
end

@testset "Tests: JuMP set_start_value for z_onoff_var variables" begin
    function input_circuit()
        # [(depth, gate)]
        return [(1, "CNot_2_1"), 
                (2, "S_1"), 
                (3, "H_2"), 
                (4, "S_2")
                ]
    end

    params = Dict{String, Any}(
    "num_qubits" => 2,
    "depth" => 4,
    "elementary_gates" => ["S_1", "S_2", "H_1", "H_2", "CNot_1_2", "CNot_2_1", "Identity"], 
    "target_gate" => QCO.MGate(),
    "input_circuit" => input_circuit(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact"
    )

    result_qc = QCO.run_QCModel(params, CBC)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][6,:]), 1, atol=tol_0)
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][5,:]), 0, atol=tol_0)
    
end

@testset "Tests: Involutory gate constraints" begin
    
    params = Dict{String, Any}(
    "num_qubits" => 2,
    "depth" => 4,
    "elementary_gates" => ["S_1", "S_2", "H_1", "H_2", "CNot_1_2", "CNot_2_1", "Identity"], 
    "target_gate" => QCO.MGate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact"
    )

    gates_dict = QCO.get_data(params)["gates_dict"]
    involutory_gates = QCO.get_involutory_gates(gates_dict)
    @test length(involutory_gates) == 4 #excluding Identity gate
     
    for i in keys(gates_dict)
        if ("S_1" in gates_dict[i]["type"]) || ("S_2" in gates_dict[i]["type"])
            @test !(parse(Int, i) in involutory_gates)
        end
        if ("H_1" in gates_dict[i]["type"]) || ("H_2" in gates_dict[i]["type"]) || ("CNot_1_2" in gates_dict[i]["type"]) || ("CNot_2_1" in gates_dict[i]["type"]) 
            @test (parse(Int, i) in involutory_gates)
        end
    end

    result_qc = QCO.run_QCModel(params, CBC)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 4, atol = tol_0)
    
end

@testset "Tests: TIME_LIMIT for building results dict and log" begin
    
    params = Dict{String, Any}(
    "num_qubits" => 2,
    "depth" => 5,
    "elementary_gates" => ["S_1", "S_2", "H_1", "H_2", "CNot_1_2", "CNot_2_1", "Identity"], 
    "target_gate" => QCO.MGate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    "time_limit" => 1
    )

    result_qc = QCO.run_QCModel(params, CBC)
    @test result_qc["termination_status"] == MOI.TIME_LIMIT
    @test result_qc["primal_status"] == MOI.NO_SOLUTION
    
end

@testset "Tests: constraint_redundant_gate_product_pairs" begin
    
    function target_gate()
        T1 = QCO.get_full_sized_gate("U3_2", 2, angle = [0,π/2,π])
        T2 = QCO.get_full_sized_gate("U3_1", 2, angle = [π/2,π/2,-π/2])
        return T1*T2
    end
    
    params = Dict{String, Any}(
               "num_qubits" => 2, 
               "depth" => 2,    
               "elementary_gates" => ["U3_1", "U3_2", "Identity"], 
               "target_gate" => target_gate(),   
               "U3_θ_discretization" => [0, π/2],
               "U3_ϕ_discretization" => [π/2],
               "U3_λ_discretization" => [-π/2, π], 
               "objective" => "minimize_depth")

    data = QCO.get_data(params)
    redundant_pairs = QCO.get_redundant_gate_product_pairs(data["gates_dict"])
    @test length(redundant_pairs) == 2
    @test redundant_pairs[1] == (1,4)
    @test redundant_pairs[2] == (5,7)
    @test isapprox(data["gates_dict"]["1"]["matrix"] * data["gates_dict"]["4"]["matrix"], data["gates_dict"]["2"]["matrix"], atol = tol_0)
    @test isapprox(data["gates_dict"]["5"]["matrix"] * data["gates_dict"]["7"]["matrix"], data["gates_dict"]["6"]["matrix"], atol = tol_0)

    result_qc = QCO.run_QCModel(params, CBC)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 2.0, atol = tol_0)

end

@testset "Tests: constraint_idempotent_gates" begin
    params = Dict{String, Any}(
    "num_qubits" => 2, 
    "depth" => 3,    
    "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "Identity"], 
    "target_gate" => QCO.CZGate(),
    "identify_real_gates" => true,
    "U3_θ_discretization" => [-π/2, 0, π/2],
    "U3_ϕ_discretization" => [0, π/2],
    "U3_λ_discretization" => [0, π/2])

    data = QCO.get_data(params)
    idempotent_pairs = QCO.get_idempotent_gates(data["gates_dict"])
    @test length(idempotent_pairs) == 2

    result_qc = QCO.run_QCModel(params, CBC, all_valid_constraints = 1)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 3.0, atol = tol_0)

    for i in keys(data["gates_dict"])
        if "Identity" in data["gates_dict"][i]["type"]
            @test isapprox(sum(result_qc["solution"]["z_onoff_var"][parse(Int64, i),:]), 0, atol = tol_0)
        end
    end

end
