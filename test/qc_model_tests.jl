@testset "Test: Minimum depth controlled-NOT gate" begin

    params = Dict{String, Any}(
    "num_qubits" => 2, 
    "depth" => 5,    

    "elementary_gates" => ["H1", "H2", "cnot_12", "Identity"],  
    "target_gate" => QCO.CNotRevGate(),

    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cbc"
    )

    result_qc = QCO.run_QCModel(params, CBC, model_type = "compact_formulation", visualize_solution=true)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 5, atol=tol_0)
    @test isapprox(result_qc["solution"]["z_onoff_var"][1,1], 1, atol=tol_0) || isapprox(result_qc["solution"]["z_onoff_var"][2,1], 1, atol=tol_0)
    @test isapprox(result_qc["solution"]["z_onoff_var"][1,2], 1, atol=tol_0) || isapprox(result_qc["solution"]["z_onoff_var"][2,2], 1, atol=tol_0)
    @test isapprox(result_qc["solution"]["z_onoff_var"][3,3], 1, atol=tol_0)
    @test isapprox(result_qc["solution"]["z_onoff_var"][1,4], 1, atol=tol_0) || isapprox(result_qc["solution"]["z_onoff_var"][2,4], 1, atol=tol_0)
    @test isapprox(result_qc["solution"]["z_onoff_var"][1,5], 1, atol=tol_0) || isapprox(result_qc["solution"]["z_onoff_var"][2,5], 1, atol=tol_0)
    
end

@testset "Test: Minimum CNOT swap gate decomposition" begin

    params = Dict{String, Any}(
        "num_qubits" => 2,
        "depth" => 4,
        
        "elementary_gates" => ["cnot_12", "cnot_21", "Identity"],
        "target_gate" => QCO.SwapGate(),
    
        "objective" => "minimize_cnot", 
        "decomposition_type" => "exact",
        "optimizer" => "cbc"                         
        )

    result_qc = QCO.run_QCModel(params, CBC, model_type = "balas_formulation", visualize_solution=true)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 3, atol = tol_0)
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][1:2,:]), 3, atol=tol_0)
    
end

@testset "Test: Minimum depth U3(0,0,π/4) gate" begin

    params = Dict{String, Any}(
    "num_qubits" => 2, 
    "depth" => 3,    

    "elementary_gates" => ["U3", "Identity"],  
    "target_gate" => QCO.kron_single_gate(2, QCO.U3Gate(0,0,π/4), "q1"),

    "U_θ_discretization" => [0, π/2],
    "U_ϕ_discretization" => [0],
    "U_λ_discretization" => [0, π/4],

    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cbc"                              
    )

    result_qc = QCO.run_QCModel(params, CBC, model_type = "compact_formulation", visualize_solution=true, eliminate_identical_gates = true)

    data = QCO.get_data(params, eliminate_identical_gates = true)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 1, atol = tol_0)
    @test "Identity" in data["gates_dict"]["5"]["type"]
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][5,:]), 2, atol = tol_0)
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][7,:]), 1, atol = tol_0) 
    @test data["gates_dict"]["7"]["qubit_loc"] == "qubit_1"
end

@testset "Test: Minimum depth RX, RY, RZ gate decomposition" begin

    params = Dict{String, Any}(
    "num_qubits" => 2, 
    "depth" => 3,    

    "elementary_gates" => ["RX", "RY", "RZ", "Identity"],  
    "target_gate" => QCO.kron_single_gate(2, QCO.RXGate(π/4), "q1") * QCO.kron_single_gate(2, QCO.RYGate(π/4), "q2") * QCO.kron_single_gate(2, QCO.RZGate(π/4), "q1"),

    "RX_discretization" => [0, π/4],
    "RY_discretization" => [π/4],
    "RZ_discretization" => [π/2, π/4],

    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cbc"                             
    )

    result_qc = QCO.run_QCModel(params, CBC, model_type = "balas_formulation", commute_matrix_constraints = true)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 3, atol = tol_0)

end

@testset "Test: 3-qubit RX gate decomposition" begin

    params = Dict{String, Any}(
    "num_qubits" => 3, 
    "depth" => 2,    

    "elementary_gates" => ["U3", "Identity"],  
    "target_gate" => QCO.kron_single_gate(3, QCO.RXGate(π/4), "q3"),

    "U_θ_discretization" => [0, π/4],
    "U_ϕ_discretization" => [0, -π/2],
    "U_λ_discretization" => [0, π/2],    

    "objective" => "minimize_cnot", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cbc"                               
    )

    result_qc = QCO.run_QCModel(params, CBC, model_type = "balas_formulation", eliminate_identical_gates = true)
    
    data = QCO.get_data(params, eliminate_identical_gates = true)

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

@testset "Feasibility objective tests" begin
    params = Dict{String, Any}(
        "num_qubits" => 2, 
        "depth" => 3,    
    
        "elementary_gates" => ["H1", "H2"],  
        "target_gate" => QCO.kron_single_gate(2, QCO.HGate(), "q1"),
    
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact",
        
        "optimizer" => "cbc"
        )

    result_qc = QCO.run_QCModel(params, CBC)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    data = QCO.get_data(params)

    for i in keys(data["gates_dict"])
        if data["gates_dict"][i]["type"] == "H1"
            @test isapprox(sum(result_qc["solution"]["z_onoff_var"][parse(Int64, i),:]), 3 , atol = tol_0)
        end
    end
  
end

@testset "JuMP set_start_value tests for on-off vars" begin
    function input_circuit()
        # [(depth, gate)]
        return [(1, "cnot_21"), 
                (2, "S1"), 
                (3, "H2"), 
                (4, "S2")
                ]
    end

    params = Dict{String, Any}(
    "num_qubits" => 2,
    "depth" => 4,
    "elementary_gates" => ["S1", "S2", "H1", "H2", "cnot_12", "cnot_21", "Identity"], 
    "target_gate" => QCO.MGate(),
    "input_circuit" => input_circuit(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    "optimizer" => "cbc"
    )

    result_qc = QCO.run_QCModel(params, CBC)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][6,:]), 1, atol=tol_0)
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][5,:]), 0, atol=tol_0)
    
end

@testset "Involutory matrix constraints tests" begin
    
    params = Dict{String, Any}(
    "num_qubits" => 2,
    "depth" => 4,
    "elementary_gates" => ["S1", "S2", "H1", "H2", "cnot_12", "cnot_21", "Identity"], 
    "target_gate" => QCO.MGate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    "optimizer" => "cbc"
    )

    gates_dict = QCO.get_data(params)["gates_dict"]
    
    num_involutory_matrices = 0 
    for i in keys(gates_dict)
        if gates_dict[i]["isInvolutory"]
            num_involutory_matrices += 1
        end

        if ("S1" in gates_dict[i]["type"]) || ("S2" in gates_dict[i]["type"])
            @test !(gates_dict[i]["isInvolutory"])
        end
    end
    @test num_involutory_matrices == 5

    result_qc = QCO.run_QCModel(params, CBC)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 4, atol = tol_0)
    
end

@testset "TIME_LIMIT test for building results dict and log" begin
    
    params = Dict{String, Any}(
    "num_qubits" => 2,
    "depth" => 5,
    "elementary_gates" => ["S1", "S2", "H1", "H2", "cnot_12", "cnot_21", "Identity"], 
    "target_gate" => QCO.MGate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    "optimizer" => "cbc",
    "time_limit" => 1
    )

    result_qc = QCO.run_QCModel(params, CBC)
    @test result_qc["termination_status"] == MOI.TIME_LIMIT
    @test result_qc["primal_status"] == MOI.NO_SOLUTION
    
end