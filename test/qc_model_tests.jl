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
    @test result_qc["objective"] == 0
    @test (result_qc["solution"]["z_onoff_var"][1,1] == 1) || (result_qc["solution"]["z_onoff_var"][2,1] == 1)
    @test (result_qc["solution"]["z_onoff_var"][1,2] == 1) || (result_qc["solution"]["z_onoff_var"][2,2] == 1)
    @test result_qc["solution"]["z_onoff_var"][3,3] == 1
    @test (result_qc["solution"]["z_onoff_var"][1,4] == 1) || (result_qc["solution"]["z_onoff_var"][2,4] == 1)
    @test (result_qc["solution"]["z_onoff_var"][1,5] == 1) || (result_qc["solution"]["z_onoff_var"][2,5] == 1)
    
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
    @test isapprox(result_qc["objective"], 3, atol = 1E-6)
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][1:2,:]), 3, atol=1E-6)
    
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
    @test isapprox(result_qc["objective"], 2, atol = 1E-6)
    @test "Identity" in data["gates_dict"]["5"]["type"]
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][5,:]), 2, atol = 1E-6)
    @test isapprox(sum(result_qc["solution"]["z_onoff_var"][7,:]), 1, atol = 1E-6) 
    @test data["gates_dict"]["7"]["qubit_location"] == "qubit_1"
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

    result_qc = QCO.run_QCModel(params, CBC, model_type = "balas_formulation", commute_matrix_cuts = true)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 0, atol = 1E-6)

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
    @test isapprox(result_qc["objective"], 1, atol = 1E-6)
    if isapprox(sum(result_qc["solution"]["z_onoff_var"][14,:]), 1, atol=1E-6)
        @test data["gates_dict"]["14"]["qubit_location"] == "qubit_3"
        @test isapprox(rad2deg(data["gates_dict"]["14"]["angle"]["θ"]),  45, atol=1E-6)
        @test isapprox(rad2deg(data["gates_dict"]["14"]["angle"]["ϕ"]), -90, atol=1E-6)
        @test isapprox(rad2deg(data["gates_dict"]["14"]["angle"]["λ"]),  90, atol=1E-6)
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
            @test isapprox(sum(result_qc["solution"]["z_onoff_var"][parse(Int64, i),:]), 3 , atol = 1E-6)
        end
    end
  
end