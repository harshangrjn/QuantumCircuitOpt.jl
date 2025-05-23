@testset "QC_model Tests: Minimize depth for controlled-NOT gate" begin

    params = Dict{String, Any}(
    "num_qubits" => 2, 
    "maximum_depth" => 5,
    "elementary_gates" => ["H_1", "H_2", "CNot_1_2", "Identity"],  
    "initial_gate" => "Identity",
    "target_gate" => QCO.CNotRevGate(),
    "set_cnot_lower_bound" => 1,
    "set_cnot_upper_bound" => 1,
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal"
    )

    model_options = Dict{Symbol, Any}(:model_type => "compact_formulation_1", # Testing incorrect model_type
                                      :all_valid_constraints => 2,            # Testing incorrect all_valid_constraints
                                      :unitary_complex_conjugate => 1,
                                     )
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 5, atol=tol_0)
    @test isapprox(result_qc["solution"]["z_bin_var"][1,1], 1, atol=tol_0) || isapprox(result_qc["solution"]["z_bin_var"][2,1], 1, atol=tol_0)
    @test isapprox(result_qc["solution"]["z_bin_var"][1,2], 1, atol=tol_0) || isapprox(result_qc["solution"]["z_bin_var"][2,2], 1, atol=tol_0)
    @test isapprox(result_qc["solution"]["z_bin_var"][3,3], 1, atol=tol_0)
    @test isapprox(result_qc["solution"]["z_bin_var"][1,4], 1, atol=tol_0) || isapprox(result_qc["solution"]["z_bin_var"][2,4], 1, atol=tol_0)
    @test isapprox(result_qc["solution"]["z_bin_var"][1,5], 1, atol=tol_0) || isapprox(result_qc["solution"]["z_bin_var"][2,5], 1, atol=tol_0)
    
end

@testset "QC_model Tests: Minimum depth T gate decomposition" begin
    target_gate = QCO.unitary("T_1", 2) * 
                  QCO.unitary("CNot_1_2", 2) * 
                  QCO.unitary("T_2", 2) * 
                  QCO.unitary("CNot_1_2", 2) * 
                  QCO.unitary("T_1", 2)

    params = Dict{String, Any}(
    "num_qubits" => 2,
    "maximum_depth" => 5,
    "elementary_gates" => ["T_1", "T_2", "Tdagger_1", "Tdagger_2", "CNot_1_2", "CNot_2_1", "Identity"],
    "target_gate" => target_gate,
    "objective" => "minimize_T",
    "decomposition_type" => "optimal_global_phase"
    )

    result_qc = QCO.run_QCModel(params, MIP_SOLVER)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 3, atol=tol_0)

    params["elementary_gates"] = ["T_1", "T_2", "S_1", "S_2", "Tdagger_1", "Tdagger_2", "CNot_1_2", "CNot_2_1", "Identity"]
    params["decomposition_type"] = "exact_optimal"

    result_qc = QCO.run_QCModel(params, MIP_SOLVER)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 1, atol=tol_0)
end


@testset "QC_model Tests: Minimum CNOT swap gate decomposition" begin

    params = Dict{String, Any}(
        "num_qubits" => 2,
        "maximum_depth" => 4,        
        "elementary_gates" => ["CNot_1_2", "CNot_2_1", "Identity"],
        "target_gate" => QCO.SwapGate(),  
        "set_cnot_lower_bound" => 2,
        "set_cnot_upper_bound" => 3,
        "objective" => "minimize_cnot", 
        "decomposition_type" => "exact_optimal"                      
        )
    
    model_options = Dict{Symbol, Any}(:model_type => "balas_formulation",
                                      :all_valid_constraints => 1,
                                      :unitary_complex_conjugate => 1,
                                      )

    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 3, atol = tol_0)
    @test isapprox(sum(result_qc["solution"]["z_bin_var"][1:2,:]), 3, atol=tol_0)
    
end

@testset "QC_model Tests: Minimum depth U3(0,0,π/4) gate" begin

    params = Dict{String, Any}(
    "num_qubits" => 2, 
    "maximum_depth" => 2,    
    "elementary_gates" => ["U3_1", "U3_2", "Identity"],  
    "target_gate" => QCO.kron_single_qubit_gate(2, QCO.U3Gate(0,0,π/4), "q1"),
    "U3_θ_discretization" => [0, π/2],
    "U3_ϕ_discretization" => [0],
    "U3_λ_discretization" => [0, π/4],
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal"                  
    )

    result_qc = QCO.run_QCModel(params, MIP_SOLVER)

    data = QCO.get_data(params)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 1, atol = tol_0)
    @test "Identity" in data["gates_dict"]["3"]["type"]
    @test isapprox(sum(result_qc["solution"]["z_bin_var"][3,:]), 1, atol = tol_0)
    @test isapprox(sum(result_qc["solution"]["z_bin_var"][4,:]), 1, atol = tol_0) 
    @test data["gates_dict"]["4"]["qubit_loc"] == [1]
end

@testset "QC_model Tests: Minimum depth CU3(0,0,π/4) gate" begin

    params = Dict{String, Any}(
    "num_qubits" => 2, 
    "maximum_depth" => 2,    
    "elementary_gates" => ["CU3_1_2", "CU3_2_1", "Identity"],  
    "target_gate" => QCO.CU3Gate(0, 0, π/4),
    "CU3_θ_discretization" => [0, π/2],
    "CU3_ϕ_discretization" => [0],
    "CU3_λ_discretization" => [0, π/4],
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal"                  
    )

    result_qc = QCO.run_QCModel(params, MIP_SOLVER)

    data = QCO.get_data(params)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 1, atol = tol_0)
    @test "Identity" in data["gates_dict"]["3"]["type"]
    @test isapprox(sum(result_qc["solution"]["z_bin_var"][3,:]), 1, atol = tol_0)
    @test isapprox(sum(result_qc["solution"]["z_bin_var"][4,:]), 1, atol = tol_0) 
    @test data["gates_dict"]["4"]["qubit_loc"] == [1,2]
end

@testset "QC_model Tests: Minimum depth Rev CU3(0,0,π/4) gate" begin

    params = Dict{String, Any}(
    "num_qubits" => 3, 
    "maximum_depth" => 2,    
    "elementary_gates" => ["CU3_3_1", "CU3_1_3", "Identity"],  
    "target_gate" => QCO.unitary("CU3_3_1", 3, angle = [0, 0, pi/4]),
    "CU3_θ_discretization" => [0, π/2],
    "CU3_ϕ_discretization" => [0],
    "CU3_λ_discretization" => [0, π/4],
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal"                  
    )

    result_qc = QCO.run_QCModel(params, MIP_SOLVER)

    data = QCO.get_data(params)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 1, atol = tol_0)
    @test "Identity" in data["gates_dict"]["3"]["type"]
    @test isapprox(sum(result_qc["solution"]["z_bin_var"][3,:]), 1, atol = tol_0)
    @test isapprox(sum(result_qc["solution"]["z_bin_var"][4,:]), 1, atol = tol_0) 
    @test data["gates_dict"]["4"]["qubit_loc"] == [3,1]
end

@testset "QC_model Tests: Minimum depth RX, RY, RZ gate decomposition" begin

    params = Dict{String, Any}(
    "num_qubits" => 2, 
    "maximum_depth" => 3,
    "elementary_gates" => ["RX_1", "RY_2", "RZ_1", "Identity"],  
    "target_gate" => QCO.kron_single_qubit_gate(2, QCO.RXGate(π/4), "q1") * QCO.kron_single_qubit_gate(2, QCO.RYGate(π/4), "q2") * QCO.kron_single_qubit_gate(2, QCO.RZGate(π/4), "q1"),
    "RX_discretization" => [0, π/4],
    "RY_discretization" => [π/4],
    "RZ_discretization" => [π/2, π/4],
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal"                            
    )

    model_options = Dict{Symbol, Any}(:model_type => "balas_formulation",
                                      :commute_gate_constraints => true,
                                      :unitary_complex_conjugate => 1,)

    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options =  model_options)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 3, atol = tol_0)

end

@testset "QC_model Tests: Minimum depth CRX, CRY, CRZ gate decomposition" begin

    params = Dict{String, Any}(
    "num_qubits" => 3, 
    "maximum_depth" => 3,
    "elementary_gates" => ["CRX_1_2", "CRY_2_3", "CRZ_3_1", "Identity"],  
    "target_gate" => QCO.unitary("CRX_1_2", 3, angle = pi/4) * QCO.unitary("CRY_2_3", 3, angle=pi/4) * QCO.unitary("CRZ_3_1", 3, angle=pi/4),
    "CRX_discretization" => [0, π/4],
    "CRY_discretization" => [π/4],
    "CRZ_discretization" => [π/2, π/4],
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal"                            
    )

    model_options = Dict{Symbol, Any}(:model_type => "balas_formulation",
                                      :commute_gate_constraints => true,
                                      :unitary_complex_conjugate => 1,)

    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options =  model_options)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 3, atol = tol_0)

end

@testset "QC_model Tests: Minimum depth Rev CRX, CRY, CRZ gate decomposition" begin

    params = Dict{String, Any}(
    "num_qubits" => 3, 
    "maximum_depth" => 3,
    "elementary_gates" => ["CRX_3_1", "CRY_3_1", "CRZ_1_3", "Identity"],  
    "target_gate" => QCO.unitary("CRX_3_1", 3, angle=pi/4) * QCO.unitary("CRY_3_1", 3, angle=pi/4) * QCO.unitary("CRZ_1_3", 3, angle = pi/2),
    "CRX_discretization" => [0, π/4],
    "CRY_discretization" => [π/4],
    "CRZ_discretization" => [π/2, π/4],
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal"                            
    )

    model_options = Dict{Symbol, Any}(:model_type => "balas_formulation",
                                      :commute_gate_constraints => true,
                                      :unitary_complex_conjugate => 1,)

    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options =  model_options)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 3, atol = tol_0)

end


@testset "QC_model Tests: 3-qubit RX gate decomposition" begin

    params = Dict{String, Any}(
    "num_qubits" => 3, 
    "maximum_depth" => 2,    
    "elementary_gates" => ["U3_1", "U3_2", "U3_3", "Identity"],  
    "target_gate" => QCO.kron_single_qubit_gate(3, QCO.RXGate(π/4), "q3"),
    "U3_θ_discretization" => [0, π/4],
    "U3_ϕ_discretization" => [0, -π/2],
    "U3_λ_discretization" => [0, π/2],    
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal"                   
    )

    model_options = Dict{Symbol, Any}(:model_type => "balas_formulation",
                                      :unitary_complex_conjugate => 0,)

    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options =  model_options)
    
    data = QCO.get_data(params)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 1, atol = tol_0)
    if isapprox(sum(result_qc["solution"]["z_bin_var"][14,:]), 1, atol=tol_0)
        @test data["gates_dict"]["14"]["qubit_loc"] == [3]
        @test isapprox(rad2deg(data["gates_dict"]["14"]["angle"]["θ"]),  45, atol=tol_0)
        @test isapprox(rad2deg(data["gates_dict"]["14"]["angle"]["ϕ"]), -90, atol=tol_0)
        @test isapprox(rad2deg(data["gates_dict"]["14"]["angle"]["λ"]),  90, atol=tol_0)
    end

end

@testset "QC_model Tests: 3-qubit CRX gate decomposition" begin

    params = Dict{String, Any}(
        "num_qubits" => 3, 
        "maximum_depth" => 2,    
        "elementary_gates" => ["CU3_1_2", "CU3_2_3", "CU3_1_3", "Identity"],  
        "target_gate" => QCO.unitary("CRX_1_3", 3, angle = pi/4),
        "CU3_θ_discretization" => [0, π/4],
        "CU3_ϕ_discretization" => [0, -π/2],
        "CU3_λ_discretization" => [0, π/2],    
        "objective" => "minimize_depth",
        "decomposition_type" => "exact_optimal"                  
    )

    model_options = Dict{Symbol, Any}(:model_type => "balas_formulation",
                                      :unitary_complex_conjugate => 1,)

    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options =  model_options)
    
    data = QCO.get_data(params)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 1, atol = tol_0)
    if isapprox(sum(result_qc["solution"]["z_bin_var"][18,:]), 1, atol=tol_0)
        @test data["gates_dict"]["18"]["qubit_loc"] == [1,3]
        @test isapprox(rad2deg(data["gates_dict"]["18"]["angle"]["θ"]),  45, atol=tol_0)
        @test isapprox(rad2deg(data["gates_dict"]["18"]["angle"]["ϕ"]), -90, atol=tol_0)
        @test isapprox(rad2deg(data["gates_dict"]["18"]["angle"]["λ"]),  90, atol=tol_0)
    end

end

@testset "QC_model Tests: feasibility problem" begin
    params = Dict{String, Any}(
        "num_qubits" => 2, 
        "maximum_depth" => 3,    
        "elementary_gates" => ["H_1", "H_2"],  
        "target_gate" => QCO.kron_single_qubit_gate(2, QCO.HGate(), "q1"),
        "objective" => "minimize_depth", 
        "decomposition_type" => "exact_optimal"
        )
    
    model_options = Dict{Symbol, Any}(:all_valid_constraints => -1)

    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options =  model_options)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    data = QCO.get_data(params)

    for i in keys(data["gates_dict"])
        if data["gates_dict"][i]["type"] == "H_1"
            @test isapprox(sum(result_qc["solution"]["z_bin_var"][parse(Int64, i),:]), 3 , atol = tol_0)
        end
    end
  
end

@testset "QC_model Tests: JuMP set_start_value for z_bin_var variables" begin
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
    "maximum_depth" => 4,
    "elementary_gates" => ["S_1", "S_2", "H_1", "H_2", "CNot_1_2", "CNot_2_1", "Identity"], 
    "target_gate" => QCO.MGate(),
    "input_circuit" => input_circuit(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal"
    )
    model_options = Dict{Symbol, Any}(:unitary_complex_conjugate => 1,)
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(sum(result_qc["solution"]["z_bin_var"][6,:]), 1, atol=tol_0)
    @test isapprox(sum(result_qc["solution"]["z_bin_var"][5,:]), 0, atol=tol_0)
    
end

@testset "QC_model Tests: Involutory gate constraints" begin
    
    params = Dict{String, Any}(
    "num_qubits" => 2,
    "maximum_depth" => 4,
    "elementary_gates" => ["S_1", "S_2", "H_1", "H_2", "CNot_1_2", "CNot_2_1", "Identity"], 
    "target_gate" => QCO.MGate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal"
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
    model_options = Dict{Symbol, Any}(:unitary_complex_conjugate => 1,)
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 4, atol = tol_0)
    
end

@testset "QC_model Tests: TIME_LIMIT for building results dict and log" begin
    
    params = Dict{String, Any}(
    "num_qubits" => 2,
    "maximum_depth" => 5,
    "elementary_gates" => ["S_1", "S_2", "H_1", "H_2", "CNot_1_2", "CNot_2_1", "Identity"], 
    "target_gate" => QCO.MGate(),
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact_optimal"
    )

    model_options = Dict{Symbol, Any}(:time_limit => 0.1,
                                      :unitary_complex_conjugate => 0,)

    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)
    @test result_qc["termination_status"] == MOI.TIME_LIMIT
    @test result_qc["primal_status"]      == MOI.NO_SOLUTION
    
end

@testset "QC_model Tests: constraint_redundant_gate_product_pairs" begin
    
    function target_gate()
        T1 = QCO.unitary("U3_2", 2, angle = [0,π/2,π])
        T2 = QCO.unitary("U3_1", 2, angle = [π/2,π/2,-π/2])
        return T1*T2
    end
    
    params = Dict{String, Any}(
               "num_qubits" => 2, 
               "maximum_depth" => 2,    
               "elementary_gates" => ["U3_1", "U3_2", "Identity"], 
               "target_gate" => target_gate(),   
               "U3_θ_discretization" => [0, π/2],
               "U3_ϕ_discretization" => [π/2],
               "U3_λ_discretization" => [-π/2, π], 
               "objective" => "minimize_depth")

    data = QCO.get_data(params)
    redundant_pairs = QCO.get_redundant_gate_product_pairs(data["gates_dict"], data["decomposition_type"])
    @test length(redundant_pairs) == 2
    @test redundant_pairs[1] == (1,4)
    @test redundant_pairs[2] == (5,7)
    @test isapprox(data["gates_dict"]["1"]["matrix"] * data["gates_dict"]["4"]["matrix"], data["gates_dict"]["2"]["matrix"], atol = tol_0)
    @test isapprox(data["gates_dict"]["5"]["matrix"] * data["gates_dict"]["7"]["matrix"], data["gates_dict"]["6"]["matrix"], atol = tol_0)

    result_qc = QCO.run_QCModel(params, MIP_SOLVER)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 2.0, atol = tol_0)

end

@testset "QC_model Tests: constraint_idempotent_gates" begin
    params = Dict{String, Any}(
    "num_qubits" => 2, 
    "maximum_depth" => 3,    
    "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "Identity"], 
    "target_gate" => QCO.CZGate(),
    "U3_θ_discretization" => [-π/2, 0, π/2],
    "U3_ϕ_discretization" => [0, π/2],
    "U3_λ_discretization" => [0, π/2])

    data = QCO.get_data(params)
    idempotent_pairs = QCO.get_idempotent_gates(data["gates_dict"], data["decomposition_type"])
    @test length(idempotent_pairs) == 2

    model_options = Dict{Symbol, Any}(:idempotent_gate_constraints => true,
                                      :unitary_complex_conjugate => 1,)
    
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 3.0, atol = tol_0)

    for i in keys(data["gates_dict"])
        if "Identity" in data["gates_dict"][i]["type"]
            @test isapprox(sum(result_qc["solution"]["z_bin_var"][parse(Int64, i),:]), 0, atol = tol_0)
        end
    end

end

@testset "QC_model Tests: constraint_convex_hull_complex_gates" begin
    
    function target_gate()
        num_qubits = 4
        CV_1_4 = QCO.unitary("CV_1_4", num_qubits);
        CV_2_4 = QCO.unitary("CV_2_4", num_qubits);
        CV_3_4 = QCO.unitary("CV_3_4", num_qubits);
        CVdagger_3_4 = QCO.unitary("CVdagger_3_4", num_qubits);        
        CNot_1_2 = QCO.unitary("CNot_1_2", num_qubits);
        CNot_2_3 = QCO.unitary("CNot_2_3", num_qubits);

        return CV_2_4 * CNot_1_2 * CV_3_4 * CV_1_4 * CNot_2_3 * CVdagger_3_4
    end

    params = Dict{String, Any}(
    "num_qubits" => 4,
    "maximum_depth" => 6,
    "elementary_gates" => ["CV_1_4", "CV_2_4", "CV_3_4", "CVdagger_3_4", "CNot_1_2", "CNot_2_3", "Identity"],
    "target_gate" => target_gate()
    )
    
    model_options = Dict{Symbol, Any}(:relax_integrality => true,
                                      :convex_hull_gate_constraints => true,
                                      :fix_unitary_variables => true,
                                      :optimizer_log => false,
                                      :unitary_complex_conjugate => 0,
                                      )

    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 0.6, atol = tol_0)
end

@testset "QC_model Tests: RGate decomposition" begin
    
    num_qubits = 2
    GR1      = QCO.GRGate(num_qubits, π/6, π/3)
    CNot_1_2 = QCO.unitary("CNot_1_2", 2)
    T        = QCO.round_complex_values(GR1 * CNot_1_2)

    params = Dict{String, Any}(
    "num_qubits" => num_qubits,
    "maximum_depth" => 3,
    "elementary_gates" => ["R_1", "R_2", "CNot_1_2", "Identity"],
    "R_θ_discretization" => [π/3, π/6],
    "R_ϕ_discretization" => [π/3, π/6],
    "target_gate" => T,
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal",
    )
    
    model_options = Dict{Symbol, Any}(:optimizer_log => false, 
                                      :fix_unitary_variables => false,
                                      :unitary_complex_conjugate => 1,)
                                      
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 3.0, atol = tol_0)
    z_sol = result_qc["solution"]["z_bin_var"]
    @test isapprox(sum(z_sol[1:8, :]), 2, atol = tol_0)
    @test isapprox(sum(z_sol[9, :]),   1, atol = tol_0)
end

@testset "QC_model Tests: GRGate decomposition for Pauli-X gate, and unitary constraints" begin

    params = Dict{String, Any}(
    "num_qubits" => 2,
    "maximum_depth" => 5, 
    "elementary_gates" => ["GR", "RZ_2", "CZ_1_2", "Identity"], 
    "target_gate" => QCO.unitary("X_1", 2),
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal",
    "GR_θ_discretization" => [-π/2, π/2],
    "GR_ϕ_discretization" => [0],
    "RZ_discretization"   => [π],
    )
    
    model_options = Dict{Symbol, Any}(:optimizer_log => false, 
                                      :unitary_constraints => true,
                                      :unitary_complex_conjugate => 1,)
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 5.0, atol = tol_0)
    z_sol = result_qc["solution"]["z_bin_var"]
    @test isapprox(sum(z_sol[4, :]), 2, atol = tol_0) # test for number of CZ gates 
end

@testset "QC_model Tests: Decomposition using exact_feasible option" begin
    params = Dict{String, Any}(
        "num_qubits" => 2,
        "maximum_depth" => 5,
        "elementary_gates" => ["CH_1_2", "CNot_2_1", "CNot_1_2", "Identity"],
        "target_gate" => QCO.WGate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "exact_feasible",
        )

    model_options = Dict{Symbol, Any}(:optimizer_log => false, :unitary_complex_conjugate => 1,)
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 0.0, atol = tol_0)
    z_sol = result_qc["solution"]["z_bin_var"]
    @test isapprox(sum(z_sol[2:4, :]), 4, atol = tol_0)
end

@testset "QC_model Tests: Approximate decomposition using outer approximation-1" begin
    params = Dict{String, Any}(
        "num_qubits" => 2,
        "maximum_depth" => 5,
        "elementary_gates" => ["H_1", "H_2", "CNot_1_2", "Identity"],
        "target_gate" => QCO.CNotRevGate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "approximate",
        )
    model_options = Dict{Symbol, Any}(:optimizer_log => false, :unitary_complex_conjugate => 1,)
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 5.0, atol = tol_0)
end

@testset "QC_model Tests: Approximate decomposition using outer approximation-2" begin
    params = Dict{String, Any}(
        "num_qubits" => 2,
        "maximum_depth" => 5,
        "elementary_gates" => ["H_1", "H_2", "CNot_1_2", "Identity"],
        "target_gate" => QCO.CNotRevGate(),
        "objective" => "minimize_depth",
        "decomposition_type" => "approximate",
        )
    model_options = Dict{Symbol, Any}(:optimizer_log => false, :unitary_complex_conjugate => 1,)
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 5.0, atol = tol_0)

    # Testing approximate decomposition for balas_formulation
    model_options = Dict{Symbol, Any}(:optimizer_log => false, :model_type => "balas_formulation", :unitary_complex_conjugate => 1,)
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 5.0, atol = tol_0)

    # Testing approximate decomposition for feasibility case
    params["elementary_gates"] = ["H_1", "H_2", "CNot_1_2"]
    model_options = Dict{Symbol, Any}(:optimizer_log => false, :unitary_complex_conjugate => 1,)
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 0.0, atol = tol_0)

    # Testing approximate decomposition for minimizing CNOT gates
    params["elementary_gates"] = ["H_1", "H_2", "CNot_1_2", "Identity"]
    params["objective"] = "minimize_cnot"
    model_options = Dict{Symbol, Any}(:optimizer_log => false, :unitary_complex_conjugate => 1,)
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 1.0, atol = tol_0)

    # Testing approximate decomposition for case when slack_var-s are fixed based on U_var-s
    params = Dict{String, Any}(
    "num_qubits" => 3,
    "maximum_depth" => 7,
    "elementary_gates" => ["CV_1_2", "CV_2_3", "CV_1_3", "CVdagger_1_2", "CVdagger_2_3", "CVdagger_1_3", "CNot_1_2", "CNot_3_2", "CNot_2_3", "CNot_1_3", "Identity"],
    "target_gate" => QCO.CSwapGate(), #also Fredkin
    "objective" => "minimize_depth",
    "decomposition_type" => "approximate"
    )
    model_options = Dict{Symbol, Any}(:optimizer_log => false, :relax_integrality => true, :fix_unitary_variables => true, :unitary_complex_conjugate => 0,)
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 0.38281250, atol = tol_0)

end

@testset "QC_model Tests: Global phase constraints" begin

    #Real elementary gates - Target is real but with a minus sign
    params = Dict{String, Any}(
        "num_qubits" => 2,
        "maximum_depth" => 5,
        "elementary_gates" => ["H_1", "H_2", "CNot_1_2", "Identity"],
        "target_gate" => - QCO.CNotRevGate(),
        "objective" => "minimize_depth",
        )
    model_options = Dict{Symbol, Any}(:optimizer_log => false, :unitary_complex_conjugate => 1,)

    # Without global phase constraints
    params["decomposition_type"] = "exact_optimal"
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)
    @test result_qc["termination_status"] == MOI.INFEASIBLE
    # @test result_qc["primal_status"]      == MOI.NO_SOLUTION
    
    # With global phase constraints
    params["decomposition_type"] = "optimal_global_phase"
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 5.0, atol = tol_0)

    # Real elementary gates - Target is complex, but real in a global phase sence
    params["target_gate"] = exp(im*pi*0.3) * QCO.CNotRevGate()

    # Without global phase constraints 
    # This is an infeasible decomposition: all elementary gates have zero imaginary parts and 
    # target is not real for exact decomposition or not real up to a global phase for 
    # `optimal_global_phase` decomposition.
    
    # With global phase constraints
    params["decomposition_type"] = "optimal_global_phase"
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 5.0, atol = tol_0)


    # Using complex-valued elementary gates
    params["elementary_gates"] = ["T_1", "H_1", "H_2", "CNot_1_2", "Identity"]

    # Without global phase constraints
    params["decomposition_type"] = "exact_optimal"
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)
    @test result_qc["termination_status"] == MOI.INFEASIBLE
    # @test result_qc["primal_status"]      == MOI.NO_SOLUTION

    # Using global phase constraints
    params["decomposition_type"] = "optimal_global_phase"
    result_qc = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    @test isapprox(result_qc["objective"], 5.0, atol = tol_0)

end