@testset "Depth-aware Tests: minimize_layers on controlled-NOT gate" begin

    params = Dict{String, Any}(
    "num_qubits" => 2,
    "maximum_depth" => 5,
    "elementary_gates" => ["H_1", "H_2", "CNot_1_2", "Identity"],
    "target_gate" => QCO.CNotRevGate(),
    "objective" => "minimize_layers",
    "decomposition_type" => "exact_optimal"
    )

    result_qc = QCO.run_QCModel(params, MIP_SOLVER)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"]      == MOI.FEASIBLE_POINT
    # H on both qubits parallelizes: 5 gates schedule into 3 layers
    @test isapprox(result_qc["objective"], 3, atol = tol_0)
end

@testset "Depth-aware Tests: maximum_layers caps the schedule" begin

    params = Dict{String, Any}(
    "num_qubits" => 2,
    "maximum_depth" => 5,
    "maximum_layers" => 3,
    "elementary_gates" => ["H_1", "H_2", "CNot_1_2", "Identity"],
    "target_gate" => QCO.CNotRevGate(),
    "objective" => "minimize_layers",
    "decomposition_type" => "exact_optimal"
    )

    result_qc = QCO.run_QCModel(params, MIP_SOLVER)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test isapprox(result_qc["objective"], 3, atol = tol_0)
end

@testset "Depth-aware Tests: norm-based approximate types build and solve" begin

    params = Dict{String, Any}(
    "num_qubits" => 2,
    "maximum_depth" => 5,
    "elementary_gates" => ["H_1", "H_2", "CNot_1_2", "Identity"],
    "target_gate" => QCO.CNotRevGate(),
    "objective" => "minimize_depth",
    "decomposition_type" => "sum-approximate"
    )

    result_qc = QCO.run_QCModel(params, MIP_SOLVER)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    # the target is reachable within the depth budget, so the deviation closes to zero
    @test isapprox(result_qc["objective"], 0, atol = 1E-4)

    params["decomposition_type"] = "infinity-approximate"
    result_qc = QCO.run_QCModel(params, MIP_SOLVER)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test isapprox(result_qc["objective"], 0, atol = 1E-4)
end

@testset "Depth-aware Tests: fidelity types build the expected model" begin

    params = Dict{String, Any}(
    "num_qubits" => 2,
    "maximum_depth" => 5,
    "elementary_gates" => ["H_1", "H_2", "CNot_1_2", "Identity"],
    "target_gate" => QCO.CNotRevGate(),
    "objective" => "minimize_depth",
    "decomposition_type" => "linear-fidelity-approximate"
    )

    # the linearized surrogate is a MILP; maximizing the real trace overlap of a
    # reachable target closes to 1
    result_qc = QCO.run_QCModel(params, MIP_SOLVER)
    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test isapprox(result_qc["objective"], 1, atol = 1E-4)

    # the exact fidelity objective is quadratic (nonconvex MIQP); build only
    params["decomposition_type"] = "fidelity-approximate"
    data = QCO.get_data(params)
    qcm = QCO.build_QCModel(data)
    @test haskey(qcm.variables, :w_fidelity_var)
end

@testset "Depth-aware Tests: symmetry-breaking options preserve the optimum" begin

    params = Dict{String, Any}(
    "num_qubits" => 2,
    "maximum_depth" => 5,
    "elementary_gates" => ["H_1", "H_2", "CNot_1_2", "Identity"],
    "target_gate" => QCO.CNotRevGate(),
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal"
    )

    base = QCO.run_QCModel(params, MIP_SOLVER)

    model_options = Dict{Symbol, Any}(
        :redundant_gate_product_constraints => true,
        :equivalent_gate_pair_constraints   => true,
    )
    cut = QCO.run_QCModel(params, MIP_SOLVER; options = model_options)

    @test base["termination_status"] == cut["termination_status"] == MOI.OPTIMAL
    @test isapprox(base["objective"], cut["objective"], atol = tol_0)
end

@testset "Depth-aware Tests: linear_cs global-phase formulation matches the default" begin

    params = Dict{String, Any}(
    "num_qubits" => 2,
    "maximum_depth" => 4,
    "elementary_gates" => ["T_1", "T_2", "Tdagger_1", "Tdagger_2", "CNot_1_2", "Identity"],
    "target_gate" => QCO.unitary("T_1", 2) * QCO.unitary("CNot_1_2", 2) * QCO.unitary("T_2", 2),
    "objective" => "minimize_depth",
    "decomposition_type" => "optimal_global_phase"
    )

    base = QCO.run_QCModel(params, MIP_SOLVER)
    cs = QCO.run_QCModel(params, MIP_SOLVER; options = Dict{Symbol, Any}(:global_phase_formulation => "linear_cs"))

    @test base["termination_status"] == cs["termination_status"] == MOI.OPTIMAL
    @test isapprox(base["objective"], cs["objective"], atol = tol_0)
end
