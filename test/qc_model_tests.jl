@testset "building controlled NOT gate tests" begin
    
    include("../examples/solver.jl")

    params = Dict{String, Any}(
    "n_qubits" => 2, 
    "depth" => 5,    

    "elementary_gates" => ["H1", "H2", "Identity", "cnot_12"],  
    "target_gate" => "cnot_21",

    "initial_gate" => "Identity",
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cbc",
    "presolve" => true,
    "optimizer_log" => false, 
    "relax_integrality" => false,                                
    )

    qcm_optimizer = get_solver(params)
    data = QuantumCircuitOpt.get_data(params)

    model_qc  = QuantumCircuitOpt.build_QCModel(data, model_type = "compact_formulation")
    result_qcm = QuantumCircuitOpt.optimize_QCModel!(model_qc, optimizer = qcm_optimizer)

    @test result_qcm["termination_status"] == MOI.OPTIMAL
    @test result_qcm["primal_status"] == MOI.FEASIBLE_POINT
    @test result_qcm["objective"] == 0
    @test (result_qcm["solution"]["z_onoff_var"][1,1] == 1) || (result_qcm["solution"]["z_onoff_var"][1,2] == 1)
    @test (result_qcm["solution"]["z_onoff_var"][2,1] == 1) || (result_qcm["solution"]["z_onoff_var"][2,2] == 1)
    @test result_qcm["solution"]["z_onoff_var"][4,3] == 1
    @test (result_qcm["solution"]["z_onoff_var"][2,4] == 1) || (result_qcm["solution"]["z_onoff_var"][2,5] == 1)
    @test (result_qcm["solution"]["z_onoff_var"][1,5] == 1) || (result_qcm["solution"]["z_onoff_var"][1,4] == 1)
    
end