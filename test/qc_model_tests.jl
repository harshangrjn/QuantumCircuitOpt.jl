@testset "decomposing controlled-NOT gate tests" begin
    
    # include("../examples/solver.jl")

    params = Dict{String, Any}(
    "num_qubits" => 2, 
    "depth" => 5,    

    "elementary_gates" => ["H1", "H2", "Identity", "cnot_12"],  
    "target_gate" => QCO.CNotRevGate(),

    "initial_gate" => "Identity",
    "objective" => "minimize_depth", 
    "decomposition_type" => "exact",
    
    "optimizer" => "cbc",
    "optimizer_presolve" => false                               
    )

    # qcm_optimizer = get_solver(params)
    result_qc = QCO.run_QCModel(params, CBC, model_type = "compact_formulation", visualize_solution=false)

    @test result_qc["termination_status"] == MOI.OPTIMAL
    @test result_qc["primal_status"] == MOI.FEASIBLE_POINT
    @test result_qc["objective"] == 0
    @test (result_qc["solution"]["z_onoff_var"][1,1] == 1) || (result_qc["solution"]["z_onoff_var"][2,1] == 1)
    @test (result_qc["solution"]["z_onoff_var"][1,2] == 1) || (result_qc["solution"]["z_onoff_var"][2,2] == 1)
    @test result_qc["solution"]["z_onoff_var"][4,3] == 1
    @test (result_qc["solution"]["z_onoff_var"][1,4] == 1) || (result_qc["solution"]["z_onoff_var"][2,4] == 1)
    @test (result_qc["solution"]["z_onoff_var"][1,5] == 1) || (result_qc["solution"]["z_onoff_var"][2,5] == 1)
    
end