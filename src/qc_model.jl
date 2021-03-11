function build_QCModel(data)
    m_qc = QuantumCircuitModel(data, JuMP.Model(), Dict{Symbol,Any}(), Dict{String,Any}())

    variable_QCModel(m_qc)
    constraint_QCModel(m_qc)
    objective_QCModel(m_qc)

    return m_qc
end

function run_QCModel(qcm::QuantumCircuitModel; optimizer=nothing)
    JuMP.set_optimizer(qcm.model, optimizer)
    # start_time = time()
    JuMP.optimize!(qcm.model)    
    qcm.result = Dict{String,Any}(
        "optimizer" => JuMP.solver_name(qcm.model),
        "termination_status" => JuMP.termination_status(qcm.model),
        "primal_status" => JuMP.primal_status(qcm.model),
        "dual_status" => JuMP.dual_status(qcm.model),
        "objective" => JuMP.objective_value(qcm.model),
        "objective_lb" => JuMP.objective_bound(qcm.model),
        "solve_time" => JuMP.solve_time(qcm.model),
        "solution" => JuMP.value.(qcm.variables[:x])
    )
    return qcm.result
end