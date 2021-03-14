function build_QCModel(data)
    m_qc = QuantumCircuitModel(data, JuMP.Model(), Dict{Symbol,Any}(), Dict{String,Any}())

    variable_QCModel(m_qc)
    constraint_QCModel(m_qc)
    objective_QCModel(m_qc)

    return m_qc
end

function variable_QCModel(qcm::QuantumCircuitModel)
    variable_matrix_per_depth(qcm)
    variable_gates_onoff(qcm)
    variable_sequential_gate_products(qcm)
    variable_gate_products_copy(qcm)
    variable_gate_products_linearization(qcm)
    return
end

function constraint_QCModel(qcm::QuantumCircuitModel)
    constraint_single_gate_per_depth(qcm)
    constraint_disjunction_of_gates_per_depth(qcm)
    constraint_gate_initial_condition(qcm)
    constraint_gate_intermediate_products(qcm)
    constraint_gate_product_linearization(qcm)
    constraint_gate_target_condition(qcm)
    constraint_complex_to_real_symmetry(qcm)
    return
end

function objective_QCModel(qcm::QuantumCircuitModel)
    if qcm.data["objective"] == "depth"
        objective_minimize_total_depth(qcm)
    elseif qcm.data["objective"][1:4] == "cnot"
        objective_minimize_cnot_gates(qcm)
    end
    return
end

""
function optimize_QCModel!(qcm::QuantumCircuitModel; optimizer=nothing)
    if qcm.data["relax_integrality"]
        JuMP.relax_integrality(qcm.model)
    end

    if JuMP.mode(qcm.model) != JuMP.DIRECT && optimizer !== nothing
        if qcm.model.moi_backend.state == MOI.Utilities.NO_OPTIMIZER
            JuMP.set_optimizer(qcm.model, optimizer)
        else
            Memento.warn(_LOGGER, "Model already contains optimizer, cannot use optimizer specified in `optimize_QCModel!`")
        end
    end

    if JuMP.mode(qcm.model) != JuMP.DIRECT && qcm.model.moi_backend.state == MOI.Utilities.NO_OPTIMIZER
        Memento.error(_LOGGER, "No optimizer specified in `optimize_QCModel!` or the given JuMP model.")
    end
    
    start_time = time()

    _, solve_time, solve_bytes_alloc, sec_in_gc = @timed JuMP.optimize!(qcm.model)

    try
        solve_time = JuMP.solve_time(qcm.model)
    catch
        Memento.warn(_LOGGER, "The given optimizer does not provide the SolveTime() attribute, falling back on @timed.  This is not a rigorous timing value.");
    end
    
    Memento.debug(_LOGGER, "JuMP model optimize time: $(time() - start_time)")
    
    qcm.result = build_QCModel_result(qcm, solve_time) 

    return qcm.result
end