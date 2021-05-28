function build_QCModel(data::Dict{String, Any}; 
                       model_type = "compact_formulation", 
                       commute_matrix_cuts = false)
    
    m_qc = QuantumCircuitModel(data)

    # convex-hull formulation per depth, but larger number of variables and constraints
    if model_type == "balas_formulation" 

        variable_QCModel(m_qc)
        constraint_QCModel(m_qc, commute_matrix_cuts)

    # minimal variables and constraints, but not a convex-hull formulation per depth
    elseif model_type == "compact_formulation" 
        
        variable_QCModel_compact(m_qc)
        constraint_QCModel_compact(m_qc, commute_matrix_cuts)

    end

    objective_QCModel(m_qc)

    return m_qc
end

function variable_QCModel(qcm::QuantumCircuitModel)
    variable_matrix_per_depth(qcm)
    variable_gates_onoff(qcm)
    variable_sequential_gate_products(qcm)
    variable_gate_products_copy(qcm)
    variable_gate_products_linearization(qcm)

    if qcm.data["decomposition_type"] == "approximate"
        variable_slack_for_feasibility(qcm)
    end

    return
end

function constraint_QCModel(qcm::QuantumCircuitModel, commute_matrix_cuts::Bool)
    constraint_single_gate_per_depth(qcm)
    constraint_disjunction_of_gates_per_depth(qcm)
    constraint_gate_initial_condition(qcm)
    constraint_gate_intermediate_products(qcm)
    constraint_gate_product_linearization(qcm)
    constraint_gate_target_condition(qcm)
    constraint_complex_to_real_symmetry(qcm)

    if commute_matrix_cuts
        constraint_commutative_gates(qcm)
    end

    return
end

function variable_QCModel_compact(qcm::QuantumCircuitModel)

    variable_gates_onoff(qcm)
    variable_sequential_gate_products(qcm)
    variable_gate_products_linearization(qcm)

    if qcm.data["decomposition_type"] == "approximate"
        variable_slack_for_feasibility(qcm)
    end

    return
end

function constraint_QCModel_compact(qcm::QuantumCircuitModel, commute_matrix_cuts::Bool)

    constraint_single_gate_per_depth(qcm)
    constraint_gate_initial_condition_compact(qcm)
    constraint_gate_intermediate_products_compact(qcm)
    constraint_gate_product_linearization(qcm)
    constraint_gate_target_condition_compact(qcm)
    constraint_complex_to_real_symmetry_compact(qcm)

    if commute_matrix_cuts
        constraint_commutative_gates(qcm)
    end

    return
end

function objective_QCModel(qcm::QuantumCircuitModel)
    
    if qcm.data["objective"] == "minimize_depth"
        objective_minimize_total_depth(qcm)
        
    elseif qcm.data["objective"] == "minimize_cnot"
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