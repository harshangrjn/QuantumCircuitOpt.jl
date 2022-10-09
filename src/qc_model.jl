#-----------------------------------------------------------#
# Build optimization model for circuit decomsposition here  #
#-----------------------------------------------------------#
import JuMP: MOI

function build_QCModel(data::Dict{String, Any}; options = nothing)
    
    m_qc = QCO.QuantumCircuitModel(data)

    # Update defaults to user-defined options
    if options !== nothing
        for i in keys(options)
            QCO.set_option(m_qc, i, options[i])
        end
    end

    QCO._catch_options_errors(m_qc)

    # convex-hull formulation per depth, but larger number of variables and constraints
    if m_qc.options.model_type == "balas_formulation" 

        QCO.variable_QCModel_balas(m_qc)
        QCO.constraint_QCModel_balas(m_qc)

    # minimal variables and constraints, but not a convex-hull formulation per depth
    elseif m_qc.options.model_type == "compact_formulation" 
        
        QCO.variable_QCModel_compact(m_qc)                                    
        QCO.constraint_QCModel_compact(m_qc)

    end

    QCO.objective_QCModel(m_qc)

    return m_qc
end

function variable_QCModel_balas(qcm::QuantumCircuitModel)
    QCO.variable_gates_per_depth(qcm)
    QCO.variable_gates_onoff(qcm)
    QCO.variable_sequential_gate_products(qcm)
    QCO.variable_gate_products_copy(qcm)
    QCO.variable_gate_products_linearization(qcm)

    if qcm.data["decomposition_type"] == "approximate"
        QCO.variable_slack_for_feasibility(qcm)
        QCO.variable_slack_var_outer_approximation(qcm)
    end

    QCO.variable_QCModel_valid(qcm)

    return
end

function constraint_QCModel_balas(qcm::QuantumCircuitModel)
    
    QCO.constraint_single_gate_per_depth(qcm)
    QCO.constraint_gates_onoff_per_depth(qcm)
    QCO.constraint_initial_gate_condition(qcm)
    QCO.constraint_intermediate_products(qcm)
    QCO.constraint_gate_product_linearization(qcm)
    QCO.constraint_target_gate_condition(qcm)
    QCO.constraint_cnot_gate_bounds(qcm)
    (qcm.data["decomposition_type"] == "approximate") && (QCO.constraint_slack_var_outer_approximation(qcm))
    (!qcm.data["are_gates_real"]) && (QCO.constraint_complex_to_real_symmetry(qcm))

    QCO.constraint_QCModel_valid(qcm)

    return
end

function variable_QCModel_compact(qcm::QuantumCircuitModel)

    QCO.variable_gates_onoff(qcm)
    QCO.variable_sequential_gate_products(qcm)
    QCO.variable_gate_products_linearization(qcm)

    if qcm.data["decomposition_type"] == "approximate"
        QCO.variable_slack_for_feasibility(qcm)
        QCO.variable_slack_var_outer_approximation(qcm)
    end

    QCO.variable_QCModel_valid(qcm)

    return
end

function variable_QCModel_valid(qcm::QuantumCircuitModel)

    if qcm.options.all_valid_constraints != -1

        if qcm.options.all_valid_constraints == 1 
            QCO.variable_binary_products(qcm)
            
        elseif qcm.options.all_valid_constraints == 0 
            qcm.options.unitary_constraints && QCO.variable_binary_products(qcm)
        end

    end

    return
end

function constraint_QCModel_compact(qcm::QuantumCircuitModel)

    QCO.constraint_single_gate_per_depth(qcm)
    QCO.constraint_initial_gate_condition_compact(qcm)
    QCO.constraint_intermediate_products_compact(qcm)
    QCO.constraint_gate_product_linearization(qcm)

    if qcm.data["decomposition_type"] == "optimal_global_phase"
        QCO.constraint_target_gate_condition_glphase(qcm)  
    else 
        QCO.constraint_target_gate_condition_compact(qcm)
    end
    
    QCO.constraint_cnot_gate_bounds(qcm)
    (qcm.data["decomposition_type"] == "approximate") && (QCO.constraint_slack_var_outer_approximation(qcm))

    # (!qcm.data["are_gates_real"]) && (QCO.constraint_complex_to_real_symmetry(qcm)) # seems to slow down MIP run times

    QCO.constraint_QCModel_valid(qcm)

    return
end

function constraint_QCModel_valid(qcm::QuantumCircuitModel)

    if qcm.options.all_valid_constraints != -1

        if qcm.options.all_valid_constraints == 1 
            QCO.constraint_commutative_gate_pairs(qcm)
            QCO.constraint_involutory_gates(qcm)
            QCO.constraint_redundant_gate_product_pairs(qcm)
            QCO.constraint_idempotent_gates(qcm)
            QCO.constraint_identity_gate_symmetry(qcm)
            QCO.constraint_convex_hull_complex_gates(qcm)
            QCO.constraint_unitary_property(qcm)
            
        elseif qcm.options.all_valid_constraints == 0 
            qcm.options.commute_gate_constraints            && QCO.constraint_commutative_gate_pairs(qcm)
            qcm.options.involutory_gate_constraints         && QCO.constraint_involutory_gates(qcm)
            qcm.options.redundant_gate_pair_constraints     && QCO.constraint_redundant_gate_product_pairs(qcm)
            qcm.options.idempotent_gate_constraints         && QCO.constraint_idempotent_gates(qcm)
            qcm.options.identity_gate_symmetry_constraints  && QCO.constraint_identity_gate_symmetry(qcm)
            qcm.options.convex_hull_gate_constraints        && QCO.constraint_convex_hull_complex_gates(qcm)
            qcm.options.unitary_constraints                 && QCO.constraint_unitary_property(qcm)
        end

    end

    return
end

function objective_QCModel(qcm::QuantumCircuitModel)
    
    if qcm.data["objective"] == "minimize_depth"
        QCO.objective_minimize_total_depth(qcm)
        
    elseif qcm.data["objective"] == "minimize_cnot"
        QCO.objective_minimize_cnot_gates(qcm)
    end

    return
end

""
function optimize_QCModel!(qcm::QuantumCircuitModel; optimizer=nothing)
    
    if qcm.options.relax_integrality
        JuMP.relax_integrality(qcm.model)
    end

    if JuMP.mode(qcm.model) != JuMP.DIRECT && optimizer !== nothing
        if JuMP.backend(qcm.model).optimizer === nothing
            JuMP.set_optimizer(qcm.model, optimizer)
        else
            Memento.warn(_LOGGER, "Model already contains optimizer, cannot use optimizer specified in `optimize_QCModel!`")
        end
    end

    JuMP.set_time_limit_sec(qcm.model, qcm.options.time_limit)

    if !qcm.options.optimizer_log
        JuMP.set_silent(qcm.model)
    end

    if JuMP.mode(qcm.model) != JuMP.DIRECT && JuMP.backend(qcm.model).optimizer === nothing
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
    
    qcm.result = QCO.build_QCModel_result(qcm, solve_time) 

    return qcm.result
end

function run_QCModel(params::Dict{String, Any}, 
                     qcm_optimizer::MOI.OptimizerWithAttributes; 
                     options = nothing)

    data = QCO.get_data(params)

    model_qc = QCO.build_QCModel(data, options = options)

    result_qc = QCO.optimize_QCModel!(model_qc, optimizer = qcm_optimizer)

    if model_qc.options.relax_integrality
        if result_qc["primal_status"] == MOI.FEASIBLE_POINT 
            Memento.info(_LOGGER, "Integrality-relaxed solutions can be found in the results dictionary")
        else
            Memento.info(_LOGGER, "Infeasible primal status for the integrality-relaxed problem")
        end
    else 
        if model_qc.options.visualize_solution
            QCO.visualize_solution(result_qc, data)
        end
    end

    return result_qc
end

function set_option(qcm::QuantumCircuitModel, s::Symbol, val)
    Base.setproperty!(qcm.options, s, val)
end

function _catch_options_errors(qcm::QuantumCircuitModel)
    
    if !(qcm.options.model_type in ["compact_formulation", "balas_formulation"])
        Memento.warn(_LOGGER, "Invalid model_type. Setting it to default value (compact_formulation).")
        QCO.set_option(qcm, :model_type, "compact_formulation")
    end

    if !(qcm.options.all_valid_constraints in [-1,0,1])
        Memento.warn(_LOGGER, "Invalid all_valid_constraints; choose a value ∈ [-1,0,1]. Setting it to default value of 0.")
        QCO.set_option(qcm, :all_valid_constraints, 0)
    end

end