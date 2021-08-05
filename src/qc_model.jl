function build_QCModel(data::Dict{String, Any}; 
                       model_type = "compact_formulation",
                       all_valid_constraints = 0, 
                       commute_gate_constraints = true,
                       involutory_gate_constraints = true,
                       redundant_gate_pair_constraints = true,
                       idempotent_gate_constraints = false)

    if !(all_valid_constraints in [-1,0,1])
        Memento.warn(_LOGGER, "Invalid all_valid_constraints; choose a value ∈ [-1,0,1]. Setting it to default value of 0.")
        all_valid_constraints = 0
    end

    if !(model_type in ["compact_formulation", "balas_formulation"])
        Memento.warn(_LOGGER, "Invalid model_type. Setting it to default value (compact_formulation).")
        model_type = "compact_formulation"
    end
    
    m_qc = QuantumCircuitModel(data)

    # convex-hull formulation per depth, but larger number of variables and constraints
    if model_type == "balas_formulation" 

        QCO.variable_QCModel(m_qc)
        QCO.constraint_QCModel(m_qc, 
                               all_valid_constraints,
                               commute_gate_constraints, 
                               involutory_gate_constraints,
                               redundant_gate_pair_constraints,
                               idempotent_gate_constraints,
                              )

    # minimal variables and constraints, but not a convex-hull formulation per depth
    elseif model_type == "compact_formulation" 
        
        QCO.variable_QCModel_compact(m_qc)
        QCO.constraint_QCModel_compact(m_qc, 
                                   all_valid_constraints,
                                   commute_gate_constraints, 
                                   involutory_gate_constraints,
                                   redundant_gate_pair_constraints,
                                   idempotent_gate_constraints
                                   )

    end

    QCO.objective_QCModel(m_qc)

    return m_qc
end

function variable_QCModel(qcm::QuantumCircuitModel)
    QCO.variable_matrix_per_depth(qcm)
    QCO.variable_gates_onoff(qcm)
    QCO.variable_sequential_gate_products(qcm)
    QCO.variable_gate_products_copy(qcm)
    QCO.variable_gate_products_linearization(qcm)

    if qcm.data["decomposition_type"] == "approximate"
        QCO.variable_slack_for_feasibility(qcm)
    end

    return
end

function constraint_QCModel(qcm::QuantumCircuitModel, 
                            all_valid_constraints::Int64,
                            commute_gate_constraints::Bool, 
                            involutory_gate_constraints::Bool,
                            redundant_gate_pair_constraints::Bool,
                            idempotent_gate_constraints::Bool)
    
    QCO.constraint_single_gate_per_depth(qcm)
    QCO.constraint_disjunction_of_gates_per_depth(qcm)
    QCO.constraint_gate_initial_condition(qcm)
    QCO.constraint_gate_intermediate_products(qcm)
    QCO.constraint_gate_product_linearization(qcm)
    QCO.constraint_gate_target_condition(qcm)
    QCO.constraint_complex_to_real_symmetry(qcm)

    QCO.constraint_QCModel_valid(qcm, 
                                all_valid_constraints,
                                commute_gate_constraints, 
                                involutory_gate_constraints,
                                redundant_gate_pair_constraints,
                                idempotent_gate_constraints)

    return
end

function variable_QCModel_compact(qcm::QuantumCircuitModel)

    QCO.variable_gates_onoff(qcm)
    QCO.variable_sequential_gate_products(qcm)
    QCO.variable_gate_products_linearization(qcm)

    if qcm.data["decomposition_type"] == "approximate"
        QCO.variable_slack_for_feasibility(qcm)
    end

    return
end

function constraint_QCModel_compact(qcm::QuantumCircuitModel, 
                                    all_valid_constraints::Int64,
                                    commute_gate_constraints::Bool, 
                                    involutory_gate_constraints::Bool,
                                    redundant_gate_pair_constraints::Bool,
                                    idempotent_gate_constraints::Bool)

    QCO.constraint_single_gate_per_depth(qcm)
    QCO.constraint_gate_initial_condition_compact(qcm)
    QCO.constraint_gate_intermediate_products_compact(qcm)
    QCO.constraint_gate_product_linearization(qcm)
    QCO.constraint_gate_target_condition_compact(qcm)
    QCO.constraint_complex_to_real_symmetry_compact(qcm)

    QCO.constraint_QCModel_valid(qcm, 
                                 all_valid_constraints,
                                 commute_gate_constraints, 
                                 involutory_gate_constraints,
                                 redundant_gate_pair_constraints,
                                 idempotent_gate_constraints)

    return
end

function constraint_QCModel_valid(qcm::QuantumCircuitModel,
                                  all_valid_constraints::Int64,
                                  commute_gate_constraints::Bool, 
                                  involutory_gate_constraints::Bool,
                                  redundant_gate_pair_constraints::Bool,
                                  idempotent_gate_constraints::Bool)

    if all_valid_constraints != -1

        if all_valid_constraints == 1 
            QCO.constraint_commutative_gate_pairs(qcm)
            QCO.constraint_involutory_gates(qcm)
            QCO.constraint_redundant_gate_product_pairs(qcm)
            QCO.constraint_idempotent_gates(qcm)
            
        elseif all_valid_constraints == 0 
            commute_gate_constraints        && QCO.constraint_commutative_gate_pairs(qcm)
            involutory_gate_constraints     && QCO.constraint_involutory_gates(qcm)
            redundant_gate_pair_constraints && QCO.constraint_redundant_gate_product_pairs(qcm)
            idempotent_gate_constraints     && QCO.constraint_idempotent_gates(qcm)
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
    
    if qcm.data["relax_integrality"]
        JuMP.relax_integrality(qcm.model)
    end

    if "time_limit" in keys(qcm.data)
        JuMP.set_time_limit_sec(qcm.model, qcm.data["time_limit"])
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
    
    qcm.result = QCO.build_QCModel_result(qcm, solve_time) 

    return qcm.result
end

function run_QCModel(params::Dict{String, Any}, 
                     qcm_optimizer::MOI.OptimizerWithAttributes; 
                     model_type = "compact_formulation", 
                     all_valid_constraints = 0,
                     commute_gate_constraints = true, 
                     involutory_gate_constraints = true, 
                     redundant_gate_pair_constraints = true,
                     idempotent_gate_constraints = false,
                     visualize_solution=true, 
                     eliminate_identical_gates = true)

    data = QCO.get_data(params, 
                        eliminate_identical_gates = eliminate_identical_gates)

    model_qc  = QCO.build_QCModel(data, 
                                  model_type = model_type, 
                                  all_valid_constraints = all_valid_constraints,
                                  commute_gate_constraints = commute_gate_constraints, 
                                  involutory_gate_constraints = involutory_gate_constraints,
                                  redundant_gate_pair_constraints = redundant_gate_pair_constraints,
                                  idempotent_gate_constraints = idempotent_gate_constraints)

    result_qc = QCO.optimize_QCModel!(model_qc, optimizer = qcm_optimizer)

    if visualize_solution
        QCO.visualize_solution(result_qc, data)
    end

    return result_qc
end