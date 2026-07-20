#------------------------------------------------------------#
# Build the objective function for QuantumCircuitModel here  #
#------------------------------------------------------------#

function objective_minimize_total_depth(qcm::QuantumCircuitModel)
    
    max_depth    = qcm.data["maximum_depth"]
    identity_idx = qcm.data["identity_idx"]
    num_gates    = size(qcm.data["gates_real"])[3]

    decomposition_type = qcm.data["decomposition_type"]

    if !isempty(identity_idx)
        if decomposition_type in ["exact_optimal", "optimal_global_phase"]
            JuMP.@objective(qcm.model, Min, sum(qcm.variables[:z_bin_var][n,d] for n = 1:num_gates, d=1:max_depth if !(n in identity_idx)))

        elseif decomposition_type == "exact_feasible"
            QCO.objective_feasibility(qcm)

        elseif decomposition_type == "approximate"
            JuMP.@objective(qcm.model, Min, sum(qcm.variables[:z_bin_var][n,d] for n = 1:num_gates, d=1:max_depth if !(n in identity_idx))
            + qcm.options.objective_slack_penalty * sum(qcm.variables[:slack_var_oa]))

        elseif decomposition_type == "sum-approximate"
            JuMP.@objective(qcm.model, Min, sum(qcm.variables[:slack_var_oa]))
        
        elseif decomposition_type == "infinity-approximate"
            JuMP.@objective(qcm.model, Min, qcm.variables[:slack_var_oa_one])

        elseif decomposition_type == "fidelity-approximate"
            JuMP.@objective(qcm.model, Max, qcm.variables[:w_fidelity_var])

        elseif decomposition_type == "linear-fidelity-approximate"
            JuMP.@objective(qcm.model, Max, qcm.variables[:c_fidelity_var])
        end
    else
        QCO.objective_feasibility(qcm)
    end                              
    
    return
end

function objective_feasibility(qcm::QuantumCircuitModel)

    decomposition_type = qcm.data["decomposition_type"]

    if decomposition_type in ["exact_optimal", "exact_feasible", "optimal_global_phase"]
        # Feasibility objective
        @_warn("Switching to a feasibility problem")
        
    elseif decomposition_type == "approximate"
        JuMP.@objective(qcm.model, Min, sum(qcm.variables[:slack_var_oa]))
    end 

    return
end

function objective_minimize_specific_gates(
    qcm::QuantumCircuitModel,
    gate_type::String
    )
    
    max_depth = qcm.data["maximum_depth"]
    decomposition_type = qcm.data["decomposition_type"]
    
    if gate_type == "cnot_gate"
        gate_idx = qcm.data["cnot_idx"]
        gate_name = "CNot"
    elseif gate_type == "T_gate"
        gate_idx = qcm.data["T_idx"]
        gate_name = "T"
    else
        @_error("Unsupported gate type: $gate_type")
    end
    
    if !isempty(gate_idx)
        if decomposition_type in ["exact_optimal", "exact_feasible", "optimal_global_phase"]
            JuMP.@objective(qcm.model, Min, sum(qcm.variables[:z_bin_var][n,d] for n in gate_idx, d=1:max_depth))
            
        elseif decomposition_type == "approximate"            
            JuMP.@objective(qcm.model, Min, sum(qcm.variables[:z_bin_var][n,d] for n in gate_idx, d=1:max_depth) 
            + (qcm.options.objective_slack_penalty * sum(qcm.variables[:slack_var_oa])))
        end
    else
        @_error("$gate_name gate not found in input elementary gates")
    end

    return
end

function objective_minimize_total_layers(qcm::QuantumCircuitModel)
    
    max_depth      = qcm.data["maximum_depth"]
    maximum_layers = qcm.data["maximum_layers"]
    identity_idx   = qcm.data["identity_idx"]
    num_gates      = size(qcm.data["gates_real"])[3]

    decomposition_type = qcm.data["decomposition_type"]

    if !isempty(identity_idx)
        if decomposition_type in ["exact_optimal", "optimal_global_phase"]
            JuMP.@objective(qcm.model, Min, sum(l*qcm.variables[:layer_bin_var][max_depth,l] for l=1:maximum_layers ))

        elseif decomposition_type == "exact_feasible"
            @_error("`exact_feasible` decomposition is not supported for this objective")
            
        elseif decomposition_type in ["approximate", "fidelity-approximate"]
            @_error("`approximate` decomposition is not supported for this objective")
            
        end
    else
        QCO.objective_feasibility(qcm)
    end                              
    
    return
end