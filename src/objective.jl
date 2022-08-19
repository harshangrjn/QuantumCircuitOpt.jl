#------------------------------------------------------------#
# Build the objective function for QuantumCircuitModel here  #
#------------------------------------------------------------#

function objective_minimize_total_depth(qcm::QuantumCircuitModel)
    
    n_r          = size(qcm.data["gates_real"])[1]
    n_c          = size(qcm.data["gates_real"])[2]
    depth        = qcm.data["maximum_depth"]
    identity_idx = qcm.data["identity_idx"]
    num_gates    = size(qcm.data["gates_real"])[3]

    decomposition_type = qcm.data["decomposition_type"]

    if !isempty(identity_idx)
        if decomposition_type == "exact_optimal"
            JuMP.@objective(qcm.model, Min, sum(qcm.variables[:z_bin_var][n,d] for n = 1:num_gates, d=1:depth if !(n in identity_idx)))

        elseif decomposition_type == "exact_feasible"
            QCO.objective_feasibility(qcm)

        elseif decomposition_type == "approximate"
            JuMP.@objective(qcm.model, Min, sum(qcm.variables[:z_bin_var][n,d] for n = 1:num_gates, d=1:depth if !(n in identity_idx)) + qcm.options.objective_slack_penalty * sum(qcm.variables[:slack_var][i,j]^2 for i=1:n_r, j=1:n_c))
        end
    else
        QCO.objective_feasibility(qcm)
    end                              
    
    return
end

function objective_feasibility(qcm::QuantumCircuitModel)

    decomposition_type = qcm.data["decomposition_type"]
    
    Memento.warn(_LOGGER, "Switching to a feasibility problem")

    if decomposition_type in ["exact_optimal", "exact_feasible"]
        # Feasibility objective
    elseif decomposition_type == "approximate"
        n_r     = size(qcm.data["gates_real"])[1]
        n_c     = size(qcm.data["gates_real"])[2]
        JuMP.@objective(qcm.model, Min, sum(qcm.variables[:slack_var][i,j]^2 for i=1:n_r, j=1:n_c))
    end 

    return
end

function objective_minimize_cnot_gates(qcm::QuantumCircuitModel)
    
    n_r      = size(qcm.data["gates_real"])[1]
    n_c      = size(qcm.data["gates_real"])[2]
    depth    = qcm.data["maximum_depth"] 
    cnot_idx = qcm.data["cnot_idx"]

    decomposition_type = qcm.data["decomposition_type"]
    
    if !isempty(qcm.data["cnot_idx"])

        if decomposition_type in ["exact_optimal", "exact_feasible"]
            JuMP.@objective(qcm.model, Min, sum(qcm.variables[:z_bin_var][n,d] for n in cnot_idx, d=1:depth))
            
        elseif decomposition_type == "approximate"            
            JuMP.@objective(qcm.model, Min, sum(qcm.variables[:z_bin_var][n,d] for n in cnot_idx, d=1:depth) + (qcm.options.objective_slack_penalty * sum(qcm.variables[:slack_var][i,j]^2 for i=1:n_r, j=1:n_c)))

        end

    elseif isempty(qcm.data["cnot_idx"]) && !isempty(qcm.data["identity_idx"])
        
        Memento.warn(_LOGGER, "Switching the objective to minimize the total depth since CNOT gate is not part of the input elementary gates")
        QCO.objective_minimize_total_depth(qcm)

    else

        Memento.warn(_LOGGER, "Switching to a feasibility problem since CNOT and Identity gates are not part of input elementary gates")
        QCO.objective_feasibility(qcm)

    end

    return
end
 