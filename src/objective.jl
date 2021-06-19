#------------------------------------------------------------#
# Build the objective function for QuantumCircuitModel here  #
#------------------------------------------------------------#

function objective_minimize_total_depth(qcm::QuantumCircuitModel)
    
    n_r          = size(qcm.data["gates_real"])[1]
    n_c          = size(qcm.data["gates_real"])[2]
    depth        = qcm.data["depth"]
    M            = qcm.data["gates_real"]
    identity_idx = qcm.data["identity_idx"]

    decomposition_type = qcm.data["decomposition_type"]

    if !isempty(identity_idx)

        if decomposition_type == "exact"

            JuMP.@objective(qcm.model, Max, sum(qcm.variables[:z_onoff_var][n,d] for n in identity_idx, d=1:depth))

        elseif decomposition_type == "approximate"

            JuMP.@objective(qcm.model, Max, sum(qcm.variables[:z_onoff_var][n,d] for n in identity_idx, d=1:depth) - sum(qcm.variables[:slack_var][i,j]^2 for i=1:n_r, j=1:n_c))

        end

    else

        QCO.objective_feasibility(qcm)
        
    end                              
    
    return
end

function objective_feasibility(qcm::QuantumCircuitModel)
    n_r     = size(qcm.data["gates_real"])[1]
    n_c     = size(qcm.data["gates_real"])[2]

    decomposition_type = qcm.data["decomposition_type"]
    
    Memento.warn(_LOGGER, "Switching to a feasibility problem since Identity gate is not identified to be part of input elementary gates")

    if decomposition_type == "exact"

        JuMP.@objective(qcm.model, Min, 1)

    elseif decomposition_type == "approximate"

        JuMP.@objective(qcm.model, Min, sum(qcm.variables[:slack_var][i,j]^2 for i=1:n_r, j=1:n_c))

    end 

    return
end

function objective_minimize_cnot_gates(qcm::QuantumCircuitModel)
    
    n_r     = size(qcm.data["gates_real"])[1]
    n_c     = size(qcm.data["gates_real"])[2]
    depth   = qcm.data["depth"] 

    decomposition_type = qcm.data["decomposition_type"]
    
    # Note: The below objective minimizes both cnot_12 and cnot_21
    if !isempty(findall(x -> startswith(x, "cnot"), qcm.data["elementary_gates"]))
        cnot_id = Int64[]

        for i in keys(qcm.data["gates_dict"])
            if qcm.data["gates_dict"][i]["type"] == "cnot_12"
                push!(cnot_id, parse(Int64, i))
            end
            if qcm.data["gates_dict"][i]["type"] == "cnot_21"
                push!(cnot_id, parse(Int64, i))
            end
        end

        if decomposition_type == "exact"
            
            JuMP.@objective(qcm.model, Min, sum(qcm.variables[:z_onoff_var][n,d] for n in cnot_id, d=1:depth))

        elseif decomposition_type == "approximate"
            
            JuMP.@objective(qcm.model, Min, sum(qcm.variables[:z_onoff_var][n,d] for n in cnot_id, d=1:depth) + sum(qcm.variables[:slack_var][i,j]^2 for i=1:n_r, j=1:n_c))

        end

    else

        Memento.warn(_LOGGER, "Switching the objective to minimize the total depth since CNOT gate is not part of the input elementary gates")
        QCO.objective_minimize_total_depth(qcm)
    end

    return
end
 