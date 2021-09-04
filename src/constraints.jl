#--------------------------------------------------------#
# Build all constraints of the QuantumCircuitModel here  #
#--------------------------------------------------------#

function constraint_single_gate_per_depth(qcm::QuantumCircuitModel)

    num_gates = size(qcm.data["gates_real"])[3]
    depth     = qcm.data["depth"]
    
    JuMP.@constraint(qcm.model, [d=1:depth], sum(qcm.variables[:z_onoff_var][n,d] for n=1:num_gates) == 1)
    
    return
end

function constraint_disjunction_of_gates_per_depth(qcm::QuantumCircuitModel)

    num_gates = size(qcm.data["gates_real"])[3]
    depth     = qcm.data["depth"]

    JuMP.@constraint(qcm.model, [d=1:depth], qcm.variables[:M_var][:,:,d] .== 
                                    sum(qcm.variables[:z_onoff_var][n,d] * qcm.data["gates_real"][:,:,n] for n=1:num_gates))
    
    return
end

function constraint_gate_initial_condition(qcm::QuantumCircuitModel)

    num_gates = size(qcm.data["gates_real"])[3]

    JuMP.@constraint(qcm.model, sum(qcm.variables[:V_var][:,:,n,1] for n=1:num_gates) .== qcm.data["initial_gate"])
    JuMP.@constraint(qcm.model, [n=1:num_gates], qcm.variables[:V_var][:,:,n,1] .== (qcm.variables[:z_onoff_var][n,1] .* qcm.data["initial_gate"]))
    
    return
end

function constraint_gate_intermediate_products(qcm::QuantumCircuitModel)

    num_gates = size(qcm.data["gates_real"])[3]
    depth     = qcm.data["depth"]

    JuMP.@constraint(qcm.model, [d=1:(depth-1)], qcm.variables[:U_var][:,:,d] .== 
                                sum(qcm.variables[:V_var][:,:,n,d] * qcm.data["gates_real"][:,:,n] for n=1:num_gates))

    JuMP.@constraint(qcm.model, [d=2:depth], sum(qcm.variables[:V_var][:,:,n,d] for n=1:num_gates) .== 
                                                qcm.variables[:U_var][:,:,(d-1)])

    JuMP.@constraint(qcm.model, [d=1:(depth-1)], qcm.variables[:V_var][:,:,:,(d+1)] .== 
                                         qcm.variables[:zU_var][:,:,:,d])
    
    return
end

function constraint_gate_target_condition(qcm::QuantumCircuitModel)

    depth   = qcm.data["depth"]
    num_gates = size(qcm.data["gates_real"])[3]
    decomposition_type = qcm.data["decomposition_type"]
    
    # For correct implementation of this, use MutableArithmetics.jl >= v0.2.11
    if decomposition_type == "exact"
        JuMP.@constraint(qcm.model, sum(qcm.variables[:V_var][:,:,n,depth] * qcm.data["gates_real"][:,:,n] for n=1:num_gates) .== qcm.data["target_gate"])
    
    elseif decomposition_type == "approximate"
        JuMP.@constraint(qcm.model, sum(qcm.variables[:V_var][:,:,n,depth] * qcm.data["gates_real"][:,:,n] for n=1:num_gates) .== qcm.data["target_gate"][:,:] + qcm.variables[:slack_var][:,:])  

    end

    return
end

function constraint_complex_to_real_symmetry(qcm::QuantumCircuitModel)

    depth  = qcm.data["depth"]
    n_r    = size(qcm.data["gates_real"])[1]
    n_c    = size(qcm.data["gates_real"])[2]

    for i=1:2:n_r
        for j=1:2:n_c
            for d=1:depth

                if d <= (depth-1)
                    JuMP.@constraint(qcm.model, qcm.variables[:U_var][i,j,d]   ==  qcm.variables[:U_var][i+1,j+1,d])
                    JuMP.@constraint(qcm.model, qcm.variables[:U_var][i,j+1,d] == -qcm.variables[:U_var][i+1,j,d])
                end

                # JuMP.@constraint(qcm.model, qcm.variables[:M_var][i,j,d]   ==  qcm.variables[:M_var][i+1,j+1,d])
                # JuMP.@constraint(qcm.model, qcm.variables[:M_var][i,j+1,d] == -qcm.variables[:M_var][i+1,j,d])

            end
        end
    end
    
    return
end

function constraint_gate_product_linearization(qcm::QuantumCircuitModel)

    depth   = qcm.data["depth"]
    n_r     = size(qcm.data["gates_real"])[1]
    n_c     = size(qcm.data["gates_real"])[2]
    num_gates = size(qcm.data["gates_real"])[3]

    for i=1:2:n_r
        for j=1:n_c
            for n=1:num_gates
                for d=1:(depth-1)
                    
                    QCO.relaxation_bilinear(qcm.model, qcm.variables[:zU_var][i,j,n,d], qcm.variables[:U_var][i,j,d], qcm.variables[:z_onoff_var][n,(d+1)])
                    if isodd(j)
                        JuMP.@constraint(qcm.model, qcm.variables[:zU_var][i,j,n,d]   ==  qcm.variables[:zU_var][i+1,j+1,n,d])
                        JuMP.@constraint(qcm.model, qcm.variables[:zU_var][i,j+1,n,d] == -qcm.variables[:zU_var][i+1,j,n,d])
                    end

                end
            end
        end
    end
    
    return
end

function constraint_gate_initial_condition_compact(qcm::QuantumCircuitModel)

    num_gates = size(qcm.data["gates_real"])[3]
    
    JuMP.@constraint(qcm.model, qcm.variables[:U_var][:,:,1] .== 
                                            qcm.data["initial_gate"] * sum(qcm.variables[:z_onoff_var][n,1] * qcm.data["gates_real"][:,:,n] for n=1:num_gates))
    
    return
end

function constraint_gate_intermediate_products_compact(qcm::QuantumCircuitModel)

    num_gates = size(qcm.data["gates_real"])[3]
    depth   = qcm.data["depth"]
    
    JuMP.@constraint(qcm.model, [d=2:(depth-1)], qcm.variables[:U_var][:,:,d] .== 
                                        sum((qcm.variables[:zU_var][:,:,n,(d-1)]) * qcm.data["gates_real"][:,:,n] for n=1:num_gates))

    return
end

function constraint_gate_target_condition_compact(qcm::QuantumCircuitModel)

    depth   = qcm.data["depth"]
    num_gates = size(qcm.data["gates_real"])[3]
    decomposition_type = qcm.data["decomposition_type"]

    zU_var = qcm.variables[:zU_var]
    
    # For correct implementation of this, use MutableArithmetics.jl >= v0.2.11
    if decomposition_type == "exact"
        JuMP.@constraint(qcm.model, sum(zU_var[:,:,n,(depth-1)] * qcm.data["gates_real"][:,:,n] for n=1:num_gates) .== qcm.data["target_gate"][:,:])  
    
    elseif decomposition_type == "approximate"
        JuMP.@constraint(qcm.model, sum(zU_var[:,:,n,(depth-1)] * qcm.data["gates_real"][:,:,n] for n=1:num_gates) .== qcm.data["target_gate"][:,:] + qcm.variables[:slack_var][:,:])    
    
    end
    
    return
end

function constraint_complex_to_real_symmetry_compact(qcm::QuantumCircuitModel)

    depth  = qcm.data["depth"]
    n_r    = size(qcm.data["gates_real"])[1]
    n_c    = size(qcm.data["gates_real"])[2]

    JuMP.@constraint(qcm.model, [i=1:2:n_r, j=1:2:n_c, d=1:(depth-1)], qcm.variables[:U_var][i,j,d]   ==  qcm.variables[:U_var][i+1,j+1,d])
    JuMP.@constraint(qcm.model, [i=1:2:n_r, j=1:2:n_c, d=1:(depth-1)], qcm.variables[:U_var][i,j+1,d] == -qcm.variables[:U_var][i+1,j,d])
    
    return
end

function constraint_commutative_gate_pairs(qcm::QuantumCircuitModel)
    
    depth  = qcm.data["depth"]
    z_onoff_var  = qcm.variables[:z_onoff_var]

    commute_pairs, commute_pairs_prodIdentity = QCO.get_commutative_gate_pairs(qcm.data["gates_dict"])

    if !isempty(commute_pairs)
        (length(commute_pairs) == 1) && (Memento.info(_LOGGER, "Detected $(length(commute_pairs)) input elementary gate pair which commutes"))
        (length(commute_pairs) > 1)  && (Memento.info(_LOGGER, "Detected $(length(commute_pairs)) input elementary gate pairs which commute"))

        for i = 1:length(commute_pairs)
            JuMP.@constraint(qcm.model, [d=1:(depth-1)], z_onoff_var[commute_pairs[i][2], d] + z_onoff_var[commute_pairs[i][1], d+1] <= 1)
        end
    end

    if !isempty(commute_pairs_prodIdentity)
        (length(commute_pairs_prodIdentity) == 1) && (Memento.info(_LOGGER, "Detected $(length(commute_pairs_prodIdentity)) input elementary gate pair whose product is Identity"))
        (length(commute_pairs_prodIdentity) > 1)  && (Memento.info(_LOGGER, "Detected $(length(commute_pairs_prodIdentity)) input elementary gate pairs whose product is Identity"))

        for i = 1:length(commute_pairs_prodIdentity)
            JuMP.@constraint(qcm.model, [d=1:(depth-1)], z_onoff_var[commute_pairs_prodIdentity[i][2], d] + z_onoff_var[commute_pairs_prodIdentity[i][1], d+1] <= 1)
            JuMP.@constraint(qcm.model, [d=1:(depth-1)], z_onoff_var[commute_pairs_prodIdentity[i][1], d] + z_onoff_var[commute_pairs_prodIdentity[i][2], d+1] <= 1)
        end
    end

    return
end

function constraint_involutory_gates(qcm::QuantumCircuitModel)

    gates_dict  = qcm.data["gates_dict"]
    depth       = qcm.data["depth"]
    z_onoff_var = qcm.variables[:z_onoff_var]

    involutory_gates = QCO.get_involutory_gates(gates_dict)
    
    if !isempty(involutory_gates)
        (length(involutory_gates) == 1) && (Memento.info(_LOGGER, "Detected $(length(involutory_gates)) involutory elementary gate"))
        (length(involutory_gates) > 1)  && (Memento.info(_LOGGER, "Detected $(length(involutory_gates)) involutory elementary gates"))

        for i = 1:length(involutory_gates)
            JuMP.@constraint(qcm.model, [d=1:(depth-1)], z_onoff_var[involutory_gates[i], d] + z_onoff_var[involutory_gates[i], d+1] <= 1)
        end
    end

    return
end

function constraint_redundant_gate_product_pairs(qcm::QuantumCircuitModel)

    gates_dict  = qcm.data["gates_dict"]
    depth       = qcm.data["depth"]
    z_onoff_var = qcm.variables[:z_onoff_var]

    redundant_pairs = QCO.get_redundant_gate_product_pairs(gates_dict)
    
    if !isempty(redundant_pairs)
        (length(redundant_pairs) == 1) && (Memento.info(_LOGGER, "Detected $(length(redundant_pairs)) redundant input elementary gate pair"))
        (length(redundant_pairs) > 1)  && (Memento.info(_LOGGER, "Detected $(length(redundant_pairs)) redundant input elementary gate pairs"))

        for i = 1:length(redundant_pairs)
            JuMP.@constraint(qcm.model, [d=1:(depth-1)], z_onoff_var[redundant_pairs[i][2], d] + z_onoff_var[redundant_pairs[i][1], d+1] <= 1)
        end
    end

    return
end

function constraint_idempotent_gates(qcm::QuantumCircuitModel)

    gates_dict  = qcm.data["gates_dict"]
    depth       = qcm.data["depth"]
    z_onoff_var = qcm.variables[:z_onoff_var]

    idempotent_gates = QCO.get_idempotent_gates(gates_dict)
    
    if !isempty(idempotent_gates)
        (length(idempotent_gates) == 1) && (Memento.info(_LOGGER, "Detected $(length(idempotent_gates)) idempotent elementary gate"))
        (length(idempotent_gates) > 1)  && (Memento.info(_LOGGER, "Detected $(length(idempotent_gates)) idempotent elementary gates"))

        for i = 1:length(idempotent_gates)
            JuMP.@constraint(qcm.model, [d=1:(depth-1)], z_onoff_var[idempotent_gates[i], d] + z_onoff_var[idempotent_gates[i], d+1] <= 1)
        end
    end

    return
end