#--------------------------------------------------------#
# Build all constraints of the QuantumCircuitModel here  #
#--------------------------------------------------------#

function constraint_single_gate_per_depth(qcm::QuantumCircuitModel)
    n_gates = size(qcm.data["M_real"])[3]
    depth   = qcm.data["depth"]
    
    JuMP.@constraint(qcm.model, [d=1:depth], sum(qcm.variables[:z_onoff_var][n,d] for n=1:n_gates) == 1)
    
    return
end

function constraint_disjunction_of_gates_per_depth(qcm::QuantumCircuitModel)
    n_gates = size(qcm.data["M_real"])[3]
    depth   = qcm.data["depth"]

    JuMP.@constraint(qcm.model, [d=1:depth], qcm.variables[:M_var][:,:,d] .== 
                                    sum(qcm.variables[:z_onoff_var][n,d] * qcm.data["M_real"][:,:,n] for n=1:n_gates))
    
    return
end

function constraint_gate_initial_condition(qcm::QuantumCircuitModel)
    n_gates = size(qcm.data["M_real"])[3]

    JuMP.@constraint(qcm.model, sum(qcm.variables[:V_var][:,:,n,1] for n=1:n_gates) .== qcm.data["M_initial"])
    JuMP.@constraint(qcm.model, [n=1:n_gates], qcm.variables[:V_var][:,:,n,1] .== (qcm.variables[:z_onoff_var][n,1] .* qcm.data["M_initial"]))
    
    return
end

function constraint_gate_intermediate_products(qcm::QuantumCircuitModel)
    n_gates = size(qcm.data["M_real"])[3]
    depth   = qcm.data["depth"]

    JuMP.@constraint(qcm.model, [d=1:(depth-1)], qcm.variables[:U_var][:,:,d] .== 
                                sum(qcm.variables[:V_var][:,:,n,d] * qcm.data["M_real"][:,:,n] for n=1:n_gates))

    JuMP.@constraint(qcm.model, [d=2:depth], sum(qcm.variables[:V_var][:,:,n,d] for n=1:n_gates) .== 
                                                qcm.variables[:U_var][:,:,(d-1)])

    JuMP.@constraint(qcm.model, [d=1:(depth-1)], qcm.variables[:V_var][:,:,:,(d+1)] .== 
                                         qcm.variables[:zU_var][:,:,:,d])
    
    return
end

function constraint_gate_target_condition(qcm::QuantumCircuitModel)

    depth   = qcm.data["depth"]
    n_gates = size(qcm.data["M_real"])[3]
    decomposition_type = qcm.data["decomposition_type"]
    
    # For correct implementation of this, use MutableArithmetics.jl >= v0.2.11
    if decomposition_type == "exact"
    
        JuMP.@constraint(qcm.model, sum(qcm.variables[:V_var][:,:,n,depth] * qcm.data["M_real"][:,:,n] for n=1:n_gates) .== qcm.data["Target_real"])  
    
    elseif decomposition_type == "approximate"

        JuMP.@constraint(qcm.model, sum(qcm.variables[:V_var][:,:,n,depth] * qcm.data["M_real"][:,:,n] for n=1:n_gates) .== qcm.data["Target_real"][:,:] + qcm.variables[:slack_var][:,:])  
        
    end

    return
end

function constraint_complex_to_real_symmetry(qcm::QuantumCircuitModel)
    depth  = qcm.data["depth"]
    n_r    = size(qcm.data["M_real"])[1]
    n_c    = size(qcm.data["M_real"])[2]

    for i=1:2:n_r
        for j=1:2:n_c
            for d=1:depth
                if d <= (depth-1)
                    JuMP.@constraint(qcm.model, qcm.variables[:U_var][i,j,d]   ==  qcm.variables[:U_var][i+1,j+1,d])
                    JuMP.@constraint(qcm.model, qcm.variables[:U_var][i,j+1,d] == -qcm.variables[:U_var][i+1,j,d])
                end
                JuMP.@constraint(qcm.model, qcm.variables[:M_var][i,j,d]   ==  qcm.variables[:M_var][i+1,j+1,d])
                JuMP.@constraint(qcm.model, qcm.variables[:M_var][i,j+1,d] == -qcm.variables[:M_var][i+1,j,d])
            end
        end
    end
    
    return
end

function constraint_gate_product_linearization(qcm::QuantumCircuitModel)
    depth   = qcm.data["depth"]
    n_r     = size(qcm.data["M_real"])[1]
    n_c     = size(qcm.data["M_real"])[2]
    n_gates = size(qcm.data["M_real"])[3]

    for i=1:2:n_r
        for j=1:n_c
            for n=1:n_gates
                for d=1:(depth-1)
                    QCO.get_mccormick_relaxation(qcm.model, qcm.variables[:zU_var][i,j,n,d], qcm.variables[:U_var][i,j,d], qcm.variables[:z_onoff_var][n,(d+1)])
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
    n_gates = size(qcm.data["M_real"])[3]
    depth   = qcm.data["depth"]
    
    JuMP.@constraint(qcm.model, qcm.variables[:U_var][:,:,1] .== 
                                            qcm.data["M_initial"] * sum(qcm.variables[:z_onoff_var][n,1] * qcm.data["M_real"][:,:,n] for n=1:n_gates))
    
    return
end

function constraint_gate_intermediate_products_compact(qcm::QuantumCircuitModel)
    n_gates = size(qcm.data["M_real"])[3]
    depth   = qcm.data["depth"]
    
    if depth > 2
        JuMP.@constraint(qcm.model, [d=2:(depth-1)], qcm.variables[:U_var][:,:,d] .== 
                                            sum((qcm.variables[:zU_var][:,:,n,(d-1)]) * qcm.data["M_real"][:,:,n] for n=1:n_gates))
    end

    return
end

function constraint_gate_target_condition_compact(qcm::QuantumCircuitModel)
    n_r     = size(qcm.data["M_real"])[1]
    n_c     = size(qcm.data["M_real"])[2]
    depth   = qcm.data["depth"]
    n_gates = size(qcm.data["M_real"])[3]
    decomposition_type = qcm.data["decomposition_type"]

    zU_var = qcm.variables[:zU_var]
    
    # For correct implementation of this, use MutableArithmetics.jl >= v0.2.11
    if decomposition_type == "exact"
        
        JuMP.@constraint(qcm.model, sum(zU_var[:,:,n,(depth-1)] * qcm.data["M_real"][:,:,n] for n=1:n_gates) .== qcm.data["Target_real"][:,:])  
    
    elseif decomposition_type == "approximate"

        JuMP.@constraint(qcm.model, sum(zU_var[:,:,n,(depth-1)] * qcm.data["M_real"][:,:,n] for n=1:n_gates) .== qcm.data["Target_real"][:,:] + qcm.variables[:slack_var][:,:])    
    
    end
    
    return
end

function constraint_complex_to_real_symmetry_compact(qcm::QuantumCircuitModel)
    depth  = qcm.data["depth"]
    n_r    = size(qcm.data["M_real"])[1]
    n_c    = size(qcm.data["M_real"])[2]

    for i=1:2:n_r
        for j=1:2:n_c
            for d=1:(depth-1)
                JuMP.@constraint(qcm.model, qcm.variables[:U_var][i,j,d]   ==  qcm.variables[:U_var][i+1,j+1,d])
                JuMP.@constraint(qcm.model, qcm.variables[:U_var][i,j+1,d] == -qcm.variables[:U_var][i+1,j,d])
            end
        end
    end
    
    return
end

function constraint_commutative_gates(qcm::QuantumCircuitModel)
    
    depth  = qcm.data["depth"]

    commute_pairs, commute_triplets = QCO.get_commutative_gates(qcm.data["M_real"])
    z = qcm.variables[:z_onoff_var]

    if !isempty(commute_pairs)
        Memento.info(_LOGGER, "Detected $(length(commute_pairs)) input elementary gate pairs which commute")

        for i = 1:length(commute_pairs)
            JuMP.@constraint(qcm.model, [d=1:(depth-1)], sum(z[commute_pairs[i][k], d] for k=1:2)  + sum(z[commute_pairs[i][k], (d+1)] for k=1:2) <= 2)
        end

    end

    if !isempty(commute_triplets)
        Memento.info(_LOGGER, "Detected $(length(commute_triplets)) input elementary gate triplets which commute")

        for i = 1:length(commute_triplets)
            JuMP.@constraint(qcm.model, [d=1:(depth-1)], sum(z[commute_triplets[i][k], d] for k=1:3)  + sum(z[commute_triplets[i][k], (d+1)] for k=1:3) <= 3)
        end

    end

    return
end