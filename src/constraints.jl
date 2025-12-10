#--------------------------------------------------------#
# Build all constraints of the QuantumCircuitModel here  #
#--------------------------------------------------------#

###################################   NECESSARY CONSTRAINTS   ###################################

function constraint_single_gate_per_depth(qcm::QuantumCircuitModel)

    num_gates = size(qcm.data["gates_real"])[3]
    max_depth = qcm.data["maximum_depth"]
    
    JuMP.@constraint(qcm.model, [d=1:max_depth], sum(qcm.variables[:z_bin_var][n,d] for n=1:num_gates) == 1)
    
    return
end

function constraint_gates_onoff_per_depth(qcm::QuantumCircuitModel)

    num_gates = size(qcm.data["gates_real"])[3]
    max_depth = qcm.data["maximum_depth"]

    JuMP.@constraint(qcm.model, [d=1:max_depth], qcm.variables[:G_var][:,:,d] .== 
                                    sum(qcm.variables[:z_bin_var][n,d] * qcm.data["gates_real"][:,:,n] for n=1:num_gates))
    
    return
end

function constraint_initial_gate_condition(qcm::QuantumCircuitModel)

    num_gates = size(qcm.data["gates_real"])[3]

    JuMP.@constraint(qcm.model, sum(qcm.variables[:V_var][:,:,n,1] for n=1:num_gates) .== qcm.data["initial_gate"])
    JuMP.@constraint(qcm.model, [n=1:num_gates], qcm.variables[:V_var][:,:,n,1] .== (qcm.variables[:z_bin_var][n,1] .* qcm.data["initial_gate"]))
    
    return
end

function constraint_intermediate_products(qcm::QuantumCircuitModel)

    num_gates = size(qcm.data["gates_real"])[3]
    max_depth = qcm.data["maximum_depth"]

    JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], qcm.variables[:U_var][:,:,d] .== 
                                sum(qcm.variables[:V_var][:,:,n,d] * qcm.data["gates_real"][:,:,n] for n=1:num_gates))

    JuMP.@constraint(qcm.model, [d=2:max_depth], sum(qcm.variables[:V_var][:,:,n,d] for n=1:num_gates) .== 
                                                qcm.variables[:U_var][:,:,(d-1)])

    JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], qcm.variables[:V_var][:,:,:,(d+1)] .== 
                                         qcm.variables[:zU_var][:,:,:,d])
    
    return
end

function constraint_gate_product_linearization(qcm::QuantumCircuitModel)

    max_depth      = qcm.data["maximum_depth"]
    n_r            = size(qcm.data["gates_real"])[1]
    n_c            = size(qcm.data["gates_real"])[2]
    num_gates      = size(qcm.data["gates_real"])[3]
    are_gates_real = qcm.data["are_gates_real"]

    U_var  = qcm.variables[:U_var]
    zU_var = qcm.variables[:zU_var]
    z_bin_var = qcm.variables[:z_bin_var]

    i_val = 2 - are_gates_real
    U_var_bound_tol = 1E-8

    for i=1:i_val:n_r, j=1:n_c, n=1:num_gates, d=1:(max_depth-1)
        U_var_l = JuMP.lower_bound(U_var[i,j,d])
        U_var_u = JuMP.upper_bound(U_var[i,j,d])
        
        if !(isapprox(abs(U_var_u - U_var_l), 0, atol = 2*U_var_bound_tol))
            QCO.relaxation_bilinear(qcm.model, zU_var[i,j,n,d], U_var[i,j,d], z_bin_var[n,(d+1)])
        else
            JuMP.@constraint(qcm.model, zU_var[i,j,n,d] == (U_var_l + U_var_u)/2 * z_bin_var[n,(d+1)])
        end

        if !are_gates_real && isodd(j)
            JuMP.@constraint(qcm.model, zU_var[i,j,n,d]   ==  zU_var[i+1,j+1,n,d])
            JuMP.@constraint(qcm.model, zU_var[i,j+1,n,d] == -zU_var[i+1,j,n,d])
        end
    end
    
    return
end

function constraint_initial_gate_condition_compact(qcm::QuantumCircuitModel)

    num_gates = size(qcm.data["gates_real"])[3]
    
    JuMP.@constraint(qcm.model, qcm.variables[:U_var][:,:,1] .== 
                                            qcm.data["initial_gate"] * sum(qcm.variables[:z_bin_var][n,1] * qcm.data["gates_real"][:,:,n] for n=1:num_gates))
    
    return
end

function constraint_intermediate_products_compact(qcm::QuantumCircuitModel)

    num_gates = size(qcm.data["gates_real"])[3]
    max_depth = qcm.data["maximum_depth"]

    JuMP.@constraint(qcm.model, [d=2:max_depth], qcm.variables[:U_var][:,:,d] .== 
                                        sum((qcm.variables[:zU_var][:,:,n,(d-1)]) * qcm.data["gates_real"][:,:,n] for n=1:num_gates))
                                        
    return
end

function constraint_layer_assigments(qcm::QuantumCircuitModel)

    max_depth      = qcm.data["maximum_depth"]
    maximum_layers = qcm.data["maximum_layers"]
    
    #Some elements have to be 0: layer number <= depth number
    for d = 1:(maximum_layers-1)
        JuMP.@constraint(qcm.model, [l=(d+1):maximum_layers], qcm.variables[:layer_bin_var][d,l] == 0)
    end

    #The first element is 1
    JuMP.@constraint(qcm.model, qcm.variables[:layer_bin_var][1, 1] == 1)
    
    #Each position assigned o exactly one layer
    JuMP.@constraint(qcm.model, [d=1:max_depth], sum(qcm.variables[:layer_bin_var][d,l] for l=1:maximum_layers) == 1)

    #Layers stepwise monotonic increase
    JuMP.@constraint(qcm.model, [d=1:max_depth-1], sum(qcm.variables[:layer_bin_var][d,l]*l for l=1:maximum_layers) <= sum(qcm.variables[:layer_bin_var][d+1,l]*l for l=1:maximum_layers))
    JuMP.@constraint(qcm.model, [d=1:max_depth-1], sum(qcm.variables[:layer_bin_var][d+1,l]*l for l=1:maximum_layers) <= 1 + sum(qcm.variables[:layer_bin_var][d,l]*l for l=1:maximum_layers))

    return
end

function constraint_single_gate_per_qubit_per_layer(qcm::QuantumCircuitModel)

    max_depth      = qcm.data["maximum_depth"]
    maximum_layers = qcm.data["maximum_layers"]
    num_gates      = size(qcm.data["gates_real"])[3]
    num_qubits     = qcm.data["num_qubits"]
    gates_dict     = qcm.data["gates_dict"]
    
    #Some elements have to be 0: layer number <= depth number
    for d = 1:(maximum_layers-1)
        JuMP.@constraint(qcm.model, [n=1:num_gates, l=(d+1):maximum_layers], qcm.variables[:delta_bin_var][n,d,l] == 0)
    end
    
    #auxiliary variable multiplication constraints 
    JuMP.@constraint(qcm.model, [n=1:num_gates, d=1:max_depth, l=1:maximum_layers], qcm.variables[:delta_bin_var][n,d,l] <= qcm.variables[:z_bin_var][n,d])
    JuMP.@constraint(qcm.model, [n=1:num_gates, d=1:max_depth, l=1:maximum_layers], qcm.variables[:delta_bin_var][n,d,l] <= qcm.variables[:layer_bin_var][d,l])
    JuMP.@constraint(qcm.model, [n=1:num_gates, d=1:max_depth, l=1:maximum_layers], qcm.variables[:delta_bin_var][n,d,l] >= qcm.variables[:z_bin_var][n,d] + qcm.variables[:layer_bin_var][d,l] -1)

    #each qubit used only once per layer
    JuMP.@constraint(qcm.model, [q=1:num_qubits, l=1:maximum_layers], sum(sum(qcm.variables[:delta_bin_var][n,d,l] for d=1:max_depth) for n in _get_all_gates_for_qubit(gates_dict, q)) <= 1)
    
    return
end

###################################   TARGET CONSTRAINTS   ###################################

function constraint_target_gate_condition(qcm::QuantumCircuitModel)

    max_depth          = qcm.data["maximum_depth"]
    num_gates          = size(qcm.data["gates_real"])[3]
    decomposition_type = qcm.data["decomposition_type"]
    
    if decomposition_type in ["exact_optimal", "exact_feasible"]
        JuMP.@constraint(qcm.model, sum(qcm.variables[:V_var][:,:,n,max_depth] * qcm.data["gates_real"][:,:,n] for n=1:num_gates) 
                                    .== qcm.data["target_gate"])
    
    elseif decomposition_type == "approximate"
        JuMP.@constraint(qcm.model, sum(qcm.variables[:V_var][:,:,n,max_depth] * qcm.data["gates_real"][:,:,n] for n=1:num_gates) 
                                    .== qcm.data["target_gate"] + qcm.variables[:slack_var][:,:])

    end

                                        
    return
end

function constraint_layer_assigments(qcm::QuantumCircuitModel)

    max_depth      = qcm.data["maximum_depth"]
    maximum_layers = qcm.data["maximum_layers"]
    
    #Some elements have to be 0: layer number <= depth number
    for d = 1:(maximum_layers-1)
        JuMP.@constraint(qcm.model, [l=(d+1):maximum_layers], qcm.variables[:layer_bin_var][d,l] == 0)
    end

    #The first element is 1
    JuMP.@constraint(qcm.model, qcm.variables[:layer_bin_var][1, 1] == 1)
    
    #Each position assigned o exactly one layer
    JuMP.@constraint(qcm.model, [d=1:max_depth], sum(qcm.variables[:layer_bin_var][d,l] for l=1:maximum_layers) == 1)

    #Layers stepwise monotonic increase
    JuMP.@constraint(qcm.model, [d=1:max_depth-1], sum(qcm.variables[:layer_bin_var][d,l]*l for l=1:maximum_layers) <= sum(qcm.variables[:layer_bin_var][d+1,l]*l for l=1:maximum_layers))
    JuMP.@constraint(qcm.model, [d=1:max_depth-1], sum(qcm.variables[:layer_bin_var][d+1,l]*l for l=1:maximum_layers) <= 1 + sum(qcm.variables[:layer_bin_var][d,l]*l for l=1:maximum_layers))

    return
end

function constraint_single_gate_per_qubit_per_layer(qcm::QuantumCircuitModel)

    max_depth      = qcm.data["maximum_depth"]
    maximum_layers = qcm.data["maximum_layers"]
    num_gates      = size(qcm.data["gates_real"])[3]
    num_qubits     = qcm.data["num_qubits"]
    gates_dict     = qcm.data["gates_dict"]
    
    #Some elements have to be 0: layer number <= depth number
    for d = 1:(maximum_layers-1)
        JuMP.@constraint(qcm.model, [n=1:num_gates, l=(d+1):maximum_layers], qcm.variables[:delta_bin_var][n,d,l] == 0)
    end
    
    #auxiliary variable multiplication constraints 
    JuMP.@constraint(qcm.model, [n=1:num_gates, d=1:max_depth, l=1:maximum_layers], qcm.variables[:delta_bin_var][n,d,l] <= qcm.variables[:z_bin_var][n,d])
    JuMP.@constraint(qcm.model, [n=1:num_gates, d=1:max_depth, l=1:maximum_layers], qcm.variables[:delta_bin_var][n,d,l] <= qcm.variables[:layer_bin_var][d,l])
    JuMP.@constraint(qcm.model, [n=1:num_gates, d=1:max_depth, l=1:maximum_layers], qcm.variables[:delta_bin_var][n,d,l] >= qcm.variables[:z_bin_var][n,d] + qcm.variables[:layer_bin_var][d,l] -1)

    #each qubit used only once per layer
    JuMP.@constraint(qcm.model, [q=1:num_qubits, l=1:maximum_layers], sum(sum(qcm.variables[:delta_bin_var][n,d,l] for d=1:max_depth) for n in _get_all_gates_for_qubit(gates_dict, q)) <= 1)
    
    return
end

###################################   TARGET CONSTRAINTS   ###################################

function constraint_target_gate_condition(qcm::QuantumCircuitModel)

    max_depth          = qcm.data["maximum_depth"]
    num_gates          = size(qcm.data["gates_real"])[3]
    decomposition_type = qcm.data["decomposition_type"]
    
    if decomposition_type in ["exact_optimal", "exact_feasible"]
        JuMP.@constraint(qcm.model, sum(qcm.variables[:V_var][:,:,n,max_depth] * qcm.data["gates_real"][:,:,n] for n=1:num_gates) 
                                    .== qcm.data["target_gate"])
    
    elseif decomposition_type == "approximate"
        JuMP.@constraint(qcm.model, sum(qcm.variables[:V_var][:,:,n,max_depth] * qcm.data["gates_real"][:,:,n] for n=1:num_gates) 
                                    .== qcm.data["target_gate"] + qcm.variables[:slack_var][:,:])

    end

    return
end

function constraint_target_gate_condition_compact(qcm::QuantumCircuitModel)

    max_depth = qcm.data["maximum_depth"]
    decomposition_type = qcm.data["decomposition_type"]

    U_var = qcm.variables[:U_var]
    
    if decomposition_type in ["exact_optimal", "exact_feasible"]
        JuMP.@constraint(qcm.model, U_var[:,:,max_depth] .== qcm.data["target_gate"][:,:])
    
    elseif decomposition_type in ["approximate", "infinity-approximate", "sum-approximate"]
        JuMP.@constraint(qcm.model, U_var[:,:,max_depth] .== qcm.data["target_gate"][:,:] + qcm.variables[:slack_var][:,:])
    end
    
    return
end

# Helper function for global phase constraints
function _real_phase_difference(U_var, target_gate, ref_r, ref_c, i, j, depth)
    return U_var[i, j, depth]*target_gate[ref_r, ref_c] -
           U_var[ref_r, ref_c, depth]*target_gate[i, j]
end

# Helper function for global phase constraints
function _complex_phase_differences(U_var, target_gate, ref_r, ref_c, i, j, depth)
    realLHS = U_var[i, j, depth]   * target_gate[ref_r, ref_c] -
              U_var[i, j+1, depth] * target_gate[ref_r, ref_c+1]
    realRHS = U_var[ref_r, ref_c, depth]   * target_gate[i, j] -
              U_var[ref_r, ref_c+1, depth] * target_gate[i, j+1]

    imagLHS = U_var[i, j, depth]   * target_gate[ref_r, ref_c+1] +
              U_var[i, j+1, depth] * target_gate[ref_r, ref_c]
    imagRHS = U_var[ref_r, ref_c, depth]   * target_gate[i, j+1] +
              U_var[ref_r, ref_c+1, depth] * target_gate[i, j]

    return realLHS, realRHS, imagLHS, imagRHS
end

function constraint_target_gate_condition_glphase1_cs(qcm::QuantumCircuitModel)

    max_depth      = qcm.data["maximum_depth"]
    n_r            = size(qcm.data["gates_real"])[1]
    n_c            = size(qcm.data["gates_real"])[2]
    target_gate    = qcm.data["target_gate"]
    U_var          = qcm.variables[:U_var]
    phase_real_var = qcm.variables[:phase_real_var]
    phase_complex_var = qcm.variables[:phase_complex_var]

    if qcm.data["are_gates_real"]

        JuMP.@constraint(qcm.model, phase_complex_var == 0)

        for i in 1:n_r, j in 1:n_c
            JuMP.@constraint(qcm.model,  U_var[i,j,max_depth] == target_gate[i,j] * phase_real_var)
        end

    else
        for i in 1:2:n_r, j in 1:2:n_c
            JuMP.@constraint(qcm.model, U_var[i,j,max_depth] == phase_real_var * target_gate[i,j] - phase_complex_var * target_gate[i+1,j])
            JuMP.@constraint(qcm.model, U_var[i,j+1,max_depth] == - phase_complex_var * target_gate[i,j] - phase_real_var * target_gate[i+1,j])
            JuMP.@constraint(qcm.model, U_var[i+1,j,max_depth] == phase_complex_var * target_gate[i,j] + phase_real_var * target_gate[i+1,j])
            JuMP.@constraint(qcm.model, U_var[i+1,j+1,max_depth] == phase_real_var * target_gate[i,j] - phase_complex_var * target_gate[i+1,j])

        end

    end

    return
end

function constraint_target_gate_condition_glphase2_cs(qcm::QuantumCircuitModel)

    max_depth      = qcm.data["maximum_depth"]
    n_r            = size(qcm.data["gates_real"])[1]
    n_c            = size(qcm.data["gates_real"])[2]
    target_gate    = qcm.data["target_gate"]
    U_var          = qcm.variables[:U_var]
    phase_real_var = qcm.variables[:phase_real_var]
    phase_complex_var = qcm.variables[:phase_complex_var]

    if qcm.data["are_gates_real"]

        JuMP.@constraint(qcm.model, phase_complex_var == 0)
        
        for i in 1:n_r, j in 1:n_c
            if !(i==j)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[j,k] for k=1:n_r) == 0)
            end    
        end

        for i in 1:n_r
            JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[i,k] for k=1:n_r) == phase_real_var)
        end


    else
        for i in 1:2:n_r, j in 1:2:n_c
            if !(i==j)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[j,k] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[j+1,k] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[j,k] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[j+1,k] for k=1:n_r) == 0)
            end
        end

        for i in 1:2:n_r
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[i,k] for k=1:n_r) == phase_real_var)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[i+1,k] for k=1:n_r) == -phase_complex_var)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[i,k] for k=1:n_r) == phase_complex_var)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[i+1,k] for k=1:n_r) == phase_real_var)
        end
            
    end

    return
end

function constraint_target_gate_condition_glphase1_cs(qcm::QuantumCircuitModel)

    max_depth      = qcm.data["maximum_depth"]
    n_r            = size(qcm.data["gates_real"])[1]
    n_c            = size(qcm.data["gates_real"])[2]
    target_gate    = qcm.data["target_gate"]
    U_var          = qcm.variables[:U_var]
    phase_real_var = qcm.variables[:phase_real_var]
    phase_complex_var = qcm.variables[:phase_complex_var]

    if qcm.data["are_gates_real"]

        JuMP.@constraint(qcm.model, phase_complex_var == 0)

        for i in 1:n_r, j in 1:n_c
            JuMP.@constraint(qcm.model,  U_var[i,j,max_depth] == target_gate[i,j] * phase_real_var)
        end

    else
        for i in 1:2:n_r, j in 1:2:n_c
            JuMP.@constraint(qcm.model, U_var[i,j,max_depth] == phase_real_var * target_gate[i,j] - phase_complex_var * target_gate[i+1,j])
            JuMP.@constraint(qcm.model, U_var[i,j+1,max_depth] == - phase_complex_var * target_gate[i,j] - phase_real_var * target_gate[i+1,j])
            JuMP.@constraint(qcm.model, U_var[i+1,j,max_depth] == phase_complex_var * target_gate[i,j] + phase_real_var * target_gate[i+1,j])
            JuMP.@constraint(qcm.model, U_var[i+1,j+1,max_depth] == phase_real_var * target_gate[i,j] - phase_complex_var * target_gate[i+1,j])

        end

    end

    return
end

function constraint_target_gate_condition_glphase2_cs(qcm::QuantumCircuitModel)

    max_depth      = qcm.data["maximum_depth"]
    n_r            = size(qcm.data["gates_real"])[1]
    n_c            = size(qcm.data["gates_real"])[2]
    target_gate    = qcm.data["target_gate"]
    U_var          = qcm.variables[:U_var]
    phase_real_var = qcm.variables[:phase_real_var]
    phase_complex_var = qcm.variables[:phase_complex_var]

    if qcm.data["are_gates_real"]

        JuMP.@constraint(qcm.model, phase_complex_var == 0)
        
        for i in 1:n_r, j in 1:n_c
            if !(i==j)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[j,k] for k=1:n_r) == 0)
            end    
        end

        for i in 1:n_r
            JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[i,k] for k=1:n_r) == phase_real_var)
        end


    else
        for i in 1:2:n_r, j in 1:2:n_c
            if !(i==j)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[j,k] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[j+1,k] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[j,k] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[j+1,k] for k=1:n_r) == 0)
            end
        end

        for i in 1:2:n_r
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[i,k] for k=1:n_r) == phase_real_var)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[i+1,k] for k=1:n_r) == -phase_complex_var)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[i,k] for k=1:n_r) == phase_complex_var)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[i+1,k] for k=1:n_r) == phase_real_var)
        end
            
    end

    return
end

function constraint_target_gate_condition_glphase2_cs_only_diag(qcm::QuantumCircuitModel)

    max_depth      = qcm.data["maximum_depth"]
    n_r            = size(qcm.data["gates_real"])[1]
    n_c            = size(qcm.data["gates_real"])[2]
    target_gate    = qcm.data["target_gate"]
    U_var          = qcm.variables[:U_var]
    phase_real_var = qcm.variables[:phase_real_var]
    phase_complex_var = qcm.variables[:phase_complex_var]

    if qcm.data["are_gates_real"]

        JuMP.@constraint(qcm.model, phase_complex_var == 0)
        
        for i in 1:n_r
            JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[i,k] for k=1:n_r) == phase_real_var)
        end


    else
        
        for i in 1:2:n_r
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[i,k] for k=1:n_r) == phase_real_var)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[i+1,k] for k=1:n_r) == -phase_complex_var)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[i,k] for k=1:n_r) == phase_complex_var)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[i+1,k] for k=1:n_r) == phase_real_var)
        end
            
    end

    return
end

function constraint_target_gate_condition_glphase1_star(qcm::QuantumCircuitModel)
    max_depth      = qcm.data["maximum_depth"]
    n_r            = size(qcm.data["gates_real"])[1]
    n_c            = size(qcm.data["gates_real"])[2]
    target_gate    = qcm.data["target_gate"]
    U_var          = qcm.variables[:U_var]

    if qcm.data["are_gates_real"]
        ref_r, ref_c = QCO._get_nonzero_idx_of_complex_matrix(Array{ComplexF64}(target_gate))

        for i in 1:n_r, j in 1:n_c
            if QCO.is_zero(target_gate[i, j])
                JuMP.@constraint(qcm.model, U_var[i, j, max_depth] == 0)

            elseif !((i == ref_r) && (j == ref_c))
                LHS = QCO._real_phase_difference(U_var, target_gate, ref_r, ref_c, i, j, max_depth)
                JuMP.@constraint(qcm.model, LHS == 0)
            end
        end

    else
        ref_r, ref_c = QCO._get_nonzero_idx_of_complex_to_real_matrix(target_gate)
        for i in 1:2:n_r, j in 1:2:n_c
            if QCO.is_zero(target_gate[i, j]) && QCO.is_zero(target_gate[i, j+1])
                JuMP.@constraint(qcm.model, U_var[i, j, max_depth] == 0.0)
                JuMP.@constraint(qcm.model, U_var[i, j+1, max_depth] == 0.0)

            elseif !((i == ref_r) && (j == ref_c))
                realLHS, realRHS, imagLHS, imagRHS = QCO._complex_phase_differences(U_var, target_gate, ref_r, ref_c, i, j, max_depth)
        
                JuMP.@constraint(qcm.model, realLHS == realRHS)
                JuMP.@constraint(qcm.model, imagLHS == imagRHS)
            end
        end
    end

    return
end

function constraint_target_gate_condition_glphase2_a2a(qcm::QuantumCircuitModel)
    max_depth      = qcm.data["maximum_depth"]
    n_r            = size(qcm.data["gates_real"])[1]
    n_c            = size(qcm.data["gates_real"])[2]
    target_gate    = qcm.data["target_gate"]
    U_var          = qcm.variables[:U_var]

    if qcm.data["are_gates_real"]
        
        for i in 1:n_r, j in 1:n_c
            if !(i==j)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[j,k] for k=1:n_r) == 0)
            end
        end

        for i in 1:n_r, j in i+1:n_r
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[i,k] for k=1:n_r) == sum( U_var[j,k,max_depth] * target_gate[j,k] for k=1:n_r))
        end

    else
        for i in 1:2:n_r, j in 1:2:n_c
            if !(i==j)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[j,k] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[j+1,k] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[j,k] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[j+1,k] for k=1:n_r) == 0)
            end
        end

        for i in 1:2:n_r, j in (i+2):2:n_c
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[i,k] for k=1:n_r) == sum( U_var[j,k,max_depth] * target_gate[j,k] for k=1:n_r))
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[i,k] for k=1:n_r) == sum( U_var[j+1,k,max_depth] * target_gate[j,k] for k=1:n_r))
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[i+1,k] for k=1:n_r) == sum( U_var[j,k,max_depth] * target_gate[j+1,k] for k=1:n_r))
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[i+1,k] for k=1:n_r) == sum( U_var[j+1,k,max_depth] * target_gate[j+1,k] for k=1:n_r))
        end
            
    end

    return
end

function constraint_target_gate_condition_glphase2_star(qcm::QuantumCircuitModel)
    max_depth      = qcm.data["maximum_depth"]
    n_r            = size(qcm.data["gates_real"])[1]
    n_c            = size(qcm.data["gates_real"])[2]
    target_gate    = qcm.data["target_gate"]
    U_var          = qcm.variables[:U_var]

    if qcm.data["are_gates_real"]
        
        for i in 1:n_r, j in 1:n_c
            if !(i==j)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[j,k] for k=1:n_r) == 0)
            end
        end

        for i in 2:n_r
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[i,k] for k=1:n_r) == sum( U_var[1,k,max_depth] * target_gate[1,k] for k=1:n_r))
        end

    else
        for i in 1:2:n_r, j in 1:2:n_c
            if !(i==j)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[j,k] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[j+1,k] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[j,k] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[j+1,k] for k=1:n_r) == 0)
            end
        end

        for i in 3:2:n_r
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[i,k] for k=1:n_r) == sum( U_var[1,k,max_depth] * target_gate[1,k] for k=1:n_r))
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[i,k] for k=1:n_r) == sum( U_var[1+1,k,max_depth] * target_gate[1,k] for k=1:n_r))
                JuMP.@constraint(qcm.model, sum( U_var[i,k,max_depth] * target_gate[i+1,k] for k=1:n_r) == sum( U_var[1,k,max_depth] * target_gate[1+1,k] for k=1:n_r))
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,max_depth] * target_gate[i+1,k] for k=1:n_r) == sum( U_var[1+1,k,max_depth] * target_gate[1+1,k] for k=1:n_r))
        end
            
    end

    return
end

function constraint_slack_var_outer_approximation(qcm::QuantumCircuitModel)
    slack_var    = qcm.variables[:slack_var]
    slack_var_oa = qcm.variables[:slack_var_oa]
    
    # Number of under-estimators for the quadratic function
    num_points = 5
    #eps = 1.0/9.0
    

    for i=1:size(slack_var)[1], j=1:size(slack_var)[2]
        lb = JuMP.lower_bound(slack_var[i,j])
        ub = JuMP.upper_bound(slack_var[i,j])
        if !(isapprox(lb, ub, atol = 1E-6))
            oa_points = range(lb, ub, num_points)
            #oa_points = range(-eps, eps, num_points)
            JuMP.@constraint(qcm.model, [k=1:length(oa_points)], 
                            slack_var_oa[i,j] >= 2*slack_var[i,j]*oa_points[k] - (oa_points[k])^2)
            if !(0 in oa_points)
                #JuMP.@constraint(qcm.model, slack_var_oa[i,j] >= 0)
            end
        else
            mid_point = (lb+ub)/2
            JuMP.@constraint(qcm.model, slack_var_oa[i,j] >= 2*slack_var[i,j]*mid_point - (mid_point)^2)
        end
    end

    return
end

function constraint_slack_var_outer_approximation_quadratic(qcm::QuantumCircuitModel)
    n_r          = size(qcm.data["gates_real"])[1]
    n_c          = size(qcm.data["gates_real"])[2]
    slack_var    = qcm.variables[:slack_var]
    slack_var_oa = qcm.variables[:slack_var_oa]
    
    JuMP.@constraint(qcm.model, [i=1:n_r, j=1:n_c], slack_var_oa[i,j] >= slack_var[i,j]^2)
 
    return
end

function constraint_slack_var_sum_approximation(qcm::QuantumCircuitModel)
    n_r          = size(qcm.data["gates_real"])[1]
    n_c          = size(qcm.data["gates_real"])[2]
    slack_var    = qcm.variables[:slack_var]
    slack_var_oa = qcm.variables[:slack_var_oa]
    
    JuMP.@constraint(qcm.model, [i=1:n_r, j=1:n_c], slack_var_oa[i,j] >= slack_var[i,j])
    JuMP.@constraint(qcm.model, [i=1:n_r, j=1:n_c], slack_var_oa[i,j] >= - slack_var[i,j])
 
    return
end

function constraint_slack_var_infinity_approximation(qcm::QuantumCircuitModel)
    n_r          = size(qcm.data["gates_real"])[1]
    n_c          = size(qcm.data["gates_real"])[2]
    slack_var    = qcm.variables[:slack_var]
    slack_var_oa_one = qcm.variables[:slack_var_oa_one]
    
    JuMP.@constraint(qcm.model, [i=1:n_r, j=1:n_c], slack_var_oa_one >= slack_var[i,j])
    JuMP.@constraint(qcm.model, [i=1:n_r, j=1:n_c], slack_var_oa_one >= - slack_var[i,j])
    
    return
end

function constraint_slack_var_infinity_approximation_quadratic(qcm::QuantumCircuitModel)
    n_r          = size(qcm.data["gates_real"])[1]
    n_c          = size(qcm.data["gates_real"])[2]
    slack_var    = qcm.variables[:slack_var]
    slack_var_oa_one = qcm.variables[:slack_var_oa_one]
    
    JuMP.@constraint(qcm.model, [i=1:n_r, j=1:n_c], slack_var_oa_one >= slack_var[i,j]^2)
    
    return
end

function constraint_fidelity(qcm::QuantumCircuitModel)

    # w = c^2 + s^2
    JuMP.@constraint(qcm.model, qcm.variables[:w_fidelity_var] ==  qcm.variables[:c_real_var] * qcm.variables[:c_real_var]  +  qcm.variables[:s_complex_var] * qcm.variables[:s_complex_var] )
    
    return
end

function constraint_trace_for_approximate(qcm::QuantumCircuitModel)

    max_depth      = qcm.data["maximum_depth"]
    n_r            = size(qcm.data["gates_real"])[1]
    n_c            = size(qcm.data["gates_real"])[2]
    target_gate    = qcm.data["target_gate"]
    U_var          = qcm.variables[:U_var]

    P = [0 -1; 1 0]
    n_blocks = Int(n_c) ÷ 2
    Jc = kron(Matrix{Float64}(I, n_blocks, n_blocks), P)   # size n_c × n_c
    
    JuMP.@constraint(qcm.model, qcm.variables[:c_real_var] == 1/n_r* sum(sum( target_gate[i,j] * U_var[i,j,max_depth]  for i=1:n_r) for j=1:n_c ))
    # JuMP.@constraint(qcm.model,  qcm.variables[:s_complex_var] == (1/n_r) * sum( U_var[i,k,max_depth] * (sum(target_gate[i,l] * Jc[l,k] for l=1:n_c)) for i=1:n_r, k=1:n_c ))
    JuMP.@constraint(qcm.model,  qcm.variables[:s_complex_var] == - (1/n_r) * sum(sum( Jc[l,k] * target_gate[i,k] * U_var[i,l,max_depth] for l=1:n_c) for i=1:n_r, k=1:n_c ))
    
    return
end

function constraint_real_trace_for_approximate(qcm::QuantumCircuitModel)

    max_depth      = qcm.data["maximum_depth"]
    n_r            = size(qcm.data["gates_real"])[1]
    n_c            = size(qcm.data["gates_real"])[2]
    target_gate    = qcm.data["target_gate"]
    U_var          = qcm.variables[:U_var]
    
    JuMP.@constraint(qcm.model, qcm.variables[:c_fidelity_var] == 1/n_r* sum(sum( U_var[i,j,max_depth] * target_gate[i,j] for i=1:n_r) for j=1:n_c ))
       
    return
end

###################################   SYMMETRIES    ###################################

function constraint_identity_gate_symmetry(qcm::QuantumCircuitModel)

    gates_dict = qcm.data["gates_dict"]
    max_depth  = qcm.data["maximum_depth"]
    z_bin_var  = qcm.variables[:z_bin_var]

    identity_idx = []
    for i=1:length(keys(gates_dict))
        if "Identity" in gates_dict["$i"]["type"]
            push!(identity_idx, i)
        end
    end
    
    if !isempty(identity_idx)
        for i = 1:length(identity_idx)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], z_bin_var[identity_idx[i], d] <= z_bin_var[identity_idx[i], d+1])
        end
    end

    return
end

function constraint_commutative_gate_pairs(qcm::QuantumCircuitModel)
    
    objective      = qcm.data["objective"]
    max_depth      = qcm.data["maximum_depth"]
    maximum_layers = qcm.data["maximum_layers"]
    z_bin_var      = qcm.variables[:z_bin_var]

    commute_pairs, commute_pairs_prodIdentity = QCO.get_commutative_gate_pairs(qcm.data["gates_dict"], qcm.data["decomposition_type"])

    if !isempty(commute_pairs)
        (length(commute_pairs) == 1) && (Memento.info(_LOGGER, "Detected $(length(commute_pairs)) input elementary gate pair which commutes"))
        (length(commute_pairs) > 1)  && (Memento.info(_LOGGER, "Detected $(length(commute_pairs)) input elementary gate pairs which commute"))

        for i = 1:length(commute_pairs)
            if objective in ["minimize_layers"]
                JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], z_bin_var[commute_pairs[i][2], d] + z_bin_var[commute_pairs[i][1], d+1] <= 1 
                +  sum(qcm.variables[:layer_bin_var][d+1,l]*l for l=1:maximum_layers) -  sum(qcm.variables[:layer_bin_var][d,l]*l for l=1:maximum_layers))
            else
                JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], z_bin_var[commute_pairs[i][2], d] + z_bin_var[commute_pairs[i][1], d+1] <= 1)
            end
        end
    end

    if !isempty(commute_pairs_prodIdentity)
        (length(commute_pairs_prodIdentity) == 1) && (Memento.info(_LOGGER, "Detected $(length(commute_pairs_prodIdentity)) input elementary gate pair whose product is Identity"))
        (length(commute_pairs_prodIdentity) > 1)  && (Memento.info(_LOGGER, "Detected $(length(commute_pairs_prodIdentity)) input elementary gate pairs whose product is Identity"))

        for i = 1:length(commute_pairs_prodIdentity)            
            JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], z_bin_var[commute_pairs_prodIdentity[i][2], d] + z_bin_var[commute_pairs_prodIdentity[i][1], d+1] <= 1)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], z_bin_var[commute_pairs_prodIdentity[i][1], d] + z_bin_var[commute_pairs_prodIdentity[i][2], d+1] <= 1)
        end
    end

    return
end

function constraint_equivalent_gate_pairs(qcm::QuantumCircuitModel)

    gates_dict = qcm.data["gates_dict"]
    max_depth  = qcm.data["maximum_depth"]
    z_bin_var  = qcm.variables[:z_bin_var]

    if qcm.data["objective"] in  ["minimize_depth",  "minimize_cnot", "minimize_T"]
        equivalent_pairs = QCO.get_equivalent_gate_pairs(gates_dict, qcm.data["decomposition_type"])
    end

    if qcm.data["objective"] in  ["minimize_layers"]
        equivalent_pairs = QCO.get_equivalent_gate_pairs_layers(gates_dict, qcm.data["decomposition_type"])
    end

    if !isempty(equivalent_pairs)
        (length(equivalent_pairs) == 1) && (Memento.info(_LOGGER, "Detected $(length(equivalent_pairs)) equivalent input elementary gate pair"))
        (length(equivalent_pairs) > 1)  && (Memento.info(_LOGGER, "Detected $(length(equivalent_pairs)) equivalent input elementary gate pairs"))

        for i = 1:length(equivalent_pairs)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], z_bin_var[equivalent_pairs[i][1], d] + z_bin_var[equivalent_pairs[i][2], d+1] <= 1)
        end
    end

    return
end

function constraint_equivalent_gate_triplets(qcm::QuantumCircuitModel)

    gates_dict = qcm.data["gates_dict"]
    max_depth  = qcm.data["maximum_depth"]
    z_bin_var  = qcm.variables[:z_bin_var]

    equivalent_triplets = QCO.get_equivalent_gate_triplets(gates_dict, qcm.data["decomposition_type"])
    
    if !isempty(equivalent_triplets)
        (length(equivalent_triplets) == 1) && (Memento.info(_LOGGER, "Detected $(length(equivalent_triplets)) equivalent input elementary gate triplet"))
        (length(equivalent_triplets) > 1)  && (Memento.info(_LOGGER, "Detected $(length(equivalent_triplets)) equivalent input elementary gate triplets"))

        for i = 1:length(equivalent_triplets)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-2)], z_bin_var[equivalent_triplets[i][1], d] + z_bin_var[equivalent_triplets[i][2], d+1] + z_bin_var[equivalent_triplets[i][3], d+2] <= 2)
        end
    end

    return
end

###################################   REDUNDANCIES   ###################################

function constraint_involutory_gates(qcm::QuantumCircuitModel)

    gates_dict = qcm.data["gates_dict"]
    max_depth  = qcm.data["maximum_depth"]
    z_bin_var  = qcm.variables[:z_bin_var]

    involutory_gates = QCO.get_involutory_gates(gates_dict)
    #println("Involuntary", involutory_gates )
    
    if !isempty(involutory_gates)
        (length(involutory_gates) == 1) && (Memento.info(_LOGGER, "Detected $(length(involutory_gates)) involutory elementary gate"))
        (length(involutory_gates) > 1)  && (Memento.info(_LOGGER, "Detected $(length(involutory_gates)) involutory elementary gates"))

        for i = 1:length(involutory_gates)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], z_bin_var[involutory_gates[i], d] + z_bin_var[involutory_gates[i], d+1] <= 1)
        end
    end

    return
end

function constraint_idempotent_gates(qcm::QuantumCircuitModel)

    gates_dict = qcm.data["gates_dict"]
    max_depth  = qcm.data["maximum_depth"]
    z_bin_var  = qcm.variables[:z_bin_var]

    idempotent_gates = QCO.get_idempotent_gates(gates_dict, qcm.data["decomposition_type"])
    
    if !isempty(idempotent_gates)
        (length(idempotent_gates) == 1) && (Memento.info(_LOGGER, "Detected $(length(idempotent_gates)) idempotent elementary gate"))
        (length(idempotent_gates) > 1)  && (Memento.info(_LOGGER, "Detected $(length(idempotent_gates)) idempotent elementary gates"))
        println("Idempotent:", idempotent_gates )

        for i = 1:length(idempotent_gates)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], z_bin_var[idempotent_gates[i], d] + z_bin_var[idempotent_gates[i], d+1] <= 1)
        for i = 1:length(idempotent_gates)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], z_bin_var[idempotent_gates[i], d] + z_bin_var[idempotent_gates[i], d+1] <= 1)
        end
    end

    return
end

function constraint_redundant_gate_products(qcm::QuantumCircuitModel)
function constraint_redundant_gate_products(qcm::QuantumCircuitModel)

    gates_dict = qcm.data["gates_dict"]
    max_depth  = qcm.data["maximum_depth"]
    z_bin_var  = qcm.variables[:z_bin_var]
    num_qubits = qcm.data["num_qubits"]

    redundant_pairs, redundant_triplets, redundant_quadruplets, redundant_quintuplets = QCO.get_redundant_gate_products(gates_dict, qcm.data["decomposition_type"])
    
    if !isempty(redundant_pairs)
        (length(redundant_pairs) == 1) && (Memento.info(_LOGGER, "Detected $(length(redundant_pairs)) redundant input elementary gate pair"))
        (length(redundant_pairs) > 1)  && (Memento.info(_LOGGER, "Detected $(length(redundant_pairs)) redundant input elementary gate pairs"))

        for i = 1:length(redundant_pairs)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], z_bin_var[redundant_pairs[i][1], d] + z_bin_var[redundant_pairs[i][2], d+1] <= 1)
        end
    end

    redundant_pairs, redundant_triplets, redundant_quadruplets, redundant_quintuplets = QCO.get_redundant_gate_products(gates_dict, qcm.data["decomposition_type"])
    
    if !isempty(redundant_pairs)
        (length(redundant_pairs) == 1) && (Memento.info(_LOGGER, "Detected $(length(redundant_pairs)) redundant input elementary gate pair"))
        (length(redundant_pairs) > 1)  && (Memento.info(_LOGGER, "Detected $(length(redundant_pairs)) redundant input elementary gate pairs"))

        for i = 1:length(redundant_pairs)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], z_bin_var[redundant_pairs[i][1], d] + z_bin_var[redundant_pairs[i][2], d+1] <= 1)
        end
    end

    if !isempty(redundant_triplets)
        (length(redundant_triplets) == 1) && (Memento.info(_LOGGER, "Detected $(length(redundant_triplets)) redundant input elementary gate triplet"))
        (length(redundant_triplets) == 1) && (Memento.info(_LOGGER, "Detected $(length(redundant_triplets)) redundant input elementary gate triplet"))
        (length(redundant_triplets) > 1)  && (Memento.info(_LOGGER, "Detected $(length(redundant_triplets)) redundant input elementary gate triplets"))

        for i = 1:length(redundant_triplets)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-2)], z_bin_var[redundant_triplets[i][1], d] + z_bin_var[redundant_triplets[i][2], d+1] + z_bin_var[redundant_triplets[i][3], d+2]<= 2)
        end
    end

    if !isempty(redundant_quadruplets) && max_depth >= 4
        (length(redundant_quadruplets) == 1) && (Memento.info(_LOGGER, "Detected $(length(redundant_quadruplets)) redundant input elementary gate quadruplet"))
        (length(redundant_quadruplets) > 1)  && (Memento.info(_LOGGER, "Detected $(length(redundant_quadruplets)) redundant input elementary gate quadruplets"))

        for i = 1:length(redundant_quadruplets)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-3)], z_bin_var[redundant_quadruplets[i][1], d] + z_bin_var[redundant_quadruplets[i][2], d+1] + z_bin_var[redundant_quadruplets[i][3], d+2] + z_bin_var[redundant_quadruplets[i][4], d+3] <= 3)
        end
    end

    if !isempty(redundant_quintuplets) && max_depth >= 5
        (length(redundant_quintuplets) == 1) && (Memento.info(_LOGGER, "Detected $(length(redundant_quintuplets)) redundant input elementary gate quintuplet"))
        (length(redundant_quintuplets) > 1)  && (Memento.info(_LOGGER, "Detected $(length(redundant_quintuplets)) redundant input elementary gate quintuplets"))

        for i = 1:length(redundant_quintuplets)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-4)], z_bin_var[redundant_quintuplets[i][1], d] + z_bin_var[redundant_quintuplets[i][2], d+1] + z_bin_var[redundant_quintuplets[i][3], d+2] + z_bin_var[redundant_quintuplets[i][4], d+3] + z_bin_var[redundant_quintuplets[i][5], d+4] <= 4)
        end
    end

    if !isempty(redundant_quadruplets) && max_depth >= 4
        (length(redundant_quadruplets) == 1) && (Memento.info(_LOGGER, "Detected $(length(redundant_quadruplets)) redundant input elementary gate quadruplet"))
        (length(redundant_quadruplets) > 1)  && (Memento.info(_LOGGER, "Detected $(length(redundant_quadruplets)) redundant input elementary gate quadruplets"))

        for i = 1:length(redundant_quadruplets)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-3)], z_bin_var[redundant_quadruplets[i][1], d] + z_bin_var[redundant_quadruplets[i][2], d+1] + z_bin_var[redundant_quadruplets[i][3], d+2] + z_bin_var[redundant_quadruplets[i][4], d+3] <= 3)
        end
    end

    if !isempty(redundant_quintuplets) && max_depth >= 5
        (length(redundant_quintuplets) == 1) && (Memento.info(_LOGGER, "Detected $(length(redundant_quintuplets)) redundant input elementary gate quintuplet"))
        (length(redundant_quintuplets) > 1)  && (Memento.info(_LOGGER, "Detected $(length(redundant_quintuplets)) redundant input elementary gate quintuplets"))

        for i = 1:length(redundant_quintuplets)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-4)], z_bin_var[redundant_quintuplets[i][1], d] + z_bin_var[redundant_quintuplets[i][2], d+1] + z_bin_var[redundant_quintuplets[i][3], d+2] + z_bin_var[redundant_quintuplets[i][4], d+3] + z_bin_var[redundant_quintuplets[i][5], d+4] <= 4)
        end
    end

    return
end

function constraint_redundant_gate_product_pairs(qcm::QuantumCircuitModel)

    gates_dict = qcm.data["gates_dict"]
    max_depth  = qcm.data["maximum_depth"]
    z_bin_var  = qcm.variables[:z_bin_var]

    redundant_pairs = QCO.get_redundant_gate_product_pairs(gates_dict, qcm.data["decomposition_type"])
    
    if !isempty(redundant_pairs)
        (length(redundant_pairs) == 1) && (Memento.info(_LOGGER, "Detected $(length(redundant_pairs)) redundant input elementary gate pair"))
        (length(redundant_pairs) > 1)  && (Memento.info(_LOGGER, "Detected $(length(redundant_pairs)) redundant input elementary gate pairs"))

        for i = 1:length(redundant_pairs)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], z_bin_var[redundant_pairs[i][1], d] + z_bin_var[redundant_pairs[i][2], d+1] <= 1)
        end
    end

    return
end

function constraint_redundant_gate_product_triplets(qcm::QuantumCircuitModel)

    gates_dict = qcm.data["gates_dict"]
    max_depth  = qcm.data["maximum_depth"]
    z_bin_var  = qcm.variables[:z_bin_var]
    num_qubits = qcm.data["num_qubits"]

    redundant_triplets = QCO.get_redundant_gate_product_triplets(gates_dict, num_qubits, qcm.data["decomposition_type"])
    
    if !isempty(redundant_triplets)
        (length(redundant_triplets) == 1) && (Memento.info(_LOGGER, "Detected $(length(redundant_triplets)) redundant input elementary gate triplets"))
        (length(redundant_triplets) > 1)  && (Memento.info(_LOGGER, "Detected $(length(redundant_triplets)) redundant input elementary gate triplets"))

        for i = 1:length(redundant_triplets)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-2)], z_bin_var[redundant_triplets[i][1], d] + z_bin_var[redundant_triplets[i][2], d+1] + z_bin_var[redundant_triplets[i][3], d+2]<= 2)
        end
    end

    return
end

###################################   OTHER   ###################################

function constraint_cnot_gate_bounds(qcm::QuantumCircuitModel)

    cnot_idx  = qcm.data["cnot_idx"]
    max_depth = qcm.data["maximum_depth"]
    z_bin_var = qcm.variables[:z_bin_var]

    if !isempty(cnot_idx)
        if "cnot_lower_bound" in keys(qcm.data)
            JuMP.@constraint(qcm.model, sum(z_bin_var[n,d] for n in cnot_idx, d=1:max_depth) >= qcm.data["cnot_lower_bound"])
            Memento.info(_LOGGER, "Applied CNot-gate lower bound constraint")
        end
        
        if "cnot_upper_bound" in keys(qcm.data)
            JuMP.@constraint(qcm.model, sum(z_bin_var[n,d] for n in cnot_idx, d=1:max_depth) <= qcm.data["cnot_upper_bound"])
            Memento.info(_LOGGER, "Applied CNot-gate upper bound constraint")
        end
    end

    return
end

function constraint_convex_hull_complex_gates(qcm::QuantumCircuitModel)

    if !qcm.data["are_gates_real"] 

        max_ex_pt  = 4 # (>= 2) A parameter which can be an user input

        z_bin_var  = qcm.variables[:z_bin_var]

        gates_real = qcm.data["gates_real"]
        gates_dict = qcm.data["gates_dict"]
        num_gates  = size(gates_real)[3]
        max_depth  = qcm.data["maximum_depth"]
        n_r        = size(gates_dict["1"]["matrix"])[1]
        n_c        = size(gates_dict["1"]["matrix"])[2]

        num_facets = 0

        vertices_dict = Dict{Set{}, Tuple{<:Number, <:Number}}()
        for I=1:n_r, J=1:n_c
            vertices_coord_IJ = Set()

            for K in keys(gates_dict)
                re = QCO.round_real_value(real(gates_dict[K]["matrix"][I,J]))
                im = QCO.round_real_value(imag(gates_dict[K]["matrix"][I,J]))

                push!(vertices_coord_IJ, (re, im))
            end
            
            if (isapprox(minimum([x[1] for x in vertices_coord_IJ]), maximum([x[1] for x in vertices_coord_IJ]), atol = 1E-6)) || (isapprox(minimum([x[2] for x in vertices_coord_IJ]), maximum([x[2] for x in vertices_coord_IJ]), atol = 1E-6))
                continue
            else 
                # Eliminates repeated set of extreme points
                vertices_dict[vertices_coord_IJ] = (I,J)
            end
        end

        for vertices_coord in keys(vertices_dict)
            I = vertices_dict[vertices_coord][1]
            J = vertices_dict[vertices_coord][2]

            vertices = Vector{Tuple{<:Number, <:Number}}()
            
            if length(vertices_coord) == 2 
                for l in vertices_coord
                    push!(vertices, (l[1], l[2]))
                end

                slope, intercept = QCO._get_constraint_slope_intercept(vertices[1], vertices[2])

                if !isinf(slope)
                    if QCO.is_zero(abs(slope))
                        JuMP.@constraint(qcm.model, [d=1:max_depth], 
                                        sum(gates_real[(2*I-1),(2*J), n_g] * z_bin_var[n_g,d] for n_g = 1:num_gates) - intercept == 0)
                    else
                        
                        JuMP.@constraint(qcm.model, [d=1:max_depth], sum(gates_real[(2*I-1),(2*J), n_g] * z_bin_var[n_g,d] for n_g = 1:num_gates) 
                                                                - slope*sum(gates_real[(2*I-1),(2*J-1), n_g] * z_bin_var[n_g,d] for n_g = 1:num_gates) - intercept == 0)
                    end
                elseif isinf(slope)
                    JuMP.@constraint(qcm.model, [d=1:max_depth], sum(gates_real[(2*I-1),(2*J-1), n_g] * z_bin_var[n_g,d] for n_g = 1:num_gates) == vertices[1][1])
                end
                
                num_facets += 1

            elseif (length(vertices_coord) > 2)
                
                for l in vertices_coord
                    push!(vertices, (l[1], l[2]))
                end
                
                vertices_convex_hull = QCO.convex_hull(vertices)
                num_ex_pt = size(vertices_convex_hull)[1]
                
                # Add a planar hull cut if num_ex_pt == 2

                if (num_ex_pt >= 3) && (num_ex_pt <= max_ex_pt)

                    for i=1:num_ex_pt
                        v1 = vertices_convex_hull[i]

                        if i == num_ex_pt
                            v2 = vertices_convex_hull[1]
                        else 
                            v2 = vertices_convex_hull[i+1]
                        end

                        # Test-vertex for half-space directionality
                        if i == (num_ex_pt - 1)
                            v3 = vertices_convex_hull[1]
                        elseif i == num_ex_pt
                            v3 = vertices_convex_hull[2]
                        else 
                            v3 = vertices_convex_hull[i+2]
                        end

                        slope, intercept = QCO._get_constraint_slope_intercept(v1, v2)
                        
                        # Facets of the hull
                        if !isinf(slope)

                            if v3[2] - slope*v3[1] - intercept <= -1E-6

                                if QCO.is_zero(abs(slope))             
                                    JuMP.@constraint(qcm.model, [d=1:max_depth], sum(gates_real[(2*I-1),(2*J), n_g] * z_bin_var[n_g,d] for n_g = 1:num_gates) - intercept <= 0)
                                else           
                                    JuMP.@constraint(qcm.model, [d=1:max_depth], sum(gates_real[(2*I-1),(2*J), n_g] * z_bin_var[n_g,d] for n_g = 1:num_gates) 
                                                                        - slope*(sum(gates_real[(2*I-1),(2*J-1), n_g] * z_bin_var[n_g,d] for n_g = 1:num_gates)) - intercept <= 0)
                                end
                                num_facets += 1

                            elseif v3[2] - slope*v3[1] - intercept >= 1E-6

                                if QCO.is_zero(abs(slope))
                                    JuMP.@constraint(qcm.model, [d=1:max_depth], sum(gates_real[(2*I-1),(2*J), n_g] * z_bin_var[n_g,d] for n_g = 1:num_gates) - intercept >= 0)
                                else
                                    JuMP.@constraint(qcm.model, [d=1:max_depth], sum(gates_real[(2*I-1),(2*J), n_g] * z_bin_var[n_g,d] for n_g = 1:num_gates) 
                                                                        - slope*(sum(gates_real[(2*I-1),(2*J-1), n_g] * z_bin_var[n_g,d] for n_g = 1:num_gates)) - intercept >= 0)
                                end
                                num_facets += 1
                                
                            else 
                                Memento.warn(_LOGGER, "Indeterminate direction for the planar-hull cut")
                            end

                        else isinf(slope)

                            if v3[1] >= v1[1] + 1E-6
                                JuMP.@constraint(qcm.model, [d=1:max_depth], sum(gates_real[(2*I-1),(2*J-1), n_g] * z_bin_var[n_g,d] for n_g = 1:num_gates) >= v1[1])
                            elseif v3[1] <= v1[1] - 1E-6
                                JuMP.@constraint(qcm.model, [d=1:max_depth], sum(gates_real[(2*I-1),(2*J-1), n_g] * z_bin_var[n_g,d] for n_g = 1:num_gates) <= v1[1])
                            else
                                Memento.warn(_LOGGER, "Indeterminate direction for the planar-hull cut")
                            end
                            num_facets += 1
                                
                        end
                    end
                end
            end
        end
        
        if num_facets > 0
            Memento.info(_LOGGER, "Applied $num_facets planar-hull cuts per max_depth of the decomposition")
        end
        
    end
    
    return
end

function constraint_unitary_property(qcm::QuantumCircuitModel)
    max_depth  = qcm.data["maximum_depth"]
    Z_var      = qcm.variables[:Z_var]
    z_bin_var  = qcm.variables[:z_bin_var]
    gates_real = qcm.data["gates_real"]
    num_gates  = size(gates_real)[3]

    n_r = size(gates_real)[1]
    n_c = size(gates_real)[2]
    gates_adjoint = zeros(n_r, n_c, num_gates)
    
    for k=1:num_gates
        if qcm.data["are_gates_real"]
            gates_adjoint[:,:,k] = gates_real[:,:,k]'
        else
            gate_complex = QCO.real_to_complex_gate(gates_real[:,:,k])
            gate_complex_adj = Matrix{ComplexF64}(adjoint(gate_complex))
            gates_adjoint[:,:,k] = QCO.complex_to_real_gate(gate_complex_adj)
        end
    end

    for d = 1:max_depth
        for i=1:num_gates
            for j=i:num_gates
                if i == j 
                    JuMP.@constraint(qcm.model, Z_var[i,j,d] == z_bin_var[i,d])
                else 
                    QCO.relaxation_bilinear(qcm.model, Z_var[i,j,d], z_bin_var[i,d], z_bin_var[j,d])
                    # JuMP.@constraint(qcm.model, Z_var[i,j,d] == 0)
                    JuMP.@constraint(qcm.model, Z_var[j,i,d] == Z_var[i,j,d])
                end
            end
            # RLT-type constraint
            JuMP.@constraint(qcm.model, sum(Z_var[i,:,d]) == z_bin_var[i,d])
        end
        # JuMP.@constraint(qcm.model, sum(Z_var[i,i,d] for i=1:num_gates) == 1)
    end

    # Unitary constraint
    JuMP.@constraint(qcm.model, [d=1:max_depth], sum(Z_var[i,j,d]* QCO.round_real_value.(gates_real[:,:,i] * gates_adjoint[:,:,j] + 
                                                 gates_real[:,:,j] * gates_adjoint[:,:,i]) for i=1:(num_gates-1), 
                                                 j=(i+1):num_gates) .== zeros(n_r, n_c))

    return
end

function constraint_unitary_complex_conjugate(qcm::QuantumCircuitModel; 
                                              num_gates_bnd = 50)

    max_depth   = qcm.data["maximum_depth"]
    z_bin_var   = qcm.variables[:z_bin_var]
    U_var       = qcm.variables[:U_var]
    gates_real  = qcm.data["gates_real"]
    n_r         = size(gates_real)[1]
    n_c         = size(gates_real)[2]
    num_gates   = size(gates_real)[3]

    decomposition_type = qcm.data["decomposition_type"]
    conj_option        = qcm.options.unitary_complex_conjugate
    target_gate        = qcm.data["target_gate"]

    # Helper for repeated sums of gates
    function sum_of_bin_vars(depth)
        sum(z_bin_var[i, depth] * (gates_real[:, :, i]') for i in 1:num_gates)
    end

    # Decomposition without global phase
    if decomposition_type !== "optimal_global_phase"

        # Linear constraints
        if (max_depth >= 2) && (conj_option in [1, 2])
            Memento.info(_LOGGER, "Applying linear unitary complex-conjugate constraints")
            JuMP.@constraint(qcm.model,
                U_var[:, :, max_depth-1] .== target_gate * sum_of_bin_vars(max_depth)
            )
        end

        # Quadratic constraints
        if (max_depth >= 3) && (num_gates <= num_gates_bnd) && (conj_option in [2])
            Memento.info(_LOGGER, "Applying quadratic unitary complex-conjugate constraints")
            JuMP.@constraint(qcm.model, U_var[:, :, max_depth-2] .== target_gate * 
                        sum_of_bin_vars(max_depth) * sum_of_bin_vars(max_depth-1))
        end
    end

    # Decomposition with global phase
    #=
    if decomposition_type == "optimal_global_phase"

        # Linear constraints
        if (max_depth >= 2) && (conj_option in [1, 2])
            Memento.info(_LOGGER, "Applying linear unitary complex-conjugate constraints (global-phase)")

            for gate in 1:num_gates
                product_gate = target_gate * gates_real[:, :, gate]'

                if qcm.data["are_gates_real"]
                    ref_r, ref_c = QCO._get_nonzero_idx_of_complex_matrix(Array{ComplexF64}(product_gate))

                    for i in 1:n_r, j in 1:n_c
                        if QCO.is_zero(product_gate[i, j])
                            JuMP.@constraint(qcm.model, U_var[i, j, max_depth-1] ≥ -(1 - z_bin_var[gate, max_depth]))
                            JuMP.@constraint(qcm.model, U_var[i, j, max_depth-1] ≤  (1 - z_bin_var[gate, max_depth]))

                        elseif !((i == ref_r) && (j == ref_c))
                            expr_diff = QCO._real_phase_difference(U_var, product_gate, ref_r, ref_c, i, j, max_depth-1)

                            big_M = 2  # assumes that all elementary and target gate elements lie within the range [-1, 1].
                            JuMP.@constraint(qcm.model, expr_diff ≥ -big_M * (1 - z_bin_var[gate, max_depth]))
                            JuMP.@constraint(qcm.model, expr_diff ≤  big_M * (1 - z_bin_var[gate, max_depth]))
                        end
                    end

                else
                    ref_r, ref_c = QCO._get_nonzero_idx_of_complex_to_real_matrix(product_gate)
                    for i in 1:2:n_r, j in 1:2:n_c
                        if QCO.is_zero(product_gate[i, j]) && QCO.is_zero(product_gate[i, j+1])
                            JuMP.@constraint(qcm.model, U_var[i, j, max_depth-1] ≥ -(1 - z_bin_var[gate, max_depth]))
                            JuMP.@constraint(qcm.model, U_var[i, j+1, max_depth-1] ≤  (1 - z_bin_var[gate, max_depth]))

                        elseif !((i == ref_r) && (j == ref_c))
                            realLHS, realRHS, imagLHS, imagRHS = QCO._complex_phase_differences(U_var, product_gate, ref_r, ref_c, i, j, max_depth-1)
                            
                            big_M = 2  # assumes that all elementary and target gate elements lie within the range [-1, 1].
                            JuMP.@constraint(qcm.model, realLHS - realRHS ≥ -big_M * (1 - z_bin_var[gate, max_depth]))
                            JuMP.@constraint(qcm.model, realLHS - realRHS ≤  big_M * (1 - z_bin_var[gate, max_depth]))
                            JuMP.@constraint(qcm.model, imagLHS - imagRHS ≥ -big_M * (1 - z_bin_var[gate, max_depth]))
                            JuMP.@constraint(qcm.model, imagLHS - imagRHS ≤  big_M * (1 - z_bin_var[gate, max_depth]))
                        end
                    end
                end
            end
        end
    end
    =#
    if decomposition_type == "optimal_global_phase"

        phase_real_var     = qcm.variables[:phase_real_var]
        phase_complex_var  = qcm.variables[:phase_complex_var]

        # Linear constraints
        if (max_depth >= 2) && (conj_option in [1, 2])
            Memento.info(_LOGGER, "Applying linear unitary complex-conjugate constraints (global-phase)")

            # For each candidate gate g at pozition P
            for gate in 1:num_gates
                # Real-embedding: R(T * g^†) = R(T) * R(g)^T
                product_gate = target_gate * gates_real[:, :, gate]'
                zc = 2 * (1 - z_bin_var[gate, max_depth])  # 0 when this gate is chosen

                if qcm.data["are_gates_real"]
                    
                    for i in 1:n_r, j in 1:n_c
                        JuMP.@constraint(qcm.model,  U_var[i,j,max_depth-1] - product_gate[i,j] * phase_real_var ≤ zc) 
                        JuMP.@constraint(qcm.model,  - U_var[i,j,max_depth-1] + product_gate[i,j] * phase_real_var ≤ zc)                        
                    end

                else

                    for i in 1:2:n_r, j in 1:2:n_c
                        # Extract A_g and B_g from the real-embedded product:
                        A = product_gate[i,   j  ]   # Re
                        B = product_gate[i+1, j  ]   # Im  (since bottom-left block stores +Im)

                        # Re(G) = r*A - s*B ; Im(G) = r*B + s*A
                        t_re = phase_real_var * A - phase_complex_var * B
                        t_im = phase_real_var * B + phase_complex_var * A

                        # --- Re constraints on (i,j) and (i+1,j+1) ---
                        JuMP.@constraint(qcm.model,  U_var[i,   j,   max_depth-1] - t_re <= zc)
                        JuMP.@constraint(qcm.model, -U_var[i,   j,   max_depth-1] + t_re <= zc)

                        JuMP.@constraint(qcm.model,  U_var[i+1, j+1, max_depth-1] - t_re <= zc)
                        JuMP.@constraint(qcm.model, -U_var[i+1, j+1, max_depth-1] + t_re <= zc)

                        # --- Im constraints ---
                        JuMP.@constraint(qcm.model,  U_var[i+1, j,   max_depth-1] - t_im <= zc)
                        JuMP.@constraint(qcm.model, -U_var[i+1, j,   max_depth-1] + t_im <= zc)

                        JuMP.@constraint(qcm.model,  U_var[i,   j+1, max_depth-1] + t_im <= zc)
                        JuMP.@constraint(qcm.model, -U_var[i,   j+1, max_depth-1] - t_im <= zc)
                    end    
                end
            end
        end
    end




    return
end

function constraint_unitary_complex_conjugate_at_depth(qcm::QuantumCircuitModel, at_depth::Int64)

    z_bin_var   = qcm.variables[:z_bin_var]
    U_var       = qcm.variables[:U_var]
    gates_real  = qcm.data["gates_real"]
    num_gates   = size(gates_real)[3]

    decomposition_type = qcm.data["decomposition_type"]

    # Helper for repeated sums of gates
    function sum_of_bin_vars(depth)
        sum(z_bin_var[i, depth] * (gates_real[:, :, i]') for i in 1:num_gates)
    end

    Memento.info(_LOGGER, "Applying bilinear unitary complex-conjugate constraints for depth $at_depth")
    JuMP.@constraint(qcm.model, U_var[:, :, at_depth-1] .== U_var[:, :, at_depth] * sum_of_bin_vars(at_depth))
    
end

function constraint_unitarity_of_U(qcm::QuantumCircuitModel,at_depth::Int64)

    n_r            = size(qcm.data["gates_real"])[1]
    n_c            = size(qcm.data["gates_real"])[2]
    target_gate    = qcm.data["target_gate"]
    U_var          = qcm.variables[:U_var]

    if qcm.data["are_gates_real"]
        
        for i in 1:n_r, j in 1:n_c
            if !(i==j)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,at_depth] * U_var[j,k,at_depth] for k=1:n_r) == 0)
            end
        end

        for i in 2:n_r
                JuMP.@constraint(qcm.model, sum( U_var[i,k,at_depth] * U_var[i,k,at_depth] for k=1:n_r) == 1)
        end

    else
        for i in 1:2:n_r, j in 1:2:n_c
            if !(i==j)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,at_depth] * U_var[j,k,at_depth] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,at_depth] * U_var[j+1,k,at_depth] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,at_depth] * U_var[j,k,at_depth] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,at_depth] * U_var[j+1,k,at_depth] for k=1:n_r) == 0)
            end
        end

        for i in 1:2:n_r
                JuMP.@constraint(qcm.model, sum( U_var[i,k,at_depth] * U_var[i,k,at_depth] for k=1:n_r) == 1)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,at_depth] * U_var[i,k,at_depth] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,at_depth] * U_var[i+1,k,at_depth] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,at_depth] * U_var[i+1,k,at_depth] for k=1:n_r) == 1)
        end
            
    end

    return
end

function constraint_unitarity_of_U_diag(qcm::QuantumCircuitModel,at_depth::Int64)

    n_r            = size(qcm.data["gates_real"])[1]
    n_c            = size(qcm.data["gates_real"])[2]
    target_gate    = qcm.data["target_gate"]
    U_var          = qcm.variables[:U_var]

    if qcm.data["are_gates_real"]
        
        for i in 2:n_r
                JuMP.@constraint(qcm.model, sum( U_var[i,k,at_depth] * U_var[i,k,at_depth] for k=1:n_r) == 1)
        end

    else
        
        for i in 1:2:n_r
                JuMP.@constraint(qcm.model, sum( U_var[i,k,at_depth] * U_var[i,k,at_depth] for k=1:n_r) == 1)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,at_depth] * U_var[i,k,at_depth] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i,k,at_depth] * U_var[i+1,k,at_depth] for k=1:n_r) == 0)
                JuMP.@constraint(qcm.model, sum( U_var[i+1,k,at_depth] * U_var[i+1,k,at_depth] for k=1:n_r) == 1)
        end
            
    end

    return
end

function constraint_complex_to_real_symmetry(qcm::QuantumCircuitModel)

    max_depth  = qcm.data["maximum_depth"]
    n_r    = size(qcm.data["gates_real"])[1]
    n_c    = size(qcm.data["gates_real"])[2]
    
    for i=1:2:n_r, j=1:2:n_c, d=1:(max_depth-1)
        JuMP.@constraint(qcm.model, qcm.variables[:U_var][i,j,d] == qcm.variables[:U_var][i+1,j+1,d])
        JuMP.@constraint(qcm.model, qcm.variables[:U_var][i,j+1,d] == -qcm.variables[:U_var][i+1,j,d]) # This seems to slow down the solution search
    end

    return
end