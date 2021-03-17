#---------------------------------------------#
# Build all constraints of the QC_model here  #
#---------------------------------------------#
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

function constraint_gate_product_linearization(qcm::QuantumCircuitModel)
    depth   = qcm.data["depth"]
    n_r     = size(qcm.data["M_real"])[1]
    n_c     = size(qcm.data["M_real"])[2]
    n_gates = size(qcm.data["M_real"])[3]

    for i=1:2:n_r
        for j=1:2:n_c
            for n=1:n_gates
                for d=1:(depth-1)
                    QCO.get_mccormick_relaxation(qcm.model, qcm.variables[:zU_var][i,j,n,d], qcm.variables[:U_var][i,j,d], qcm.variables[:z_onoff_var][n,(d+1)])
                    JuMP.@constraint(qcm.model, qcm.variables[:zU_var][i+1,j+1,n,d] ==  qcm.variables[:zU_var][i,j,n,d])
                    JuMP.@constraint(qcm.model, qcm.variables[:zU_var][i,j+1,n,d]   == -qcm.variables[:zU_var][i+1,j,n,d])
                end
            end
        end
    end
    
    return
end

function constraint_gate_target_condition(qcm::QuantumCircuitModel)
    depth   = qcm.data["depth"]
    n_gates = size(qcm.data["M_real"])[3]
    
    # For correct implementation of this, use MutableArithmetics.jl >= v0.2.11
    JuMP.@constraint(qcm.model, sum(qcm.variables[:V_var][:,:,n,depth] * qcm.data["M_real"][:,:,n] for n=1:n_gates) .== qcm.data["Target_real"])  
    
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