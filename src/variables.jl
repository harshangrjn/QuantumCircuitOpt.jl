#-----------------------------------------------------------#
# Initialize all variables of the QuantumCircuitModel here  #
#-----------------------------------------------------------#

function variable_gates_per_depth(qcm::QuantumCircuitModel)
    
    tol_0 = 1E-6
    n_r     = size(qcm.data["gates_real"])[1]
    n_c     = size(qcm.data["gates_real"])[2]
    depth   = qcm.data["maximum_depth"]

    M_l, M_u = QCO.gate_element_bounds(qcm.data["gates_real"])

    qcm.variables[:G_var] = JuMP.@variable(qcm.model, M_l[i,j] <= G_var[i=1:n_r, j=1:n_c, 1:depth] <= M_u[i,j])

    num_vars_fixed = 0

    for i=1:n_r
        for j=1:n_c

            if isapprox(M_l[i,j], M_u[i,j], atol=tol_0)
                for d=1:depth
                    JuMP.fix(G_var[i,j,d], M_l[i,j]; force=true)
                end
                num_vars_fixed += 1
            end

        end
    end

    if num_vars_fixed > 0
        Memento.info(_LOGGER, "$num_vars_fixed of $(n_r * n_c) entries are constants in the given set of gates.")
    end
    
    return
end

function variable_gates_onoff(qcm::QuantumCircuitModel)
    num_gates = size(qcm.data["gates_real"])[3]
    depth   = qcm.data["maximum_depth"]

    qcm.variables[:z_bin_var] = JuMP.@variable(qcm.model, z_bin_var[1:num_gates,1:depth], Bin)

    if "input_circuit" in keys(qcm.data)
        
        gates_dict  = qcm.data["gates_dict"]
        start_value = qcm.data["input_circuit"]
        
        for d in keys(start_value)
            circuit_depth = start_value[d]["depth"]
            gate_start = start_value[d]["gate"]
            
            for n in keys(gates_dict)    
                if gate_start in gates_dict[n]["type"]                 
                    JuMP.set_start_value(z_bin_var[parse(Int64, n), circuit_depth], 1)
                end
            end
        end
    end
        
    return
end

function variable_sequential_gate_products(qcm::QuantumCircuitModel)
    depth  = qcm.data["maximum_depth"]
    n_r    = size(qcm.data["gates_real"])[1]
    n_c    = size(qcm.data["gates_real"])[2]

    qcm.variables[:U_var] = JuMP.@variable(qcm.model, -1 <= U_var[1:n_r, 1:n_c, 1:(depth-1)] <= 1)
    
    return
end

function variable_gate_products_copy(qcm::QuantumCircuitModel)
    depth   = qcm.data["maximum_depth"]
    n_r     = size(qcm.data["gates_real"])[1]
    n_c     = size(qcm.data["gates_real"])[2]
    num_gates = size(qcm.data["gates_real"])[3]

    qcm.variables[:V_var] = JuMP.@variable(qcm.model, -1 <= V_var[1:n_r, 1:n_c, 1:num_gates, 1:depth] <= 1)
    
    return
end

function variable_gate_products_linearization(qcm::QuantumCircuitModel)
    n_r     = size(qcm.data["gates_real"])[1]
    n_c     = size(qcm.data["gates_real"])[2]
    depth   = qcm.data["maximum_depth"]
    num_gates = size(qcm.data["gates_real"])[3]

    qcm.variables[:zU_var] = JuMP.@variable(qcm.model, -1 <= zU_var[1:n_r, 1:n_c, 1:num_gates, 1:(depth-1)] <= 1)
    
    return
end

function variable_slack_for_feasibility(qcm::QuantumCircuitModel)
    n_r     = size(qcm.data["gates_real"])[1]
    n_c     = size(qcm.data["gates_real"])[2]

    qcm.variables[:slack_var] = JuMP.@variable(qcm.model, -1 <= slack_var[1:n_r, 1:n_c] <= 1)
    
    return
end