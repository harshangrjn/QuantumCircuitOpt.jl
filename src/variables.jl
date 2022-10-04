#-----------------------------------------------------------#
# Initialize all variables of the QuantumCircuitModel here  #
#-----------------------------------------------------------#

function variable_gates_per_depth(qcm::QuantumCircuitModel)
    
    tol_0 = 1E-6
    n_r     = size(qcm.data["gates_real"])[1]
    n_c     = size(qcm.data["gates_real"])[2]
    max_depth   = qcm.data["maximum_depth"]

    M_l, M_u = QCO.gate_element_bounds(qcm.data["gates_real"])

    qcm.variables[:G_var] = JuMP.@variable(qcm.model, M_l[i,j] <= G_var[i=1:n_r, j=1:n_c, 1:max_depth] <= M_u[i,j])

    num_vars_fixed = 0

    for i=1:n_r, j=1:n_c
        if isapprox(M_l[i,j], M_u[i,j], atol=tol_0)
            for d=1:max_depth
                JuMP.fix(G_var[i,j,d], M_l[i,j]; force=true)
            end
            num_vars_fixed += 1
        end
    end

    if num_vars_fixed > 0
        Memento.info(_LOGGER, "$num_vars_fixed of $(n_r * n_c) entries are constants in the given set of gates.")
    end
    
    return
end

function variable_gates_onoff(qcm::QuantumCircuitModel)
    num_gates = size(qcm.data["gates_real"])[3]
    max_depth   = qcm.data["maximum_depth"]

    qcm.variables[:z_bin_var] = JuMP.@variable(qcm.model, z_bin_var[1:num_gates,1:max_depth], Bin)

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
    max_depth  = qcm.data["maximum_depth"]
    n_r    = size(qcm.data["gates_real"])[1]
    n_c    = size(qcm.data["gates_real"])[2]

    if !qcm.options.fix_unitary_variables
        qcm.variables[:U_var] = JuMP.@variable(qcm.model, -1 <= U_var[1:n_r, 1:n_c, 1:max_depth] <= 1)
        return 
    end

    qcm.variables[:U_var] = JuMP.@variable(qcm.model, U_var[1:n_r, 1:n_c, 1:max_depth])
    U_fixed_idx = QCO._get_unitary_variables_fixed_indices(qcm.data["gates_real"], max_depth)
    # U_var_bound_tol = 1E-8

    for depth = 1:(max_depth-1)
        U_fix = U_fixed_idx[depth]
        for ii = 1:n_r, jj = 1:n_c
            if ((ii, jj) in keys(U_fix))
                if isapprox(U_fix[(ii, jj)]["value"], -1, atol = 1E-6)
                    JuMP.set_lower_bound(U_var[ii, jj, depth], -1)
                    JuMP.set_upper_bound(U_var[ii, jj, depth], -1)
                elseif isapprox(U_fix[(ii, jj)]["value"], 0, atol = 1E-6)
                    JuMP.set_lower_bound(U_var[ii, jj, depth], 0)
                    JuMP.set_upper_bound(U_var[ii, jj, depth], 0)
                elseif isapprox(U_fix[(ii, jj)]["value"], 1, atol = 1E-6)
                    JuMP.set_lower_bound(U_var[ii, jj, depth], 1)
                    JuMP.set_upper_bound(U_var[ii, jj, depth], 1)
                else 
                    JuMP.set_lower_bound(U_var[ii, jj, depth], U_fix[(ii, jj)]["value"])
                    JuMP.set_upper_bound(U_var[ii, jj, depth], U_fix[(ii, jj)]["value"])
                end
            else
                JuMP.set_lower_bound(U_var[ii, jj, depth], -1)
                JuMP.set_upper_bound(U_var[ii, jj, depth],  1)
            end
        end
    end

    for ii = 1:n_r, jj = 1:n_c
        JuMP.set_lower_bound(U_var[ii, jj, max_depth], -1)
        JuMP.set_upper_bound(U_var[ii, jj, max_depth],  1)
    end
    
    return
end

function variable_gate_products_copy(qcm::QuantumCircuitModel)
    max_depth   = qcm.data["maximum_depth"]
    n_r     = size(qcm.data["gates_real"])[1]
    n_c     = size(qcm.data["gates_real"])[2]
    num_gates = size(qcm.data["gates_real"])[3]

    qcm.variables[:V_var] = JuMP.@variable(qcm.model, -1 <= V_var[1:n_r, 1:n_c, 1:num_gates, 1:max_depth] <= 1)
    
    return
end

function variable_gate_products_linearization(qcm::QuantumCircuitModel)
    n_r     = size(qcm.data["gates_real"])[1]
    n_c     = size(qcm.data["gates_real"])[2]
    max_depth   = qcm.data["maximum_depth"]
    num_gates = size(qcm.data["gates_real"])[3]

    if !qcm.options.fix_unitary_variables
        qcm.variables[:zU_var] = JuMP.@variable(qcm.model, -1 <= zU_var[1:n_r, 1:n_c, 1:num_gates, 1:(max_depth-1)] <= 1)
        return 
    end

    U_var = qcm.variables[:U_var]
    qcm.variables[:zU_var] = JuMP.@variable(qcm.model, zU_var[1:n_r, 1:n_c, 1:num_gates, 1:(max_depth-1)])

    for d = 1:(max_depth-1), i = 1:n_r, j = 1:n_c, k = 1:num_gates
        U_var_l = JuMP.lower_bound(U_var[i,j,d])
        U_var_u = JuMP.upper_bound(U_var[i,j,d])

        if isapprox(U_var_l, 0, atol = 1E-6) && isapprox(U_var_u, 0, atol = 1E-6)
            JuMP.set_lower_bound(zU_var[i,j,k,d], 0)
            JuMP.set_upper_bound(zU_var[i,j,k,d], 0)
        else
            JuMP.set_lower_bound(zU_var[i,j,k,d], -1)
            JuMP.set_upper_bound(zU_var[i,j,k,d],  1)
        end
    end 
    
    return
end

function variable_slack_for_feasibility(qcm::QuantumCircuitModel)
    n_r     = size(qcm.data["gates_real"])[1]
    n_c     = size(qcm.data["gates_real"])[2]
    max_depth   = qcm.data["maximum_depth"]
    U_var = qcm.variables[:U_var]
    
    qcm.variables[:slack_var] = JuMP.@variable(qcm.model, slack_var[1:n_r, 1:n_c])

    for i=1:n_r, j=1:n_c
        lb = JuMP.lower_bound(U_var[i,j,max_depth-1])
        ub = JuMP.upper_bound(U_var[i,j,max_depth-1])
        if isapprox(lb, 0, atol = 1E-6) && isapprox(ub, 0, atol = 1E-6)
            JuMP.set_lower_bound(slack_var[i,j], 0)
            JuMP.set_upper_bound(slack_var[i,j], 0)
        else 
            JuMP.set_lower_bound(slack_var[i,j], -1)
            JuMP.set_upper_bound(slack_var[i,j],  1)
        end
    end
    return
end

function variable_slack_var_outer_approximation(qcm::QuantumCircuitModel)
    n_r     = size(qcm.data["gates_real"])[1]
    n_c     = size(qcm.data["gates_real"])[2]

    qcm.variables[:slack_var_oa] = JuMP.@variable(qcm.model, -1 <= slack_var_oa[1:n_r, 1:n_c] <= 1)
    return
end

function variable_binary_products(qcm::QuantumCircuitModel)
    num_gates = size(qcm.data["gates_real"])[3]
    max_depth = qcm.data["maximum_depth"]
    
    qcm.variables[:Z_var] = JuMP.@variable(qcm.model, 0 <= Z[1:num_gates, 1:num_gates, 1:max_depth] <= 1)
    return
end