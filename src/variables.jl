#------------------------------------------------#
# Initialize all variables of the QC_model here  #
#------------------------------------------------#

function variable_matrix_per_depth(qcm::QuantumCircuitModel)
    tol_0 = 1E-6
    n_r     = size(qcm.data["M_real"])[1]
    n_c     = size(qcm.data["M_real"])[2]
    n_gates = size(qcm.data["M_real"])[3]
    depth   = qcm.data["depth"]

    if n_r != n_c
        Memento.warn(_LOGGER, "number of rows and columns have to be equal for unitary quantum gates")
    end

    M_real_l, M_real_u = QCO.get_gate_element_bounds(qcm.data["M_real"])

    qcm.variables[:M_var] = JuMP.@variable(qcm.model, M_real_l[i,j] <= M_var[i=1:n_r, j=1:n_c, 1:depth] <= M_real_u[i,j])

    num_vars_fixed = 0
    for i=1:n_r
        for j=1:n_c
            if (abs(M_real_l[i,j] - M_real_u[i,j])) <= tol_0
                for d=1:depth
                    JuMP.fix(M_var[i,j,d], M_real_l[i,j]; force=true)
                end
                num_vars_fixed += 1
            end
        end
    end
    if num_vars_fixed > 1 
        Memento.info(_LOGGER, "$num_vars_fixed of $(n_r * n_c) number of entries in the given set of gates (in real form) have equal lower and upper bounds.")
    end
    return
end

function variable_gates_onoff(qcm::QuantumCircuitModel)
    n_gates = size(qcm.data["M_real"])[3]
    depth   = qcm.data["depth"]

    qcm.variables[:z_onoff_var] = JuMP.@variable(qcm.model, z_onoff[1:n_gates,1:depth], Bin)
    return
end

function variable_sequential_gate_products(qcm::QuantumCircuitModel)
    depth  = qcm.data["depth"]
    n_r    = size(qcm.data["M_real"])[1]
    n_c    = size(qcm.data["M_real"])[2]

    qcm.variables[:U_var] = JuMP.@variable(qcm.model, -1 <= U_var[1:n_r, 1:n_c, 1:(depth-1)] <= 1)
    return
end

function variable_gate_products_copy(qcm::QuantumCircuitModel)
    depth   = qcm.data["depth"]
    n_r     = size(qcm.data["M_real"])[1]
    n_c     = size(qcm.data["M_real"])[2]
    n_gates = size(qcm.data["M_real"])[3]

    qcm.variables[:V_var] = JuMP.@variable(qcm.model, -1 <= V_var[1:n_r, 1:n_c, 1:n_gates, 1:depth] <= 1)
    return
end

function variable_gate_products_linearization(qcm::QuantumCircuitModel)
    depth   = qcm.data["depth"]
    n_r     = size(qcm.data["M_real"])[1]
    n_c     = size(qcm.data["M_real"])[2]
    n_gates = size(qcm.data["M_real"])[3]

    qcm.variables[:zU_var] = JuMP.@variable(qcm.model, -1 <= zU_var[1:n_r, 1:n_c, 1:n_gates, 1:(depth-1)] <= 1)
    return
end