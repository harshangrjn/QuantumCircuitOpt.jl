#--------------------------------------------------------#
# Build all constraints of the QuantumCircuitModel here  #
#--------------------------------------------------------#

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

function constraint_target_gate_condition(qcm::QuantumCircuitModel)

    max_depth          = qcm.data["maximum_depth"]
    num_gates          = size(qcm.data["gates_real"])[3]
    decomposition_type = qcm.data["decomposition_type"]
    
    if decomposition_type in ["exact_optimal", "exact_feasible"]
        JuMP.@constraint(qcm.model, sum(qcm.variables[:V_var][:,:,n,max_depth] * qcm.data["gates_real"][:,:,n] for n=1:num_gates) .== qcm.data["target_gate"])
    
    elseif decomposition_type == "approximate"
        JuMP.@constraint(qcm.model, sum(qcm.variables[:V_var][:,:,n,max_depth] * qcm.data["gates_real"][:,:,n] for n=1:num_gates) .== qcm.data["target_gate"][:,:] + qcm.variables[:slack_var][:,:])  

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

function constraint_target_gate_condition_compact(qcm::QuantumCircuitModel)

    max_depth = qcm.data["maximum_depth"]
    decomposition_type = qcm.data["decomposition_type"]

    U_var = qcm.variables[:U_var]
    
    if decomposition_type in ["exact_optimal", "exact_feasible"]
        JuMP.@constraint(qcm.model, U_var[:,:,max_depth] .== qcm.data["target_gate"][:,:])
    
    elseif decomposition_type == "approximate"
        JuMP.@constraint(qcm.model, U_var[:,:,max_depth] .== qcm.data["target_gate"][:,:] + qcm.variables[:slack_var][:,:])
    end
    
    return
end

function constraint_target_gate_condition_glphase(qcm::QuantumCircuitModel)
    
    max_depth      = qcm.data["maximum_depth"]
    n_r            = size(qcm.data["gates_real"])[1]
    n_c            = size(qcm.data["gates_real"])[2]
    target_gate    = qcm.data["target_gate"]

    U_var = qcm.variables[:U_var]    

    if qcm.data["are_gates_real"] 
        ref_r, ref_c = QCO._get_nonzero_idx_of_complex_matrix(Array{Complex{Float64},2}(target_gate))
    
        for i=1:n_r, j=1:n_c
            if QCO.is_zero(target_gate[i,j])
                JuMP.@constraint(qcm.model, U_var[i,j,max_depth] == 0)
            else
                if !((i == ref_r) && (j == ref_c))
                    JuMP.@constraint(qcm.model, U_var[i,j,max_depth] * target_gate[ref_r,ref_c] == 
                                                U_var[ref_r,ref_c,max_depth] * target_gate[i,j])
                end
            end
        end
    
    else
        ref_r, ref_c = QCO._get_nonzero_idx_of_complex_to_real_matrix(target_gate)
    
        for i=1:2:n_r, j=i:2:n_c
            isapprox(target_gate[i,j], target_gate[j,i], atol = 1E-6) 
            if QCO.is_zero(target_gate[i,j]) && QCO.is_zero(target_gate[i,j+1])
                JuMP.@constraint(qcm.model, U_var[i,j,max_depth] == 0.0)
                JuMP.@constraint(qcm.model, U_var[i,(j+1),max_depth] == 0.0)
                
            else
                if !((i == ref_r) && (j == ref_c))
                    
                    #real parts equal
                    JuMP.@constraint(qcm.model, U_var[i,j,max_depth]*target_gate[ref_r,ref_c] - 
                                                U_var[i,(j+1),max_depth]*target_gate[ref_r,(ref_c+1)] ==
                                                U_var[ref_r,ref_c,max_depth]*target_gate[i,j] - 
                                                U_var[ref_r,(ref_c+1),max_depth]*target_gate[i,(j+1)])       
                    #complex parts equal
                    JuMP.@constraint(qcm.model, U_var[i,j,max_depth]*target_gate[ref_r,(ref_c+1)] + 
                                                U_var[i,(j+1),max_depth]*target_gate[ref_r,ref_c] ==
                                                U_var[ref_r,ref_c,max_depth]*target_gate[i,(j+1)] + 
                                                U_var[ref_r,(ref_c+1),max_depth]*target_gate[i,j]) 
                end
            end
        end
    end

    return    
end

function constraint_slack_var_outer_approximation(qcm::QuantumCircuitModel)
    slack_var    = qcm.variables[:slack_var]
    slack_var_oa = qcm.variables[:slack_var_oa]
    
    # Number of under-estimators for the quadratic function
    num_points = 9 

    for i=1:size(slack_var)[1], j=1:size(slack_var)[2]
        lb = JuMP.lower_bound(slack_var[i,j])
        ub = JuMP.upper_bound(slack_var[i,j])
        if !(isapprox(lb, ub, atol = 1E-6))
            oa_points = range(lb, ub, num_points)
            JuMP.@constraint(qcm.model, [k=1:length(oa_points)], 
                            slack_var_oa[i,j] >= 2*slack_var[i,j]*oa_points[k] - (oa_points[k])^2)
            if !(0 in oa_points)
                JuMP.@constraint(qcm.model, slack_var_oa[i,j] >= 0)
            end
        else
            mid_point = (lb+ub)/2
            JuMP.@constraint(qcm.model, slack_var_oa[i,j] >= 2*slack_var[i,j]*mid_point - (mid_point)^2)
        end
    end

    return
end

function constraint_commutative_gate_pairs(qcm::QuantumCircuitModel)
    
    max_depth = qcm.data["maximum_depth"]
    z_bin_var = qcm.variables[:z_bin_var]

    commute_pairs, commute_pairs_prodIdentity = QCO.get_commutative_gate_pairs(qcm.data["gates_dict"], qcm.data["decomposition_type"])

    if !isempty(commute_pairs)
        (length(commute_pairs) == 1) && (Memento.info(_LOGGER, "Detected $(length(commute_pairs)) input elementary gate pair which commutes"))
        (length(commute_pairs) > 1)  && (Memento.info(_LOGGER, "Detected $(length(commute_pairs)) input elementary gate pairs which commute"))

        for i = 1:length(commute_pairs)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], z_bin_var[commute_pairs[i][2], d] + z_bin_var[commute_pairs[i][1], d+1] <= 1)
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

function constraint_involutory_gates(qcm::QuantumCircuitModel)

    gates_dict = qcm.data["gates_dict"]
    max_depth  = qcm.data["maximum_depth"]
    z_bin_var  = qcm.variables[:z_bin_var]

    involutory_gates = QCO.get_involutory_gates(gates_dict)
    
    if !isempty(involutory_gates)
        (length(involutory_gates) == 1) && (Memento.info(_LOGGER, "Detected $(length(involutory_gates)) involutory elementary gate"))
        (length(involutory_gates) > 1)  && (Memento.info(_LOGGER, "Detected $(length(involutory_gates)) involutory elementary gates"))

        for i = 1:length(involutory_gates)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], z_bin_var[involutory_gates[i], d] + z_bin_var[involutory_gates[i], d+1] <= 1)
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
            JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], z_bin_var[redundant_pairs[i][2], d] + z_bin_var[redundant_pairs[i][1], d+1] <= 1)
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

        for i = 1:length(idempotent_gates)
            JuMP.@constraint(qcm.model, [d=1:(max_depth-1)], z_bin_var[idempotent_gates[i], d] + z_bin_var[idempotent_gates[i], d+1] <= 1)
        end
    end

    return
end

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
                                              quadratic_constraint = true,
                                              num_gates_bnd = 50)

    max_depth  = qcm.data["maximum_depth"]
    z_bin_var  = qcm.variables[:z_bin_var]
    U_var  = qcm.variables[:U_var]
    gates_real = qcm.data["gates_real"]

    num_gates  = size(gates_real)[3]

    if max_depth >= 2
        JuMP.@constraint(qcm.model, U_var[:, :, (max_depth-1)] .== 
                                    qcm.data["target_gate"] *
                                    sum(z_bin_var[i, (max_depth)] *  
                                    (gates_real[:,:,i]') for i=1:num_gates)
                        )
    end

    # Quadratic constraint
    if (max_depth >= 3) && (num_gates <= num_gates_bnd) && (quadratic_constraint)
        Memento.info(_LOGGER, "Applying quadratic unitary complex-congugate constraints")
        
        JuMP.@constraint(qcm.model, U_var[:, :, (max_depth-2)] .== 
                                    qcm.data["target_gate"] *
                                    sum(z_bin_var[i, (max_depth)] *  
                                    (gates_real[:,:,i]') for i=1:num_gates) * 
                                    sum(z_bin_var[i, (max_depth-1)] *  
                                    (gates_real[:,:,i]') for i=1:num_gates)
    )

    end

    return
end