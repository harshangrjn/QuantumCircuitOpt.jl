#--------------------------------------------------------#
# Build all constraints of the QuantumCircuitModel here  #
#--------------------------------------------------------#

# import QHull: chull 

function constraint_single_gate_per_depth(qcm::QuantumCircuitModel)

    num_gates = size(qcm.data["gates_real"])[3]
    depth     = qcm.data["depth"]
    
    JuMP.@constraint(qcm.model, [d=1:depth], sum(qcm.variables[:z_onoff_var][n,d] for n=1:num_gates) == 1)
    
    return
end

function constraint_gates_onoff_per_depth(qcm::QuantumCircuitModel)

    num_gates = size(qcm.data["gates_real"])[3]
    depth     = qcm.data["depth"]

    JuMP.@constraint(qcm.model, [d=1:depth], qcm.variables[:G_var][:,:,d] .== 
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
    are_gates_real = qcm.data["are_gates_real"]

    i_val = 2 - are_gates_real

    for i=1:i_val:n_r
        for j=1:n_c
            for n=1:num_gates
                for d=1:(depth-1)
                    
                    QCO.relaxation_bilinear(qcm.model, qcm.variables[:zU_var][i,j,n,d], qcm.variables[:U_var][i,j,d], qcm.variables[:z_onoff_var][n,(d+1)])
                    if !are_gates_real
                        if isodd(j)
                            JuMP.@constraint(qcm.model, qcm.variables[:zU_var][i,j,n,d]   ==  qcm.variables[:zU_var][i+1,j+1,n,d])
                            JuMP.@constraint(qcm.model, qcm.variables[:zU_var][i,j+1,n,d] == -qcm.variables[:zU_var][i+1,j,n,d])
                        end
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

function constraint_identity_gate_symmetry(qcm::QuantumCircuitModel)

    gates_dict  = qcm.data["gates_dict"]
    depth       = qcm.data["depth"]
    z_onoff_var = qcm.variables[:z_onoff_var]

    identity_idx = []
    for i=1:length(keys(gates_dict))
        if "Identity" in gates_dict["$i"]["type"]
            push!(identity_idx, i)
        end
    end
    
    if !isempty(identity_idx)
        for i = 1:length(identity_idx)
            JuMP.@constraint(qcm.model, [d=1:(depth-1)], z_onoff_var[identity_idx[i], d] <= z_onoff_var[identity_idx[i], d+1])
        end
    end

    return
end

function constraint_convex_hull_complex_gates(qcm::QuantumCircuitModel)

#     if !qcm.data["are_gates_real"] 

#         max_ex_pt = 10 # (>= 2) A parameter which can be an user input

#         z_onoff_var = qcm.variables[:z_onoff_var]

#         gates_real = qcm.data["gates_real"]
#         gates_dict = qcm.data["gates_dict"]

#         num_gates = size(gates_real)[3]
#         depth     = qcm.data["depth"]
#         n_r       = size(gates_dict["1"]["matrix"])[1]
#         n_c       = size(gates_dict["1"]["matrix"])[2]

#         num_facets = 0

#         for I=1:n_r
#             for J=1:n_c
                
#                 vertices_coord = Set()

#                 for K in keys(gates_dict)

#                     re = QCO.round_real_value(real(gates_dict[K]["matrix"][I,J]))
#                     im = QCO.round_real_value(imag(gates_dict[K]["matrix"][I,J]))

#                     push!(vertices_coord, (re, im))
#                 end

#                 if (isapprox(minimum([x[1] for x in vertices_coord]), maximum([x[1] for x in vertices_coord]), atol = 1E-6)) || (isapprox(minimum([x[2] for x in vertices_coord]), maximum([x[2] for x in vertices_coord]), atol = 1E-6))
#                     continue
#                 end

#                 if length(vertices_coord) == 2 
                    
#                     vertices = Vector{Vector{Float64}}()

#                     for l in vertices_coord
#                         push!(vertices, [l[1], l[2]])
#                     end

#                     m,c = QCO._get_constraint_slope_intercept(vertices[1], vertices[2])
                    
#                     if isapprox(abs(m), 0, atol=1E-6)

#                         JuMP.@constraint(qcm.model, [d=1:depth], 
#                                          sum(gates_real[(2*I-1),(2*J), n_g] * z_onoff_var[n_g,d] for n_g = 1:num_gates) - c == 0)
#                     else
                        
#                         JuMP.@constraint(qcm.model, [d=1:depth], sum(gates_real[(2*I-1),(2*J), n_g] * z_onoff_var[n_g,d] for n_g = 1:num_gates) 
#                                                                  - slope*sum(gates_real[(2*I-1),(2*J-1), n_g] * z_onoff_var[n_g,d] for n_g = 1:num_gates) - c == 0)
#                     end
#                     num_facets += 1

#                 elseif (length(vertices_coord) > 2) && (length(vertices_coord) <= max_ex_pt)
                
#                     vertices = zeros(length(vertices_coord),2)

#                     ctr = 1
#                     for l in vertices_coord
#                         vertices[ctr,:] = [l[1] l[2]]
#                         ctr += 1
#                     end
                    
#                     # vertices_convex_hull = Polyhedra.planar_hull(Polyhedra.vrep(vertices)).points.points
#                     vertices_convex_hull = vertices[QHull.chull(vertices).vertices, :]
#                     num_ex_pt = size(vertices_convex_hull)[1]

#                     for i=1:num_ex_pt
#                         v1 = vertices_convex_hull[i,:]

#                         if i == num_ex_pt
#                             v2 = vertices_convex_hull[1,:]
#                         else 
#                             v2 = vertices_convex_hull[i+1,:]
#                         end

#                         # Test-vertex for half-space directionality
#                         if i == (num_ex_pt - 1)
#                             v3 = vertices_convex_hull[1,:]
#                         elseif i == num_ex_pt
#                             v3 = vertices_convex_hull[2,:]
#                         else 
#                             v3 = vertices_convex_hull[i+2,:]
#                         end

#                         slope, intercept = QCO._get_constraint_slope_intercept(v1, v2)
                        
#                         # Facets of the hull
#                         if !isinf(slope)

#                             if v3[2] - slope*v3[1] - intercept <= -1E-6

#                                 if isapprox(abs(slope), 0, atol=1E-6)
                                    
#                                     JuMP.@constraint(qcm.model, [d=1:depth], sum(gates_real[(2*I-1),(2*J), n_g] * z_onoff_var[n_g,d] for n_g = 1:num_gates) - intercept <= 0)
#                                 else
                                    
#                                     JuMP.@constraint(qcm.model, [d=1:depth], sum(gates_real[(2*I-1),(2*J), n_g] * z_onoff_var[n_g,d] for n_g = 1:num_gates) 
#                                                                         - slope*(sum(gates_real[(2*I-1),(2*J-1), n_g] * z_onoff_var[n_g,d] for n_g = 1:num_gates)) - intercept <= 0)
#                                 end
#                                 num_facets += 1

#                             elseif v3[2] - slope*v3[1] - intercept >= 1E-6

#                                 if isapprox(abs(slope), 0, atol=1E-6)
                                    
#                                     JuMP.@constraint(qcm.model, [d=1:depth], sum(gates_real[(2*I-1),(2*J), n_g] * z_onoff_var[n_g,d] for n_g = 1:num_gates) - intercept >= 0)
#                                 else
                                    
#                                     JuMP.@constraint(qcm.model, [d=1:depth], sum(gates_real[(2*I-1),(2*J), n_g] * z_onoff_var[n_g,d] for n_g = 1:num_gates) 
#                                                                         - slope*(sum(gates_real[(2*I-1),(2*J-1), n_g] * z_onoff_var[n_g,d] for n_g = 1:num_gates)) - intercept >= 0)
#                                 end
#                                 num_facets += 1
                                
#                             else 
#                                 Memento.warn(_LOGGER, "Indeterminate direction for the planar-hull cut")
#                             end

#                         else isinf(slope)

#                             if v3[1] >= v1[1] + 1E-6
                                
#                                 JuMP.@constraint(qcm.model, [d=1:depth], sum(gates_real[(2*I-1),(2*J-1), n_g] * z_onoff_var[n_g,d] for n_g = 1:num_gates) >= v1[1])
#                             elseif v3[1] <= v1[1] - 1E-6
                                
#                                 JuMP.@constraint(qcm.model, [d=1:depth], sum(gates_real[(2*I-1),(2*J-1), n_g] * z_onoff_var[n_g,d] for n_g = 1:num_gates) <= v1[1])
#                             else
#                                 Memento.warn(_LOGGER, "Indeterminate direction for the convex-hull cut")
#                             end
#                             num_facets += 1
                                
#                         end
#                     end
#                 end
#             end
#         end
        
#         if num_facets > 0
#             Memento.info(_LOGGER, "Applied $num_facets planar-hull cuts per depth of the decomposition")
#         end
        
#     end
    
#     return
end