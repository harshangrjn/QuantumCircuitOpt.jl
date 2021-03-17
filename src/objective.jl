#------------------------------------------------------------#
# Build the objective function for QuantumCircuitModel here  #
#------------------------------------------------------------#

function objective_minimize_total_depth(qcm::QuantumCircuitModel)
    n_gates = size(qcm.data["M_real"])[3]
    depth   = qcm.data["depth"]
    
    if !isempty(findall(x -> startswith(x, "Identity"), qcm.data["elementary_gates"]))
        
        identity_id = Int64[]

        for i in keys(qcm.data["M_complex_dict"])
            if qcm.data["M_complex_dict"][i]["type"] == "Identity"
                push!(identity_id , parse(Int64, i))
            end
        end

        JuMP.@objective(qcm.model, Max, sum(qcm.variables[:z_onoff_var][n,d] for n in identity_id, d=1:depth))

    else

        Memento.warn(_LOGGER, "Switching to a feasibility problem since 
                               Identity gate is not part of the input elementary gates")
        JuMP.@objective(qcm.model, Min, 1)          
    end                              
    
    return
end

function objective_minimize_cnot_gates(qcm::QuantumCircuitModel)
    n_gates = size(qcm.data["M_real"])[3]
    depth   = qcm.data["depth"] 
    
    # Note: The below objective minimizes both cnot_12 and cnot_21
    if !isempty(findall(x -> startswith(x, "cnot"), qcm.data["elementary_gates"]))
        cnot_id = Int64[]

        for i in keys(qcm.data["M_complex_dict"])
            if qcm.data["M_complex_dict"][i]["type"] == "cnot_12"
                push!(cnot_id, parse(Int64, i))
            end
            if qcm.data["M_complex_dict"][i]["type"] == "cnot_21"
                push!(cnot_id, parse(Int64, i))
            end
        end

        JuMP.@objective(qcm.model, Min, sum(qcm.variables[:z_onoff_var][n,d] for n in cnot_id, d=1:depth))

    else

        Memento.warn(_LOGGER, "Switching the objective to minimize the total depth since 
                                        CNOT gate is not part of the input elementary gates")
        QCO.objective_minimize_total_depth(qcm)
    end

    return
end
 