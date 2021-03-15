#------------------------------------------------------------#
# Build the objective function for QuantumCircuitModel here  #
#------------------------------------------------------------#

function objective_minimize_total_depth(qcm::QuantumCircuitModel)
    n_gates = size(qcm.data["M_real"])[3]
    depth   = qcm.data["depth"]
    identity_ids = findall(x -> startswith(x, "Identity"), qcm.data["elementary_gates"])
    
    if !isempty(identity_ids)
        JuMP.@objective(qcm.model, Max, sum(qcm.variables[:z_onoff_var][n,d] for n in identity_ids, d=1:depth))
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
    cnot_ids = findall(x -> startswith(x, "cnot"), qcm.data["elementary_gates"])    
    
    if !isempty(cnot_ids)
        JuMP.@objective(qcm.model, Min, sum(qcm.variables[:z_onoff_var][n,d] for n in cnot_ids, d=1:depth))
    else 
        Memento.warn(_LOGGER, "Switching the objective to minimize the total depth since 
                                        CNOT gate is not part of the input elementary gates")
        QCO.objective_minimize_total_depth(qcm)
    end
    return
end
 