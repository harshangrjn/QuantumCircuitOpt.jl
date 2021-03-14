#------------------------------------------------------------#
# Build the objective function for QuantumCircuitModel here  #
#------------------------------------------------------------#

function objective_QCModel(qcm::QuantumCircuitModel)
    JuMP.@objective(qcm.model, Min, 1)
    return
end

