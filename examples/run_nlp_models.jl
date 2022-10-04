import QuantumCircuitOpt as QCOpt
using JuMP
using Gurobi
using Ipopt

include("optimizers.jl")
include("2qubit_gates.jl")
include("3qubit_gates.jl")
include("4qubit_gates.jl")
include("5qubit_gates.jl")
include("decompose_all_gates.jl")

decompose_gates = ["revcnot"]

#----------------------------------------------#
#      Quantum Circuit Optimization model      #
#----------------------------------------------#
# qcopt_optimizer = get_gurobi_nonconvex()
# qcopt_optimizer = get_gurobi()
qcopt_optimizer = get_ipopt()

result = Dict{String,Any}()
num_runs = 5
result_tab = []
for k = 1:num_runs
    println("Run number = ", k)
    for gates in decompose_gates 
        
        params = getfield(Main, Symbol(gates))()

        model_options = Dict{Symbol, Any}(:convex_hull_gate_constraints => false,
                                          :unitary_constraints          => false,
                                          :fix_unitary_variables        => false,
                                          :num_dummy_vars               => 5,
                                          :model_type                   => "nlp_relaxation_1")
        
        global result = QCOpt.run_QCModel(params, qcopt_optimizer; options = model_options)
    end
    if !(result["termination_status"] in [MOI.LOCALLY_INFEASIBLE, 
                                          MOI.UNKNOWN_RESULT_STATUS, 
                                          MOI.ITERATION_LIMIT, 
                                          MOI.INVALID_MODEL])        
        push!(result_tab, [result["primal_status"], result["objective"], result["solve_time"]])
        @show result["termination_status"]
    end
end
