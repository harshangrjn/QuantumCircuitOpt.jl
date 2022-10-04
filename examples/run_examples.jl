import QuantumCircuitOpt as QCOpt
using JuMP
using Gurobi
# using HiGHS

include("optimizers.jl")
include("2qubit_gates.jl")
include("3qubit_gates.jl")
include("4qubit_gates.jl")
include("5qubit_gates.jl")
include("decompose_all_gates.jl")

<<<<<<< HEAD
decompose_gates = ["decompose_hadamard",
                   "decompose_controlled_Z",
                   "decompose_controlled_V",
                   "decompose_controlled_H",
                   "decompose_controlled_H_with_R",
                   "decompose_magic", 
                   "decompose_magic_using_SHCnot",
                   "decompose_S", 
                   "decompose_revcnot", 
                   "decompose_revcnot_with_U", 
                   "decompose_swap",
                   "decompose_W",
                   "decompose_W_using_HCnot",
                   "decompose_GroverDiffusion_using_Clifford",
                   "decompose_GroverDiffusion_using_U3",
                   "decompose_iSwap",
                   "decompose_qft2_using_HT",
                   "decompose_RX_on_q3"]

decompose_gates = ["decompose_controlled_V"]
=======
# decompose_gates = ["iSwap"]
>>>>>>> 381b589ae3f07d25e43e78a63dfa40d287114925

#----------------------------------------------#
#      Quantum Circuit Optimization model      #
#----------------------------------------------#
qcopt_optimizer = get_gurobi()

result = Dict{String,Any}()
times = zeros(length(decompose_gates), 1)

for gates = 1:length(decompose_gates)
    params = getfield(Main, Symbol(decompose_gates[gates]))()

    model_options = Dict{Symbol, Any}(
        :convex_hull_gate_constraints => false,
        :idempotent_gate_constraints  => false,
        :unitary_constraints          => false,
        :fix_unitary_variables        => true,
    )

    global result = QCOpt.run_QCModel(params, qcopt_optimizer; options = model_options)
    times[gates] = result["solve_time"]
end
