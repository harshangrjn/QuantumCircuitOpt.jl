import QuantumCircuitOpt as QCOpt
using LinearAlgebra

using JuMP
using Gurobi
# using CPLEX
# using HiGHS

include("optimizers.jl")
[include("$(i)qubit_gates.jl") for i in 2:5]
include("parametrized_gates.jl")
include("decompose_all_gates.jl")

qcopt_optimizer = get_gurobi(solver_log = false)

result = Dict{String,Any}()

function make_gates(gates::Vector{String}, num_qubits::Int)
    for gate_name in gates
        @eval $(Symbol(gate_name)) = QCOpt.unitary($gate_name, $num_qubits)
    end
    return nothing
end

function target_gate(gates::Vector{String})
    gate_product = eval(Symbol(gates[1])) # start with first gate
    for gate_name in gates[2:end]
        gate_product *= eval(Symbol(gate_name))
    end
    return gate_product
end

function gate_order(gates::Vector{String}, z::Array{Float64,2})
    N, D = size(z)
    result = Vector{String}(undef, D)
    for d in 1:D
        idx = findfirst(n -> z[n, d] == 1.0, 1:N)
        result[d] = gates[idx]
    end
    return result
end

function toffoli_left()

    return Dict{String, Any}(
    "num_qubits" => 3,
    "maximum_depth" => 7,
    "elementary_gates" => ["H_3", "T_1", "T_2", "T_3", "Tdagger_1", "Tdagger_2", "Tdagger_3", "CNot_1_2", "CNot_2_3", "CNot_1_3", "Identity"],
    #"target_gate" => target_gate("H_3" , "CNot_2_3" , "Tdagger_3", "CNot_1_3" , "T_3" , "CNot_2_3" , "Tdagger_3"),
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal"
    )
end


params = toffoli_left()

make_gates(params["elementary_gates"],params["num_qubits"])
params["target_gate"] = target_gate(["H_3" , "CNot_2_3" , "Tdagger_3", "CNot_1_3" , "T_3" , "CNot_2_3" , "Tdagger_3"])
println(params)


model_options = Dict{Symbol, Any}(
        :model_type => "compact_formulation",
        :convex_hull_gate_constraints => false,
        :idempotent_gate_constraints  => true,
        :fix_unitary_variables        => false,
        :unitary_complex_conjugate    => 2,
)

global result = QCOpt.run_QCModel(params, qcopt_optimizer; options = model_options)
z_bin_vars = result["solution"]["z_bin_var"]
elementary_gates = params["elementary_gates"]
println(params["elementary_gates"])
println(z_bin_vars)
println(gate_order(elementary_gates,z_bin_vars))


#=
function parameters()
    return Dict{String, Any}(
        "num_qubits" => 3,
        "maximum_depth" => 7,
        "elementary_gates" => ["H_3", "T_1", "T_2", "T_3", 
                               "Tdagger_1", "Tdagger_2", "Tdagger_3", 
                               "CNot_1_2", "CNot_2_3", "CNot_1_3", "Identity"],
        "objective" => "minimize_depth",
        "decomposition_type" => "exact_optimal"
    )
end


params = parameters()

# Make the elementary gates
gate_dict = Dict{String, Any}()
make_gates(params["elementary_gates"], params["num_qubits"], gate_dict)
full_target = ["H_3", "CNot_2_3", "Tdagger_3", "CNot_1_3", "T_3", "CNot_2_3", "Tdagger_3"]

# Add the target gate using the gate dictionary
params["target_gate"] = target_gate(
    #["H_3", "CNot_2_3", "Tdagger_3", "CNot_1_3", "T_3", "CNot_2_3", "Tdagger_3"], 
    ["CNot_1_2", "Identity", "T_2"],
    gate_dict
)

println(params)
=#