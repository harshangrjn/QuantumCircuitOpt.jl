import QuantumCircuitOpt as QCOpt
using LinearAlgebra
using Printf

using JuMP
using Gurobi

include("optimizers.jl")

function print_matrix_pretty(mat::Matrix{ComplexF64}; digits=3)
    for row in eachrow(mat)
        for val in row
            real_part = round(real(val), digits=digits)
            imag_part = round(imag(val), digits=digits)
            sign = imag_part < 0 ? "-" : "+"
            @printf("%7.3f %s %6.3fi  ", real_part, sign, abs(imag_part))
        end
        println()
    end
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

function push_identities(gates::Vector{String})
    others     = filter(x -> x != "Identity", gates)
    identities = filter(x -> x == "Identity", gates)
    return vcat(others, identities)
end

function delete_identities(gates::Vector{String})
    return filter(x -> x != "Identity", gates)
end

function make_gates(gates::Vector{String}, num_qubits::Int)
    gate_dict = Dict{String, Any}()
    for gate_name in gates
        gate_dict[gate_name] = QCOpt.unitary(gate_name, num_qubits)
    end
    return gate_dict
end

function make_target_gate(gates::Vector{String}, gate_dict::Dict{String,Any})
    gate_product = gate_dict[gates[1]] # start with first gate
    for gate_name in gates[2:end]
        gate_product *= gate_dict[gate_name]
    end
    #return gate_product
    return round.(gate_product; digits=8)
end

function extract_qubits(gate::String)
    # This assumes gates are named like "CNot_2_3" or "Tdagger_3"
    # Splits at underscores, parses numbers
    parts = split(gate, "_")
    return [parse(Int, s) for s in parts[2:end] if occursin(r"^\d+$", s)]
end

function extract_qubits_from_all(gates::Vector{String})
    qubits = Set{Int}()
    for gate in gates
        union!(qubits, extract_qubits(gate))
    end
    return sort(collect(qubits))
end

function gates_on_only_these_qubits(elementary_gates::Vector{String}, qubits::Vector{Int})
    selected = String[]
    qubit_set = Set(qubits)
    for gate in elementary_gates
        # Always add Identity at the end, not in the loop
        if gate == "Identity"
            continue
        end
        gate_qubits = extract_qubits(gate)
        # Include the gate if all its qubits are in the allowed set
        if all(q -> q in qubit_set, gate_qubits)
            push!(selected, gate)
        end
    end
    push!(selected, "Identity")
    return selected
end

function gates_on_qubits_up_to(t::Vector{String}, q::Vector{Int}, i::Int)
    # Returns all gates from t up to t[i] that act on qubits from q
    selected_gates = String[]
    for j in 1:i
        gate = t[j]
        gate_qubits = extract_qubits(gate)
        if !isempty(intersect(gate_qubits, q))
            push!(selected_gates, gate)
        end
    end
    return selected_gates
end

function recursive_gates_on_qubits_up_to(t::Vector{String}, q::Vector{Int}, i::Int)
    selected_gates = gates_on_qubits_up_to(t, q, i)
    new_q = extract_qubits_from_all(selected_gates)
    if length(new_q) > length(q)
        # If we found new qubits, recurse with expanded q
        return recursive_gates_on_qubits_up_to(t, new_q, i)
    else
        return selected_gates
    end
end

function find_first_block(t::Vector{String}, window_length::Int, max_qubits::Int)
    if isempty(t)
        error("Input sequence t is empty.")
    end
    # Step 1: get qubits from the first gate
    q = extract_qubits(t[1])

    i = 1
    best_sequence = []
    reached_end = false

    while true
        if i > length(t)
            i = length(t)
            reached_end = true
        end

        seq = recursive_gates_on_qubits_up_to(t, q, i)

        if length(seq) > window_length
            # Too big, return the previous sequence
            return best_sequence
        elseif length(extract_qubits_from_all(seq)) > max_qubits
            # Too many qubits in the new sequance
            return best_sequence
        elseif length(seq) == window_length
            return seq
        else
            best_sequence = seq
            if reached_end
                println("Reached the end of t before filling the window (length=$(length(best_sequence)) < window_length=$window_length).")
                return best_sequence
            end
            i += 1
        end
    end
end

function remove_block(full_target::Vector{String}, current_target::Vector{String})
    remaining = copy(full_target)
    for gate in current_target
        idx = findfirst(==(gate), remaining)
        if idx !== nothing
            deleteat!(remaining, idx)
        else 
            println("You are trying to delete something that is not there")
        end
    end
    return remaining
end

function translate_gate_list(t::Vector{String}, mapping::Dict{Int,Int})
    translated = String[]
    for gate in t
        parts = split(gate, "_")
        if length(parts) == 1
            push!(translated, gate)
        else
            new_parts = [parts[1]]
            for s in parts[2:end]
                if occursin(r"^\d+$", s)
                    orig_qubit = parse(Int, s)
                    push!(new_parts, string(mapping[orig_qubit]))
                else
                    push!(new_parts, s)
                end
            end
            push!(translated, join(new_parts, "_"))
        end
    end
    return translated
end

function back_translate_gate_list(t::Vector{String}, mapping::Dict{Int,Int})
    # Build inverse mapping: new_index => original_qubit
    inv_mapping = Dict(v => k for (k, v) in mapping)
    translated = String[]
    for gate in t
        parts = split(gate, "_")
        if length(parts) == 1
            push!(translated, gate)
        else
            new_parts = [parts[1]]
            for s in parts[2:end]
                if occursin(r"^\d+$", s)
                    new_index = parse(Int, s)
                    push!(new_parts, string(inv_mapping[new_index]))
                else
                    push!(new_parts, s)
                end
            end
            push!(translated, join(new_parts, "_"))
        end
    end
    return translated
end

function one_optimization(target_sequence::Vector{String}, # target for this optimization
                          elementary_gates::Vector{String}, # array of elementary gates
                          optimizer;
                          options=Dict{Symbol,Any}())

    num_qubits = length(extract_qubits_from_all(target_sequence))
    gate_dict = make_gates(elementary_gates, num_qubits)   

    # Build target unitary
    targetU = make_target_gate(target_sequence, gate_dict)

    # If target is identity on num_qubits, skip optimization
    d = size(targetU, 1)
    if d == (1 << num_qubits)
        idmat = Matrix{eltype(targetU)}(LinearAlgebra.I, d, d)
        if isapprox(targetU, idmat; atol=1e-12, rtol=1e-12)
            return ["Identity"]
        end
    end

    params=Dict{String, Any}(
        "num_qubits" => num_qubits,
        "maximum_depth" => length(target_sequence),
        "elementary_gates" => elementary_gates,
        "target_gate" => targetU,
        "objective" => "minimize_depth",
        "decomposition_type" => "exact_optimal",

    )

    # Run optimization
    result = Dict{String,Any}()
    result = QCOpt.run_QCModel(params, optimizer; options = options)

    # Process output 
    z_bin_vars = round.(result["solution"]["z_bin_var"]; digits=4)
    println(z_bin_vars)
    ordered_gates = gate_order(elementary_gates, z_bin_vars)

    return ordered_gates
    
end

function rolling_horizon(target_sequence::Vector{String}, # full target
                         window_length::Int,              # window length
                         max_qubits::Int,                 # max number of qubits in a window
                         elementary_gates::Vector{String}, # array of elementary gates
                         optimizer;
                         accept_window::Int = 1,
                         options=Dict{Symbol,Any}())

    
    num_qubits = length(extract_qubits_from_all(target_sequence))
    gate_dict = make_gates(elementary_gates, num_qubits)   

    output_gates = String[]
    full_target = push_identities(target_sequence)

    while true 
        println("current output: ", output_gates)
        println("full_target remaining: ", full_target)

        if length(full_target) < 1 || all(x -> x == "Identity", full_target)
            # done optimization
            println("Done rolling!")
            break
        end

        current_target = find_first_block(push_identities(full_target), window_length, max_qubits)
        println("window gates: ", current_target)

        current_qubits = extract_qubits_from_all(current_target)
        if length(current_qubits) == 1
            full_target = remove_block(full_target, current_target)
            append!(output_gates, current_target)
        elseif length(current_target) == 1
            full_target = remove_block(full_target, current_target)
            append!(output_gates, current_target)
        else
            mapping = Dict(q => i for (i, q) in enumerate(current_qubits))
            current_elementary_gates = gates_on_only_these_qubits(elementary_gates, current_qubits)
            println("current_elementary_gates:", translate_gate_list(current_elementary_gates, mapping))

            println("translated window gates: ",translate_gate_list(current_target, mapping))

            current_optimized = one_optimization(translate_gate_list(current_target, mapping), translate_gate_list(current_elementary_gates, mapping), qcopt_optimizer; options = model_options)
            current_optimized = back_translate_gate_list(current_optimized, mapping)
            println("ordered_result: ", current_optimized)

            if length(current_target) == length(full_target)
                # done optimization
                println("Done rolling!")
                append!(output_gates, current_optimized[1:end])
                break
            else
                full_target = remove_block(full_target, current_target)
                if length(current_optimized) > accept_window
                    append!(output_gates, current_optimized[1:accept_window])
                    full_target = delete_identities(vcat(current_optimized[accept_window+1:end], full_target))
                else
                    append!(output_gates, current_optimized)
                    # full_target stays the same
                end
            end
        end
    end

    output_gates = delete_identities(output_gates)

    println("compilation successfull: ", (isapprox(make_target_gate(target_sequence, gate_dict), make_target_gate(output_gates, gate_dict); atol=1e-7)))

    return output_gates

end

###########################  some TESTS   ###################################

function test_extract_qubits_from_all()
    gates = ["H_3", "H_1", "CNot_2_3", "T_6", "CNot_1_2"]
    result = extract_qubits_from_all(gates)
    println(result)  
end

function test_recursive_gates_on_qubits_up_to()
    t = ["H_3", "H_1", "H_2", "H_5", "CNot_2_3", "Tdagger_3", "CNot_1_3", "T_3", "CNot_2_3", "Tdagger_3", "CNot_1_3", "T_3", "T_2", "H_3", "CNot_1_2", "T_1", "Tdagger_2", "CNot_1_2"]
    q = [3]
    i = 7
    selected = recursive_gates_on_qubits_up_to(t, q, i)
    println(selected)
end

function test_find_first_block()
    # Example usage
    t             = [  "CNot_5_4", "H_3", "X_4", "H_1", "Identity", "H_2", "H_5", "CNot_2_3", "Tdagger_3", "CNot_1_3", "T_3", "CNot_2_3", "Tdagger_3", 
                        "CNot_1_3", "T_3", "T_2", "H_3", "CNot_1_2", "T_1", "Tdagger_2", "CNot_1_2"]
    window_length = 10
    max_qubits    = 4

    block = find_first_block(t, window_length, max_qubits)
    println("Block: ", block)
end

function test_gates_on_only_these_qubits()
    elementary_gates = ["H_1", "H_2", "H_3", 
                    "T_1", "T_2", "T_3", 
                    "Tdagger_1", "Tdagger_2", "Tdagger_3", 
                    "CNot_1_2", "CNot_2_3", "CNot_1_3", "CNot_2_1",
                    "CNot_3_2", "CNot_3_1", "Identity"]

    qubits = [2, 1]
    filtered = gates_on_only_these_qubits(elementary_gates, qubits)
    println(filtered)
    # Output: ["H_1", "H_2", "T_1", "T_2", "Tdagger_1", "Tdagger_2", "CNot_1_2", "CNot_2_1", "Identity"]
end

########################################     RUNS     ######################################

qcopt_optimizer = get_gurobi(solver_log = false)

model_options = Dict{Symbol, Any}(
        :model_type => "compact_formulation",
        :convex_hull_gate_constraints => false,
        :idempotent_gate_constraints  => true,
        :fix_unitary_variables        => false,
        :unitary_complex_conjugate    => 2,
)

elementary_gates = ["H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7",
                    "X_1", "X_2", "X_3", "X_4", "X_5", "X_6", "X_7",
                    "T_1", "T_2", "T_3", "T_4", "T_5", "T_6", "T_7",
                    "Tdagger_1", "Tdagger_2", "Tdagger_3", "Tdagger_4", "Tdagger_5", "Tdagger_6", "Tdagger_7", 
                    "CNot_1_2", "CNot_1_3", "CNot_1_4",  "CNot_1_5", "CNot_1_6", "CNot_1_7",
                    "CNot_2_1", "CNot_2_3", "CNot_2_4",  "CNot_2_5", "CNot_2_6", "CNot_2_7",
                    "CNot_3_2", "CNot_3_1", "CNot_3_4",  "CNot_3_5", "CNot_3_6", "CNot_3_7",
                    "CNot_4_2", "CNot_4_3", "CNot_4_1",  "CNot_4_5", "CNot_4_6", "CNot_4_7",
                    "CNot_5_2", "CNot_5_3", "CNot_5_4",  "CNot_5_1", "CNot_5_6", "CNot_5_7",
                    "CNot_6_2", "CNot_6_3", "CNot_6_4",  "CNot_6_5", "CNot_6_1", "CNot_6_7",
                    "CNot_7_2", "CNot_7_3", "CNot_7_4",  "CNot_7_5", "CNot_7_1", "CNot_7_6",
                    "Identity"]

elementary_gates = [
                    "S_1", "S_2", "S_3", "S_4",
                    "Sdagger_1", "Sdagger_2", "Sdagger_3", "Sdagger_4",  
                    "CNot_1_2", "CNot_1_3", "CNot_1_4", # "CNot_1_5", "CNot_1_6", "CNot_1_7",
                    "CNot_2_1", "CNot_2_3", "CNot_2_4", # "CNot_2_5", "CNot_2_6", "CNot_2_7",
                    "CNot_3_2", "CNot_3_1", "CNot_3_4", # "CNot_3_5", "CNot_3_6", "CNot_3_7",
                    "CNot_4_2", "CNot_4_3", "CNot_4_1", # "CNot_4_5", "CNot_4_6", "CNot_4_7",
                  #  "CNot_5_2", "CNot_5_3", "CNot_5_4",  "CNot_5_1", "CNot_5_6", "CNot_5_7",
                  #  "CNot_6_2", "CNot_6_3", "CNot_6_4",  "CNot_6_5", "CNot_6_1", "CNot_6_7",
                  #  "CNot_7_2", "CNot_7_3", "CNot_7_4",  "CNot_7_5", "CNot_7_1", "CNot_7_6",
                    "Identity"]

elementary_gates = [
                    "S_1", "S_2", "S_3", "S_4",
                    "Sdagger_1", "Sdagger_2", "Sdagger_3", "Sdagger_4",  
                    "CNot_1_2", "CNot_1_3", "CNot_1_4", # "CNot_1_5", "CNot_1_6", "CNot_1_7",
                    "CNot_2_1", "CNot_2_3", "CNot_2_4", # "CNot_2_5", "CNot_2_6", "CNot_2_7",
                    "CNot_3_2", "CNot_3_1", "CNot_3_4", # "CNot_3_5", "CNot_3_6", "CNot_3_7",
                    "CNot_4_2", "CNot_4_3", "CNot_4_1", # "CNot_4_5", "CNot_4_6", "CNot_4_7",
                  #  "CNot_5_2", "CNot_5_3", "CNot_5_4",  "CNot_5_1", "CNot_5_6", "CNot_5_7",
                  #  "CNot_6_2", "CNot_6_3", "CNot_6_4",  "CNot_6_5", "CNot_6_1", "CNot_6_7",
                  #  "CNot_7_2", "CNot_7_3", "CNot_7_4",  "CNot_7_5", "CNot_7_1", "CNot_7_6",
                    "Identity"]

elementary_gates = [
    # Single-qubit
    "S_1","S_2","S_3","S_4","S_5","S_6","S_7","S_8",
    "Sdagger_1","Sdagger_2","Sdagger_3","Sdagger_4","Sdagger_5","Sdagger_6","Sdagger_7","Sdagger_8",
    "H_1","H_2","H_3","H_4","H_5","H_6","H_7","H_8",

    # Two-qubit CNOTs on the 8-qubit ring (both directions)
    "CNot_1_2","CNot_2_1",
    "CNot_2_3","CNot_3_2",
    "CNot_3_4","CNot_4_3",
    "CNot_4_5","CNot_5_4",
    "CNot_5_6","CNot_6_5",
    "CNot_6_7","CNot_7_6",
    "CNot_7_8","CNot_8_7",
    "CNot_8_1","CNot_1_8",
]

elementary_gates = ["H_1", "H_2", "H_3", "H_4", #"H_5", "H_6", "H_7",
                    #"X_1", "X_2", "X_3", "X_4", "X_5", "X_6", "X_7",
                    "T_1", "T_2", "T_3", "T_4", #"T_5", "T_6", "T_7",
                    "Tdagger_1", "Tdagger_2", "Tdagger_3", "Tdagger_4", #"Tdagger_5", "Tdagger_6", "Tdagger_7", 
                    "CNot_1_2", "CNot_1_3", "CNot_1_4", # "CNot_1_5", "CNot_1_6", "CNot_1_7",
                    "CNot_2_1", "CNot_2_3", "CNot_2_4", # "CNot_2_5", "CNot_2_6", "CNot_2_7",
                    "CNot_3_2", "CNot_3_1", "CNot_3_4", # "CNot_3_5", "CNot_3_6", "CNot_3_7",
                    "CNot_4_2", "CNot_4_3", "CNot_4_1", # "CNot_4_5", "CNot_4_6", "CNot_4_7",
                  #  "CNot_5_2", "CNot_5_3", "CNot_5_4",  "CNot_5_1", "CNot_5_6", "CNot_5_7",
                  #  "CNot_6_2", "CNot_6_3", "CNot_6_4",  "CNot_6_5", "CNot_6_1", "CNot_6_7",
                  #  "CNot_7_2", "CNot_7_3", "CNot_7_4",  "CNot_7_5", "CNot_7_1", "CNot_7_6",
                    "Identity"]

#full_target = ["H_3", "H_4", "CNot_2_3", "Identity", "Tdagger_3", "H_3", "H_1", "CNot_3_1", "H_3", "H_1", "T_3", "CNot_2_3", "Tdagger_3", "CNot_1_3", "T_3", "T_2", "H_3", "CNot_1_2", "T_1", "Tdagger_2", "CNot_1_2"]
#full_target = ["T_2", "H_3", "H_1", "CNot_2_3", "H_3", "H_2"]
#full_target = ["H_3", "CNot_2_3", "Tdagger_3", "CNot_1_3", "T_3", "CNot_2_3", "Tdagger_3", "CNot_1_3", "T_3", "T_2" ]



full_target = [
    "H_5",
    "H_6",
    "CNot_2_6",
    "Tdagger_6",
    "CNot_1_6",
    "T_6",
    "CNot_2_6",
    "T_2",
    "Tdagger_6",
    "CNot_1_6",
    "CNot_1_2",
    "T_1",
    "Tdagger_2",
    "CNot_1_2",
    "T_6",
    "H_6",
    "H_7",
    "CNot_4_7",
    "Tdagger_7",
    "CNot_3_7",
    "T_7",
    "CNot_4_7",
    "T_4",
    "Tdagger_7",
    "CNot_3_7",
    "CNot_3_4",
    "T_3",
    "Tdagger_4",
    "CNot_3_4",
    "T_7",
    "H_7",
    "CNot_7_5",
    "Tdagger_5",
    "CNot_6_5",
    "T_5",
    "CNot_7_5",
    "Tdagger_5",
    "CNot_6_5",
    "T_5",
    "H_5",
    "T_7",
    "CNot_6_7",
    "T_6",
    "Tdagger_7",
    "CNot_6_7",
    "H_6",
    "CNot_2_6",
    "Tdagger_6",
    "CNot_1_6",
    "T_6",
    "CNot_2_6",
    "T_2",
    "Tdagger_6",
    "CNot_1_6",
    "CNot_1_2",
    "T_1",
    "Tdagger_2",
    "CNot_1_2",
    "T_6",
    "H_6",
    "H_7",
    "CNot_4_7",
    "Tdagger_7",
    "CNot_3_7",
    "T_7",
    "CNot_4_7",
    "T_4",
    "Tdagger_7",
    "CNot_3_7",
    "CNot_3_4",
    "T_3",
    "Tdagger_4",
    "CNot_3_4",
    "T_7",
    "H_7",
]

full_target = [ 
    "CNot_1_3",
    "CNot_2_3",
    "S_3",
    "CNot_2_3",
    "CNot_1_3",

    "CNot_1_4",
    "CNot_2_4",
    "S_4",
    "CNot_2_4",
    "CNot_1_4",

    "CNot_1_4",
    "CNot_3_4",
    "S_4",
    "CNot_3_4",
    "CNot_1_4",

    "CNot_2_4",
    "CNot_3_4",
    "S_4",
    "CNot_3_4",
    "CNot_2_4"
]

full_target = [
  # --- Layer 1: ZZ (color A)
  "CNot_1_2","S_2","CNot_1_2",
  "CNot_3_4","S_4","CNot_3_4",
  "CNot_5_6","S_6","CNot_5_6",
  "CNot_7_8","S_8","CNot_7_8",
  # --- Layer 1: ZZ (color B)
  "CNot_2_3","S_3","CNot_2_3",
  "CNot_4_5","S_5","CNot_4_5",
  "CNot_6_7","S_7","CNot_6_7",
  "CNot_8_1","S_1","CNot_8_1",
  # --- Layer 1: Mixer (V = H S H on each qubit)
  "H_1","S_1","H_1","H_2","S_2","H_2","H_3","S_3","H_3","H_4","S_4","H_4",
  "H_5","S_5","H_5","H_6","S_6","H_6","H_7","S_7","H_7","H_8","S_8","H_8",

  # --- Layer 2: ZZ (color A)
  "CNot_1_2","S_2","CNot_1_2",
  "CNot_3_4","S_4","CNot_3_4",
  "CNot_5_6","S_6","CNot_5_6",
  "CNot_7_8","S_8","CNot_7_8",
  # --- Layer 2: ZZ (color B)
  "CNot_2_3","S_3","CNot_2_3",
  "CNot_4_5","S_5","CNot_4_5",
  "CNot_6_7","S_7","CNot_6_7",
  "CNot_8_1","S_1","CNot_8_1",
  # --- Layer 2: Mixer
  "H_1","S_1","H_1","H_2","S_2","H_2","H_3","S_3","H_3","H_4","S_4","H_4",
  "H_5","S_5","H_5","H_6","S_6","H_6","H_7","S_7","H_7","H_8","S_8","H_8",

  # --- Layer 3: ZZ (color A)
  "CNot_1_2","S_2","CNot_1_2",
  "CNot_3_4","S_4","CNot_3_4",
  "CNot_5_6","S_6","CNot_5_6",
  "CNot_7_8","S_8","CNot_7_8",
  # --- Layer 3: ZZ (color B)
  "CNot_2_3","S_3","CNot_2_3",
  "CNot_4_5","S_5","CNot_4_5",
  "CNot_6_7","S_7","CNot_6_7",
  "CNot_8_1","S_1","CNot_8_1",
  # --- Layer 3: Mixer
  "H_1","S_1","H_1","H_2","S_2","H_2","H_3","S_3","H_3","H_4","S_4","H_4",
  "H_5","S_5","H_5","H_6","S_6","H_6","H_7","S_7","H_7","H_8","S_8","H_8",
]

full_target = ["CNot_1_3", "CNot_3_4", "CNot_4_1", "H_3", "T_2", "T_4", "T_3", "CNot_4_3", "Tdagger_3", "CNot_2_3", "T_3", "CNot_4_3", "Tdagger_3", "CNot_2_3", "H_3", "CNot_2_4", "Tdagger_2", "CNot_2_4", "CNot_1_2", "H_2", "T_3", "T_4", "T_2", "CNot_4_2", "Tdagger_2", "CNot_3_2", "T_2", "CNot_4_2", "Tdagger_2", "CNot_3_2", "H_2", "CNot_3_4", "Tdagger_3", "CNot_3_4", "H_4", "CNot_1_4", "Tdagger_4", "CNot_1_4", "T_1", "T_4", "CNot_1_2", "Tdagger_4", "Tdagger_2", "CNot_2_4", "T_4", "CNot_2_4", "CNot_1_2", "CNot_2_4", "Tdagger_4", "CNot_2_4", "T_2", "T_4", "CNot_2_3", "Tdagger_4", "Tdagger_3", "CNot_3_4", "T_4", "CNot_3_4", "CNot_1_3", "CNot_3_4", "Tdagger_4", "CNot_3_4", "T_3", "T_4", "CNot_2_3", "Tdagger_4", "Tdagger_3", "CNot_3_4", "T_4", "CNot_3_4", "CNot_1_3", "CNot_3_4", "Tdagger_4", "CNot_3_4", "T_3", "T_4", "H_4", "CNot_3_1", "H_2", "T_2", "T_2", "T_2", "T_2", "H_2", "H_3", "T_3", "T_3", "T_3", "T_3", "H_3", "CNot_1_4", "H_3", "T_2", "T_4", "T_3", "CNot_4_3", "Tdagger_3", "CNot_2_3", "T_3", "CNot_4_3", "Tdagger_3", "CNot_2_3", "H_3", "CNot_2_4", "Tdagger_2", "CNot_2_4", "H_1", "T_2", "T_3", "T_1", "CNot_3_1", "Tdagger_1", "CNot_2_1", "T_1", "CNot_3_1", "Tdagger_1", "CNot_2_1", "H_1", "CNot_2_3", "Tdagger_2", "CNot_2_3", "H_2", "T_1", "T_3", "T_2", "CNot_3_2", "Tdagger_2", "CNot_1_2", "T_2", "CNot_3_2", "Tdagger_2", "CNot_1_2", "H_2", "CNot_1_3", "Tdagger_1", "CNot_1_3", "H_3", "T_3", "T_3", "T_3", "T_3", "H_3"]

window_length   = 10  # how large is a block that I optimize
accept_window   = 7  # how many optimized gates we accept from a block
max_qubits      = 3  # max qubits in one window

optimized_target = rolling_horizon(full_target, window_length, max_qubits, elementary_gates, qcopt_optimizer;  accept_window = accept_window, options = model_options )