"""
    visualize_solution(results::Dict{String, Any}, data::Dict{String, Any}; gate_sequence = false)

Given dictionaries of results and data, and assuming that the optimization model had a feasible solution, 
this function aids in visualizing the optimal circuit decomposition.
"""
function visualize_solution(results::Dict{String, Any}, data::Dict{String, Any}; gate_sequence = false)

    global _header_color = :cyan 
    global _main_color   = :White

    if !(results["primal_status"] in [MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT]) || 
        (results["termination_status"] == MOI.INFEASIBLE)
         msg = results["termination_status"] == MOI.TIME_LIMIT ? 
               "Optimizer hits time limit with an infeasible primal status. Gate decomposition may be inaccurate" : 
               "Infeasible primal status. Gate decomposition may be inaccurate"
         Memento.warn(_LOGGER, msg)
         return
     end
    
    id_sequence = QCO._gate_id_sequence(results["solution"]["z_bin_var"], data["maximum_depth"])
    
    cirq_precompress, cirq_postcompress = QCO.build_circuit_layers(id_sequence, data["gates_dict"])

    cirq_precompress_depth = length(cirq_precompress)
    cirq_postcompress_depth = length(cirq_postcompress)
    
    if data["decomposition_type"] in ["exact_optimal", "exact_feasible", "optimal_global_phase"]
        # Validate pre-compressed circuit
        QCO.validate_circuit(data, cirq_precompress, global_phase = false, circuit_type = "precompress")
        # Validate post-compressed circuit only if compression occurred
        cirq_postcompress_depth < cirq_precompress_depth && 
            QCO.validate_circuit(data, cirq_postcompress, global_phase = false, circuit_type = "postcompress")
    end
    
    if !isempty(cirq_postcompress)

        printstyled("\n","=============================================================================","\n"; color = _main_color)
        printstyled("QuantumCircuitOpt version: ", Pkg.TOML.parse(read(string(pkgdir(QCO), "/Project.toml"), String))["version"], "\n"; color = _header_color, bold = true)

        printstyled("\n","Quantum Circuit Data"; color = _header_color, bold = true)
        
        printstyled("\n","  ","Number of qubits      : ", data["num_qubits"], "\n"; color = _main_color)

        printstyled("  ","Maximum circuit depth : ", data["maximum_depth"],"\n"; color = _main_color)

        printstyled("  ","Decomposition type    : ", data["decomposition_type"],"\n"; color = _main_color)

        printstyled("  ","MIP optimizer         : ", results["optimizer"],"\n"; color = _main_color)
        
        printstyled("  ","Number of input gates : ",size(data["gates_real"])[3]," (presolved)","\n"; color = _main_color)
        
        QCO.print_elementary_gates(data["elementary_gates"]; color = _main_color)

        if "discretization" in keys(data)
            for i in keys(data["discretization"])
                printstyled("    ","$(replace(i, "_discretization" => "")) : ", ceil.(rad2deg.(data["discretization"][i]), digits = 1),"\n"; color = _main_color)
            end
        end

        printstyled("\n","Optimal Circuit","\n"; color = _header_color, bold = true)

        QCO.print_circuit(cirq_postcompress_depth < cirq_precompress_depth ? cirq_postcompress : cirq_precompress, data["num_qubits"])
        printstyled("\n")
        
        if data["decomposition_type"] == "approximate"
            printstyled("  ","||Decomposition error||₂: ", round(LA.norm(results["solution"]["slack_var"]), digits = 10),"\n"; color = _main_color)
        end

        if data["objective"] == "minimize_depth"

            if length(data["identity_idx"]) >= 1 && (data["decomposition_type"] !== "exact_feasible") && !(results["termination_status"] == MOI.TIME_LIMIT)
                printstyled("  ","Optimal number of gates : ", cirq_precompress_depth,"\n"; color = _main_color)
                printstyled("  ","Optimal circuit depth   : ", cirq_postcompress_depth,"\n"; color = _main_color)
            else 
                printstyled("  ","Circuit depth: ", cirq_postcompress_depth,"\n"; color = _main_color)
            end

        elseif data["objective"] in ["minimize_cnot", "minimize_T"]
            gate_type = data["objective"] == "minimize_cnot" ? "CNOT" : "T"
            idx_key = data["objective"] == "minimize_cnot" ? "cnot_idx" : "T_idx"

            if !isempty(data[idx_key])
                if data["decomposition_type"] in ["exact_optimal", "exact_feasible", "optimal_global_phase"]
                    printstyled("  ","Minimum number of $gate_type gates: ", round(results["objective"], digits = 6),"\n"; color = _main_color)
                
                elseif data["decomposition_type"] == "approximate"
                    printstyled("  ","Minimum number of $gate_type gates: ", round((results["objective"] - 
                        results["objective_slack_penalty"]*LA.norm(results["solution"]["slack_var"])^2), digits = 6),"\n"; color = _main_color)
                end
            end

        end

        printstyled("  ","Optimizer run time      : ", ceil(results["solve_time"], digits=2)," sec.","\n"; color = _main_color)
            
        if results["termination_status"] == MOI.TIME_LIMIT
            printstyled("  ","Termination status : TIME_LIMIT", "\n"; color = _main_color)
        end

        printstyled("=============================================================================","\n"; color = _main_color)      

    else
        Memento.warn(_LOGGER, "Valid integral feasible solutions could not be found to visualize the solution")
    end

    if gate_sequence
        return gates_sol
    end

end

function print_elementary_gates(
    gates::Vector{String}; 
    color=:White
    )
    # Parse gates and group by type
    gate_map = Dict{String, Set{Vector{Int}}}()
    num_qubits = 1
    
    for gate in gates
        name, qubits = QCO._parse_gate_string(gate; type=true, qubits=true)
        name == "Identity" && continue
        
        push!(get!(gate_map, name, Set{Vector{Int}}()), qubits)
        if !isempty(qubits)
            num_qubits = max(num_qubits, maximum(qubits))
        end
    end
    
    # Header
    maxlen = maximum(length.(keys(gate_map)))
    printstyled("  ", lpad("Input gates", maxlen) * " : " * "Qubits\n"; color = color, bold = false)
    
    for name in sort(collect(keys(gate_map)))
        qubits_str = if name == "GR"
            num_qubits <= 3 ? join(1:num_qubits, ",") : "1,...,$num_qubits"
        else
            formatted = [isempty(q) ? "–" : 
                         length(q) == 1 ? string(q[1]) : 
                         "{" * join(q, ",") * "}" for q in gate_map[name]]
            join(unique(sort(formatted)), ",")
        end
        
        printstyled("      ", lpad(name, maxlen) * " : " * qubits_str * "\n"; color=color)
    end
end

"""
    print_circuit(circuit_layers::Vector{Vector{Gate}}, num_qubits::Int)

Print a quantum circuit in a human-readable format, showing each layer of gates.
"""
function print_circuit(
    circuit_layers::Vector{Vector{Gate}}, 
    num_qubits::Int,
    _main_color::Symbol = :White
    )    
    if isempty(circuit_layers)
        printstyled("  Empty circuit\n"; color = _main_color)
        return
    end     
    for (d, layer) in enumerate(circuit_layers)
        printstyled("  ", "depth $d : ", QCO.kron_layer(layer, num_qubits), "\n"; color = _main_color)
    end
end

_gate_id_sequence(z_val::Matrix{<:Number}, maximum_depth::Int64) = 
[findall(isone.(round.(abs.(z_val[:,d]), digits=3)))[1] for d = 1:maximum_depth]

"""
    kron_layer(layer::Vector{Gate}, num_qubits::Int)::String

Converts a layer of quantum gates into a human-readable string representation.
The function formats the layer as a tensor product (⊗) of gates, with identity gates (I)
inserted for qubits that don't have an explicit gate in the layer.
"""
function kron_layer(
    layer::Vector{Gate}, 
    num_qubits::Int
    )::String
    # Sort gates by leftmost qubit (left→right)
    sorted_gates = sort(layer; by = g -> minimum(g.qubits))
    factors      = String[]
    
    qubit_idx = 1
    gate_idx  = 1
    
    while qubit_idx ≤ num_qubits
        if gate_idx > length(sorted_gates) ||
           minimum(sorted_gates[gate_idx].qubits) > qubit_idx

            push!(factors, "I")
            qubit_idx += 1
            continue
        end

        # Gate begins at this wire
        current_gate   = sorted_gates[gate_idx]
        qubit_indices  = collect(current_gate.qubits)

        if length(qubit_indices) == 1
            push!(factors, current_gate.label)
        else
            push!(factors, string(current_gate.label, "_{", join(qubit_indices, ","), "}"))
        end

        lastq = maximum(qubit_indices)   # largest wire index touched
        qubit_idx = lastq + 1
        gate_idx  += 1
    end
    
    return join(factors, " ⊗ ")
end

"""
    circuit_unitary(layers::Vector{Vector{Gate}})

Returns the overall unitary matrix representing a quantum circuit by multiplying
the unitary matrices of each layer (left to right).
"""
function circuit_unitary(layers::Vector{Vector{Gate}})
    dim = size(layers[1][1].matrix,1)        # 2^n
    U = Matrix{ComplexF64}(LA.I, dim, dim)
    for layer in layers
        G = Matrix{ComplexF64}(LA.I, dim, dim)
        for g in layer
            G = g.matrix * G 
        end
        U *= G
    end
    return U
end

function build_circuit_layers(
    id_seq::Vector{Int}, 
    gates_dict::Dict{String,Any}
    )

    gates = Gate[]
    for id in id_seq
        info   = gates_dict[string(id)]
        ("Identity" in info["type"]) && continue
        typestr = info["type"][1]
        
        head, params = occursin('(', typestr) ? split(typestr,'(',limit=2) : (typestr,"")
        parts = split(head, '_')
        label = parts[1]
        
        qubits = Vector{Int}(info["qubit_loc"])
        matrix = convert(Array{ComplexF64,2}, info["matrix"])

        # round angles if present
        if haskey(info,"angle")
            angle = info["angle"]
            order = ["θ", "ϕ", "λ"]
            vals = [round(rad2deg(angle[k]), digits=3) for k in order if haskey(angle, k)]
            label *= "(" * join(vals,",") * ")"
        elseif params != ""
            label *= "(" * params
        end

        push!(gates, Gate(label, qubits, matrix))
    end

    num_qubits = round(Int, log2(size(gates[1].matrix,1)))
    layers_before = [[g] for g in gates]
    layers_after  = QCO.compress_circuit(gates, num_qubits)
    return layers_before, layers_after
end

"""
    compress_circuit(gates::Vector{Gate}, num_qubits::Int)

This function compresses a quantum circuit by placing gates into the earliest possible layer where 
they can be executed in parallel. Gates can be placed in the same layer if they operate on disjoint 
sets of qubits. If no suitable layer is found, a new layer is created. This function reduces 
circuit depth while preserving the total number of gates and the logical operation.
"""
function compress_circuit(gates::Vector{Gate}, num_qubits::Int)
    layers, masks = Vector{Vector{Gate}}(), Vector{Vector{Int}}()
    full = collect(1:num_qubits)

    for g in gates
        (all(q -> q in full, g.qubits)) || Memento.error(_LOGGER, "Gate $(g.label) uses invalid qubit (n=$num_qubits)")

        pos = length(layers)+1
        for l = length(layers):-1:1
            if !isempty(intersect(g.qubits, masks[l])); break; end
            pos = l
        end
        if pos > length(layers)
            push!(layers, [g]); push!(masks, copy(g.qubits))
        else
            push!(layers[pos], g); append!(masks[pos], setdiff(g.qubits, masks[pos]))
        end
    end
    return layers
end

function validate_circuit(
    data::Dict{String,Any},
    layers::Vector{Vector{Gate}};
    atol = 1e-4,
    error_message = true,
    global_phase = false,
    circuit_type = nothing
    )

    U_sol = QCO.circuit_unitary(layers)

    tgt = data["are_gates_real"] ?
            real(data["target_gate"]) :
            QCO.real_to_complex_gate(data["target_gate"])

    target_gate = convert(Array{ComplexF64,2}, tgt)

    valid = if (data["decomposition_type"] in ["exact_optimal", "exact_feasible"]) 
        QCO.isapprox(U_sol, target_gate; atol = atol)
    elseif (data["decomposition_type"] == "optimal_global_phase") || global_phase
        QCO.isapprox_global_phase(U_sol, target_gate; tol_0 = atol)
    else
        false
    end

    # Report error if validation fails
    if !valid && error_message
        msg = "Decomposition is not valid: Problem may be infeasible"
        circuit_type !== nothing && (msg = "Decomposition is not valid: Circuit may be infeasible for $circuit_type")
        Memento.error(_LOGGER, msg)
    end
    
    return valid
end