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
    @show id_sequence
    
    cirq_precompress, cirq_postcompress = QCO.build_circuit_layers(id_sequence, data["gates_dict"])

    QCO.print_circuit(cirq_precompress, data["num_qubits"])

    QCO.print_circuit(cirq_postcompress, data["num_qubits"])

    cirq_postcompress_depth = length(cirq_postcompress)

    @show length(cirq_precompress), length(cirq_postcompress)
    
    (data["decomposition_type"] in ["exact_optimal", "exact_feasible", "optimal_global_phase"]) && 
                                        QCO.validate_circuit(data, cirq_precompress, global_phase = false) && 
                                        QCO.validate_circuit(data, cirq_postcompress, global_phase = false)

    # gates_sol, gates_sol_compressed = QCO.get_postprocessed_circuit(results, data)
    
    if !isempty(cirq_postcompress)

        printstyled("\n","=============================================================================","\n"; color = _main_color)
        printstyled("QuantumCircuitOpt version: ", Pkg.TOML.parse(read(string(pkgdir(QCO), "/Project.toml"), String))["version"], "\n"; color = _header_color, bold = true)

        printstyled("\n","Quantum Circuit Model Data"; color = _header_color, bold = true)
        
        printstyled("\n","  ","Number of qubits: ", data["num_qubits"], "\n"; color = _main_color)
        
        printstyled("  ","Total number of elementary gates (after presolve): ",size(data["gates_real"])[3],"\n"; color = _main_color)
        
        printstyled("  ","Maximum depth of decomposition: ", data["maximum_depth"],"\n"; color = _main_color)
        
        printstyled("  ","Elementary gates: ", data["elementary_gates"],"\n"; color = _main_color)

        if "discretization" in keys(data)
            for i in keys(data["discretization"])
                printstyled("    ","$i discretization: ", ceil.(rad2deg.(data["discretization"][i]), digits = 1),"\n"; color = _main_color)
            end
        end
                
        printstyled("  ","Type of decomposition: ", data["decomposition_type"],"\n"; color = _main_color)

        printstyled("  ","MIP optimizer: ", results["optimizer"],"\n"; color = _main_color)

        printstyled("\n","Optimal Circuit Decomposition","\n"; color = _header_color, bold = true)

        QCO.print_circuit(cirq_postcompress, data["num_qubits"])
        
        # for i=1:cirq_postcompress_depth
        #     if i != cirq_postcompress_depth
        #         printstyled(gates_sol_compressed[i], " * "; color = _main_color)
        #     else    
        #         if data["decomposition_type"] in ["exact_optimal", "exact_feasible", "optimal_global_phase"]
        #             printstyled(gates_sol_compressed[i], " = ", "Target gate","\n"; color = _main_color)
        #         elseif data["decomposition_type"] == "approximate"
        #             printstyled(gates_sol_compressed[i], " ≈ ", "Target gate","\n"; color = _main_color)
        #         end
        #     end
        # end

        if data["decomposition_type"] == "approximate"
            printstyled("  ","||Decomposition error||₂: ", round(LA.norm(results["solution"]["slack_var"]), digits = 10),"\n"; color = _main_color)
        end

        if data["objective"] == "minimize_depth"

            if length(data["identity_idx"]) >= 1 && (data["decomposition_type"] !== "exact_feasible") && !(results["termination_status"] == MOI.TIME_LIMIT)
                printstyled("  ","Minimum optimal depth: ", cirq_postcompress_depth,"\n"; color = _main_color)
            else 
                printstyled("  ","Decomposition depth: ", cirq_postcompress_depth,"\n"; color = _main_color)
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

        printstyled("  ","Optimizer run time: ", ceil(results["solve_time"], digits=2)," sec.","\n"; color = _main_color)
            
        if results["termination_status"] == MOI.TIME_LIMIT
            printstyled("  ","Termination status: TIME_LIMIT", "\n"; color = _main_color)
        end

        printstyled("=============================================================================","\n"; color = _main_color)      

    else
        Memento.warn(_LOGGER, "Valid integral feasible solutions could not be found to visualize the solution")
    end

    if gate_sequence
        return gates_sol
    end

end

"""
    print_circuit(circuit_layers::Vector{Vector{Gate}}, num_qubits::Int)

Print a quantum circuit in a human-readable format, showing each layer of gates.
"""
function print_circuit(
    circuit_layers::Vector{Vector{Gate}}, 
    num_qubits::Int
    )    
    if isempty(circuit_layers)
        printstyled("  Empty circuit\n"; color = _main_color)
        return
    end    
    for (d, layer) in enumerate(circuit_layers)
        printstyled("  ", "depth $d : ", QCO.kron_layer(layer, num_qubits), "\n"; color = _main_color)
    end
end


function get_postprocessed_circuit_old(results::Dict{String, Any}, data::Dict{String, Any})

    gates_sol = Array{String,1}()
    id_sequence = QCO._gate_id_sequence(results["solution"]["z_bin_var"], data["maximum_depth"])
    (data["decomposition_type"] in ["exact_optimal", "exact_feasible", "optimal_global_phase"]) && QCO.validate_circuit(data, id_sequence)

    for d = 1:data["maximum_depth"]
        
        gate_id = data["gates_dict"]["$(id_sequence[d])"]

        if !("Identity" in gate_id["type"])
            
            s1 = gate_id["type"][1]

            if occursin(kron_symbol, s1)
                push!(gates_sol, s1)
            elseif !(QCO._parse_gate_string(s1, type = true) in union(QCO.ONE_QUBIT_GATES_ANGLE_PARAMETERS, QCO.TWO_QUBIT_GATES_ANGLE_PARAMETERS, QCO.MULTI_QUBIT_GATES_ANGLE_PARAMETERS))
                push!(gates_sol, s1) 
            else
                if "angle" in keys(gate_id)

                    if length(keys(gate_id["angle"])) == 1 
                        θ = round(rad2deg(gate_id["angle"]["θ"]), digits = 3)
                        s3 = "$(θ)"
                        push!(gates_sol, string(s1,"(", s3, ")"))

                    elseif length(keys(gate_id["angle"])) == 2
                        θ = round(rad2deg(gate_id["angle"]["θ"]), digits = 3)
                        ϕ = round(rad2deg(gate_id["angle"]["ϕ"]), digits = 3)
                        s3 = string("(","$(θ)",",","$(ϕ)",")")
                        push!(gates_sol, string(s1, s3))

                    elseif length(keys(gate_id["angle"])) == 3
                        θ = round(rad2deg(gate_id["angle"]["θ"]), digits = 3)
                        ϕ = round(rad2deg(gate_id["angle"]["ϕ"]), digits = 3)
                        λ = round(rad2deg(gate_id["angle"]["λ"]), digits = 3)
                        s3 = string("(","$(θ)",",","$(ϕ)", ",","$(λ)",")")
                        push!(gates_sol, string(s1, s3))
                    end
                end
            end
        end
    end

    @show gates_sol

    gates_sol_compressed = QCO.get_depth_compressed_circuit(data["num_qubits"], gates_sol)

    return gates_sol, gates_sol_compressed
end

"""
    validate_circuit(data::Dict{String, Any}, id_sequence::Array{Int64,1})

This function validates the circuit decomposition if it is indeed exact with respect to the specified target gate. 
"""
function validate_circuit_old(
    data::Dict{String, Any}, 
    id_sequence::Array{Int64,1}; 
    error_message = true
    )
    # Multiply all gates in sequence
    M_sol = reduce(*, [data["gates_dict"]["$i"]["matrix"] for i in id_sequence], 
                  init=Array{Complex{Float64},2}(Matrix(LA.I, 2^(data["num_qubits"]), 2^(data["num_qubits"]))))
    
    target_gate = data["are_gates_real"] ? real(data["target_gate"]) : 
                                          QCO.real_to_complex_gate(data["target_gate"])
    target_gate = convert(Array{Complex{Float64},2}, target_gate)
    
    # Check if solution matches target
    valid_status = if data["decomposition_type"] in ["exact_optimal", "exact_feasible"]
        QCO.isapprox(M_sol, target_gate, atol = 1E-4)
    elseif data["decomposition_type"] in ["optimal_global_phase"]
        QCO.isapprox_global_phase(M_sol, target_gate)
    else
        false
    end

    (!valid_status && error_message) && Memento.error(_LOGGER, "Decomposition is not valid: Problem may be infeasible")
    
    return valid_status
end

"""
    get_depth_compressed_circuit(num_qubits::Int64, gates_sol::Array{String,1})

Given the number of qubits and the sequence of gates from the solution, this function returns a 
decomposition of gates after compressing adjacent pair of gates represented on two separate qubits. 
For example, gates H1 and H2 appearing in a sequence will be compressed to H1xH2 (kron(H1,H2)). 
This functionality is currently supported only for two qubit circuits and gates without angle parameters. 
"""
function get_depth_compressed_circuit(num_qubits::Int64, gates_sol::Array{String,1})
    # This part of the code may be hacky. This needs to be updated once the input format gets cleaned up for elementary gates with U and R gates.     

    if (length(gates_sol) == 1) || (num_qubits > 2)
        return gates_sol
    end
    
    gates_sol_compressed = String[]

    angle_param_gate = false
    for i=1:length(gates_sol)
        if !occursin(kron_symbol, gates_sol[i])
            gates_sol_type = !occursin("GR", gates_sol[i]) ? QCO._parse_gate_string(gates_sol[i], type = true) : "GR"
            if gates_sol_type in union(QCO.ONE_QUBIT_GATES_ANGLE_PARAMETERS, QCO.TWO_QUBIT_GATES_ANGLE_PARAMETERS, QCO.MULTI_QUBIT_GATES_ANGLE_PARAMETERS)
                angle_param_gate = true
                break
            end
        end       
    end

    if angle_param_gate
        return gates_sol
    end
    
    i = 1
    while i <= length(gates_sol)
        if i < length(gates_sol) && 
           !QCO.is_multi_qubit_gate(gates_sol[i]) && !QCO.is_multi_qubit_gate(gates_sol[i+1]) &&
           ((occursin('1', gates_sol[i]) && occursin('2', gates_sol[i+1])) || 
            (occursin('2', gates_sol[i]) && occursin('1', gates_sol[i+1])))
            
            # Determine order for kronecker product
            gate_string = occursin('1', gates_sol[i]) ? 
                          string(gates_sol[i], "x", gates_sol[i+1]) : 
                          string(gates_sol[i+1], "x", gates_sol[i])
            
            push!(gates_sol_compressed, gate_string)
            i += 2  # Skip the next gate as we've already processed it
        else
            push!(gates_sol_compressed, gates_sol[i])
            i += 1
        end
    end

    if isempty(gates_sol_compressed)
        Memento.error(_LOGGER, "Compressed gates solution is empty")
    end

    return gates_sol_compressed
end

_gate_id_sequence(z_val::Matrix{<:Number}, maximum_depth::Int64) = 
[findall(isone.(round.(abs.(z_val[:,d]), digits=3)))[1] for d = 1:maximum_depth]

# function kron_layer(
#     layer::Vector{Gate}, 
#     n::Int
#     )::String

#     sorted = sort(layer; by = g -> minimum(g.qubits))
#     facs, qptr = String[], 1 # tensor factors, next qubit index

#     for g in sorted
#         firstq, lastq = minimum(g.qubits), maximum(g.qubits)
#         for _ in qptr:firstq-1
#             push!(facs, "I")
#         end
#         if length(g.qubits) == 1
#             push!(facs, g.label)
#         else
#             sub  = join(sort(collect(g.qubits)), ",")      
#             push!(facs, string(g.label, "_{", sub, "}"))   # CNot_{2,3}
#         end
#         qptr = lastq + 1
#     end
#     for _ in qptr:n
#         push!(facs,"I")
#     end

#     return join(facs, " ⊗ ")
# end
# function kron_layer(
#     layer::Vector{Gate}, 
#     n::Int
#     )::String

#     # one factor per qubit, qubit 1 on the left
#     ops = fill("I", n)                          

#     # union of all multi-qubit gate names known to QCO
#     multiset = union(QCO.TWO_QUBIT_GATES, QCO.MULTI_QUBIT_GATES)

#     for g in layer
#         qs   = sort(collect(g.qubits))          # e.g. [1,2]
#         base = match(r"^[A-Za-z]+", g.label).match   # strip "_…" or "(…"
        
#         if length(qs) == 1                          # 1-qubit gate
#             ops[qs[1]] = g.label

#         elseif base in multiset                     # true multi-qubit gate
#             ops[qs[1]] = string(base, "_{", join(qs, ","), "}")
#             # leave remaining qubits as "I"

#         else                                        # identical 1-qubit gates
#             for q in qs
#                 ops[q] = g.label
#             end
#         end
#     end

#     return join(ops, " ⊗ ")
# end


function kron_layer(
    layer::Vector{Gate}, 
    n::Int
    )::String

    # sort by left-most qubit so we can scan left → right
    gates   = sort(layer; by = g -> minimum(g.qubits))
    facs    = String[]          # tensor factors we’ll return
    idxgate = 1                 # position inside `gates`
    q       = 1                 # current wire we are printing

    while q ≤ n
        if idxgate > length(gates) || minimum(gates[idxgate].qubits) > q
            # no gate starts on this wire → identity here
            push!(facs, "I")
            q += 1
            continue
        end

        g  = gates[idxgate]
        qs = sort(collect(g.qubits))
        firstq, lastq = qs[1], qs[end]

        if length(qs) == 1
            push!(facs, g.label)            # ordinary 1-qubit gate
        else
            push!(facs, string(g.label, "_{", join(qs, ","), "}"))
        end

        q = lastq + 1                       # skip wires the gate already covers
        idxgate += 1
    end

    return join(facs, " ⊗ ")
end


function circuit_unitary(layers::Vector{Vector{Gate}})
    dim = size(layers[1][1].mat,1)        # 2^n
    U = Matrix{ComplexF64}(LA.I, dim, dim)
    for layer in layers
        G = Matrix{ComplexF64}(LA.I, dim, dim)
        for g in layer
            G = g.mat * G                 # gates in same layer commute (disjoint qubits)
        end
        U = G * U
    end
    U
end

read_int_prefix(s::AbstractString) = parse(Int, match(r"^\d+", s).match)

function build_circuit_layers(
    id_seq::Vector{Int}, 
    gates_dict::Dict{String,Any}
    )

    gates = Gate[]
    for id in id_seq
        info   = gates_dict[string(id)]
        typestr = info["type"][1]
        typestr == "Identity" && continue          # skip

        head, params = occursin('(', typestr) ? split(typestr,'(',limit=2) : (typestr,"")
        parts = split(head, '_'); op = parts[1]
        qs    = BitSet(QCO.read_int_prefix.(parts[2:end]))
        mat   = convert(Array{ComplexF64,2}, info["matrix"])

        # round angles if present
        label = op
        if haskey(info,"angle")
            a = info["angle"]; order = ["θ","ϕ","λ"]
            vals = [round(rad2deg(a[k]), digits=3) for k in order if haskey(a,k)]
            label *= "(" * join(vals,",") * ")"
        elseif params!=""; label *= "(" * params; end

        push!(gates, Gate(label, qs, mat))
    end

    num_qubits = round(Int, log2(size(gates[1].mat,1)))
    layers_before = [[g] for g in gates]
    layers_after  = QCO.compress_circuit(gates, num_qubits)
    return layers_before, layers_after
end

function compress_circuit(gates::Vector{Gate}, num_qubits::Int)
    layers, masks = Vector{Vector{Gate}}(), BitSet[]
    full = BitSet(1:num_qubits)

    for g in gates
        (g.qubits ⊆ full) || Memento.error(_LOGGER, "Gate $(g.label) uses invalid qubit (n=$num_qubits)")

        pos = length(layers)+1
        for l = length(layers):-1:1
            if !isempty(intersect(g.qubits, masks[l])); break; end
            pos = l
        end
        if pos > length(layers)
            push!(layers, [g]); push!(masks, copy(g.qubits))
        else
            push!(layers[pos], g); union!(masks[pos], g.qubits)
        end
    end
    return layers
end

function validate_circuit(
    data::Dict{String,Any},
    layers::Vector{Vector{Gate}};
    atol = 1e-4,
    error_message = true,
    global_phase = false
    )

    U_sol = QCO.circuit_unitary(layers)

    tgt = data["are_gates_real"] ?
            real(data["target_gate"]) :
            QCO.real_to_complex_gate(data["target_gate"])

    target_gate = convert(Array{ComplexF64,2}, tgt)

    valid = if (data["decomposition_type"] in ["exact_optimal", "exact_feasible"]) 
        QCO.isapprox(U_sol, target_gate; atol = atol)
    elseif (data["decomposition_type"] == "optimal_global_phase") || global_phase
        QCO.isapprox_global_phase(U_sol, target_gate; atol = atol)
    else
        false
    end

    if !valid && error_message
        Memento.error(_LOGGER, "Decomposition is not valid: Problem may be infeasible")
    end
    return valid
end

function _parse_gate(s::AbstractString)::Gate
    head, params = occursin('(', s) ? split(s, '(', limit = 2) : (s,"")
    parts = split(head, '_')
    op    = parts[1]
    qs    = BitSet(parse.(Int, parts[2:end]))
    label = params=="" ? op : op * "(" * params
    Gate(label, qs, Matrix{ComplexF64}(I,1,1))   # matrix overwritten later
end

_explode_kron(seq::Vector{String}) = [QCO._parse_gate(term) for item in seq for term in split(item,'x')]