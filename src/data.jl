import LinearAlgebra: I

"""
    get_data(params::Dict{String, Any}; eliminate_identical_gates = true)

Given the user input `params` dictionary, this function returns a dictionary of processed data which contains all the 
necessary information to formulate the optimization model for the circuit design problem. 
"""
function get_data(params::Dict{String, Any}; eliminate_identical_gates = true)
    
    # Number of qubits
    num_qubits = haskey(params, "num_qubits") ?
        params["num_qubits"] :
        Memento.error(_LOGGER, "Number of qubits has to be specified by the user")
    num_qubits < 2 && Memento.error(_LOGGER, "Minimum of 2 qubits is necessary")

    # Depth
    maximum_depth = haskey(params, "maximum_depth") ?
        params["maximum_depth"] :
        Memento.error(_LOGGER, "Maximum depth of the decomposition has to be specified by the user")
    maximum_depth < 2 && Memento.error(_LOGGER, "Minimum depth of 2 is necessary")

    # Elementary gates
    (!haskey(params, "elementary_gates") || isempty(params["elementary_gates"])) &&
        Memento.error(_LOGGER, "Specify at least two unique unitary elementary gates")

    # Input circuit
    input_circuit = get(params, "input_circuit", [])
    input_circuit_dict = Dict{String, Any}()

    if !isempty(input_circuit)
        if length(input_circuit) <= maximum_depth
            input_circuit_dict = QCO.get_input_circuit_dict(input_circuit, params)
        else
            Memento.warn(_LOGGER, "Neglecting the input circuit: its depth exceeds the allowable depth")
        end
    end

    # Decomposition type
    decomposition_type = get(params, "decomposition_type", "exact_optimal")
    !(decomposition_type in ["exact_optimal", "exact_feasible", "optimal_global_phase", "approximate"]) &&
        Memento.error(_LOGGER, "Decomposition type not supported")

    # Objective function
    objective = get(params, "objective", "minimize_depth")

    elementary_gates = unique(params["elementary_gates"])
    
    if length(elementary_gates) < length(params["elementary_gates"])
        Memento.warn(_LOGGER, "Eliminating non-unique gates in the input elementary gates")
    end

    gates_dict, are_elementary_gates_real = QCO.get_elementary_gates_dictionary(params, elementary_gates)

    target_real, is_target_real = QCO.get_target_gate(params, are_elementary_gates_real, decomposition_type)

    gates_dict_unique, M_real_unique, identity_idx, cnot_idx, T_idx = 
            QCO.eliminate_nonunique_gates(gates_dict, are_elementary_gates_real, eliminate_identical_gates = eliminate_identical_gates)

    # Pick "Identity" if user didn't specify
    initial_gate_str = get(params, "initial_gate", "Identity")
    initial_gate_str == "Identity" ||
        Memento.error(_LOGGER, "Currently, only \"Identity\" is supported as an initial gate")

    # Build the initial gate
    if are_elementary_gates_real && is_target_real
        initial_gate = real(QCO.IGate(num_qubits))
    else
        initial_gate = QCO.complex_to_real_gate(QCO.IGate(num_qubits))
    end
    
    data = Dict{String, Any}("num_qubits" => num_qubits,
                             "maximum_depth" => maximum_depth,
                             "gates_dict" => gates_dict_unique,
                             "gates_real" => M_real_unique,
                             "initial_gate" => initial_gate,
                             "identity_idx" => identity_idx,
                             "cnot_idx" => cnot_idx,
                             "T_idx" => T_idx,
                             "elementary_gates" => elementary_gates,
                             "target_gate" => target_real,
                             "are_gates_real" => (are_elementary_gates_real && is_target_real),
                             "objective" => objective,
                             "decomposition_type" => decomposition_type
                             )

    if data["are_gates_real"]
        Memento.info(_LOGGER, "Detected all-real elementary and target gates")
    end
    
    # Rotation and Universal gate angle discretizations
    data = QCO._populate_data_angle_discretization!(data, params)

    # Determinant test for input gates
    (data["decomposition_type"] in ["exact_optimal", "exact_feasible"]) && QCO._determinant_test_for_infeasibility(data)

    # Input circuit
    if length(keys(input_circuit_dict)) > 0
        data["input_circuit"] = input_circuit_dict
    end

    # CNOT lower/upper bound
    data = QCO._get_cnot_bounds!(data, params)
                         
    return data
end

"""
    eliminate_nonunique_gates(gates_dict::Dict{String, Any})

"""
function eliminate_nonunique_gates(gates_dict::Dict{String, Any}, are_elementary_gates_real::Bool; eliminate_identical_gates = false)

    num_gates = length(keys(gates_dict))

    # This identifies if all the elementary gates have only real entries and returns compact matrices
    if are_elementary_gates_real
        M_real = zeros(size(gates_dict["1"]["matrix"])[1], size(gates_dict["1"]["matrix"])[2], num_gates)
    else
        M_real = zeros(2*size(gates_dict["1"]["matrix"])[1], 2*size(gates_dict["1"]["matrix"])[2], num_gates)
    end
    
    for i=1:num_gates
        if are_elementary_gates_real
            M_real[:,:,i] = real(gates_dict["$i"]["matrix"])
        else 
            M_real[:,:,i] = QCO.complex_to_real_gate(gates_dict["$i"]["matrix"])
        end
    end

    M_real_unique = M_real
    M_real_idx = collect(1:size(M_real)[3]) 

    if eliminate_identical_gates
        M_real_unique, M_real_idx = QCO.unique_matrices(M_real)
    end
    
    gates_dict_unique = Dict{String, Any}()

    if size(M_real_unique)[3] < size(M_real)[3]

        Memento.info(_LOGGER, "Detected $(size(M_real)[3]-size(M_real_unique)[3]) non-unique gates (after discretization)")

        for i = 1:length(M_real_idx)
            gates_dict_unique["$i"] = gates_dict["$(M_real_idx[i])"]
        end
    
    else
        gates_dict_unique = gates_dict
    end

    identity_idx = QCO._get_identity_idx(M_real_unique)

    for i_id = 1:length(identity_idx)
        if !("Identity" in gates_dict_unique["$(identity_idx[i_id])"]["type"])
            push!(gates_dict_unique["$(identity_idx[i_id])"]["type"], "Identity")
        end
    end

    cnot_idx = QCO._get_cnot_idx(gates_dict_unique)
    T_idx = QCO._get_T_idx(gates_dict_unique)
    
    return gates_dict_unique, M_real_unique, identity_idx, cnot_idx, T_idx
end

function _populate_data_angle_discretization!(data::Dict{String, Any}, params::Dict{String, Any})

    one_angle_gates, two_angle_gates, three_angle_gates = QCO._get_angle_gates_idx(data["elementary_gates"])

    if !isempty(union(one_angle_gates, two_angle_gates, three_angle_gates)) 
        if length(findall(x -> occursin("discretization", x), collect(keys(params)))) == 0
            Memento.error(_LOGGER, "Enter valid discretization angles in the input params")
        end

        data["discretization"] = Dict{String, Any}()
        angle_types = [one_angle_gates   => ["_discretization"],
                       two_angle_gates   => ["_θ_discretization", "_ϕ_discretization"],
                       three_angle_gates => ["_θ_discretization", "_ϕ_discretization", "_λ_discretization"]]

        for (gates, suffixes) in angle_types
            for i in gates
                gate_type = QCO._parse_gate_string(data["elementary_gates"][i], type = true)
                for suffix in suffixes
                    key = "$(gate_type)$(suffix)"
                    data["discretization"][key] = Float64.(params[key])
                end
            end
        end

    end

    return data
end

"""
    get_target_gate(params::Dict{String, Any}, are_elementary_gates_real::Bool)

Given the user input `params` dictionary and a boolean if all the input elementary gates are real, 
this function returns the corresponding real version of the target gate. 
""" 
function get_target_gate(params::Dict{String, Any}, are_elementary_gates_real::Bool, decomposition_type::String)

    if !("target_gate" in keys(params)) || isempty(params["target_gate"])
        Memento.error(_LOGGER, "Target gate not found in the input data")
    end 

    if (size(params["target_gate"])[1] != size(params["target_gate"])[2]) || (size(params["target_gate"])[1] != 2^params["num_qubits"])
        Memento.error(_LOGGER, "Dimensions of target gate do not match the input num_qubits")
    end
    
    # Identify if the target gate has only real entries or if it is real up to a global phase and returns a compact matrix
    is_target_real = QCO.is_gate_real(params["target_gate"])

    if are_elementary_gates_real

        global_phase = 0

        if (decomposition_type in ["optimal_global_phase"]) && !is_target_real
            ref_nonzero_r, ref_nonzero_c = QCO._get_nonzero_idx_of_complex_matrix(params["target_gate"])
            global_phase = angle(params["target_gate"][ref_nonzero_r, ref_nonzero_c])
        end
        
        is_target_real_up_to_phase = QCO.is_gate_real(exp(-im*global_phase)*params["target_gate"])
            
        if is_target_real_up_to_phase
            return real(exp(-im*global_phase)*params["target_gate"]), is_target_real_up_to_phase
        else    
            Memento.error(_LOGGER, "Infeasible decomposition: all elementary gates have zero imaginary parts and target is not real for exact decomposition or not real up to a global phase for optimal_global_phase decomposition.")
        end
        
    else
        return QCO.complex_to_real_gate(params["target_gate"]), is_target_real
    end

end

function get_elementary_gates_dictionary(params::Dict{String, Any}, elementary_gates::Array{String,1})

    num_qubits = params["num_qubits"]

    one_angle_gates, two_angle_gates, three_angle_gates = QCO._get_angle_gates_idx(elementary_gates)
    angle_gates = [one_angle_gates, two_angle_gates, three_angle_gates]
    num_angle_gates = [1, 2, 3]
    angle_gates_dicts = [Dict{String, Any}() for _ in 1:length(num_angle_gates)]

    for (i, gates) in enumerate(angle_gates)
        if !isempty(gates)
            angle_gates_dicts[i] = QCO.get_all_angle_gates(params, elementary_gates, gates, num_angle_gates[i])
        end
    end
    all_angle_gates_dict = merge(angle_gates_dicts...)

    gates_dict = Dict{String, Any}()
    counter = 1
    angle_keys = Dict(1 => ["θ"], 2 => ["θ", "ϕ"], 3 => ["θ", "ϕ", "λ"])

    for i=1:length(elementary_gates)

        if i in union(one_angle_gates, two_angle_gates, three_angle_gates)       
            M_elementary_dict = all_angle_gates_dict[elementary_gates[i]]
            num_angles = i in one_angle_gates ? 1 : i in two_angle_gates ? 2 : i in three_angle_gates ? 3 : 0
            
            for k in keys(M_elementary_dict) # Angle
                for l in keys(M_elementary_dict[k]["$(num_qubits)qubit_rep"]) # qubits (which will now be 1)

                    gates_dict["$counter"] = Dict{String, Any}("type"      => [elementary_gates[i]],
                                                               "angle"     => Dict(key => M_elementary_dict[k][key] for key in angle_keys[num_angles]),
                                                               "qubit_loc" => l,
                                                               "matrix"    => M_elementary_dict[k]["$(num_qubits)qubit_rep"][l])
                    counter += 1
                end
            end
        
        else 
            M = QCO.unitary(elementary_gates[i], num_qubits)
            gates_dict["$counter"] = Dict{String, Any}("type"   => [elementary_gates[i]],
                                                       "matrix" => M)
            counter += 1
        end

    end

    are_elementary_gates_real = all(QCO.is_gate_real(gate["matrix"]) for gate in values(gates_dict))
    return gates_dict, are_elementary_gates_real
end

function get_all_angle_gates(params::Dict{String, Any}, 
    elementary_gates::Array{String,1}, 
    gates_idx::Vector{Int64}, 
    num_angles::Int)

    gates_complex = Dict{String, Any}()
    angles = ["θ", "ϕ", "λ"][1:num_angles]

    for i in gates_idx
        input_gate = elementary_gates[i]
        gate_name = QCO._parse_gate_string(input_gate, type=true)

        disc_angles = [Float64.(params["$(gate_name)_$(num_angles == 1 ? "" : "$(angle)_")discretization"]) for angle in angles]

        if any(isempty, disc_angles)
            Memento.error(_LOGGER, "Empty discretization angle for $input_gate gate. Input a valid angle")
        end

        gates_complex[input_gate] = QCO.get_discretized_angle_gates(input_gate, Dict{String, Any}(), disc_angles, params["num_qubits"])
    end

    return gates_complex
end

function get_discretized_angle_gates(gate_type::String, 
                                     M::Dict{String, Any}, 
                                     discretizations::Array{Array{Float64,1},1}, 
                                     num_qubits::Int64)
  
    counter = 1
    product_iterator = Base.Iterators.product(reverse(discretizations)...)

    for p_iter in product_iterator
        angles = collect(reverse(p_iter))

        angle_dict = Dict{String, Any}(["θ", "ϕ", "λ"][1:length(angles)] .=> angles)
        angle_dict["$(num_qubits)qubit_rep"] = Dict{String, Any}()

        qubit_loc = QCO._parse_gate_string(gate_type, qubits=true)
        qubit_loc_str = join(qubit_loc, QCO.qubit_separator)

        if !(gate_type in QCO.MULTI_QUBIT_GATES)
            angle_dict["$(num_qubits)qubit_rep"]["qubit_$(qubit_loc_str)"] = QCO.unitary(gate_type, num_qubits, angle = angles)
        else
            angle_dict["$(num_qubits)qubit_rep"]["multi_qubits"] = QCO.unitary(gate_type, num_qubits, angle = angles)
        end

        M["angle_$(counter)"] = angle_dict
        counter += 1
    end

    return M
end

"""
    unitary(input::String, num_qubits::Int64; angle = nothing)

Given an input string representing the gate and number of qubits of the circuit, this function returns a full-sized 
gate with respect to the input number of qubits. For example, if `num_qubits = 3` and the input gate in `H_3` 
(Hadamard on third qubit), then this function returns `IGate ⨂ IGate ⨂ HGate`, where IGate and HGate are single 
qubit Identity and Hadamard gates, respectively. Note that `angle` vector is an optional input which is 
necessary when the input gate is parametrized by Euler angles.
"""
function unitary(
    input::String, 
    num_qubits::Int64; 
    angle = nothing
    )

    if num_qubits > 10
        Memento.error(_LOGGER, "Greater than 10 qubits is currently not supported")
    end

    if occursin(QCO.kron_symbol, input)
        return QCO.get_full_sized_kron_gate(input, num_qubits)
    end

    if input == "Identity"
        return QCO.IGate(num_qubits)
    end

    gate_type, qubit_loc = QCO._parse_gate_string(input, type = true, qubits = true)

    if !(gate_type in union(QCO.ONE_QUBIT_GATES, QCO.TWO_QUBIT_GATES, QCO.MULTI_QUBIT_GATES))
        Memento.error(_LOGGER, "Specified $input gate does not exist in the predefined set of gates")
    end

    QCO._catch_input_gate_errors(gate_type, qubit_loc, num_qubits, input; angle = angle)
    
    #----------------------;
    #   One qubit gates    ;
    #----------------------; 
    if length(qubit_loc) == 1 

        if gate_type in QCO.ONE_QUBIT_GATES_CONSTANTS
            
            return QCO.kron_single_qubit_gate(num_qubits, getfield(QCO, Symbol(gate_type, "Gate"))(), "q$(qubit_loc[1])")

        elseif gate_type in QCO.ONE_QUBIT_GATES_ANGLE_PARAMETERS
            
            if (angle !== nothing) && (length(angle) > 0)
                
                if length(angle) == 1 
                    return QCO.kron_single_qubit_gate(num_qubits, getfield(QCO, Symbol(gate_type, "Gate"))(angle[1]), "q$(qubit_loc[1])")
                elseif length(angle) == 2
                    return QCO.kron_single_qubit_gate(num_qubits, getfield(QCO, Symbol(gate_type, "Gate"))(angle[1], angle[2]), "q$(qubit_loc[1])")
                elseif length(angle) == 3
                    return QCO.kron_single_qubit_gate(num_qubits, getfield(QCO, Symbol(gate_type, "Gate"))(angle[1], angle[2], angle[3]), "q$(qubit_loc[1])")
                end

            else 
                Memento.error(_LOGGER, "Enter a valid angle parameter for the input $input gate")
            end

        end
    
    #----------------------;
    #   Two qubit gates    ;
    #----------------------; 
    elseif length(qubit_loc) == 2 

        if gate_type in QCO.TWO_QUBIT_GATES_CONSTANTS

            if (qubit_loc[1] < qubit_loc[2]) || (gate_type in QCO.TWO_QUBIT_GATES_CONSTANTS_SYMMETRIC)
                return QCO.kron_two_qubit_gate(num_qubits, getfield(QCO, Symbol(gate_type, "Gate"))(), "q$(qubit_loc[1])", "q$(qubit_loc[2])")
            else
                return QCO.kron_two_qubit_gate(num_qubits, getfield(QCO, Symbol(gate_type, "RevGate"))(), "q$(qubit_loc[1])", "q$(qubit_loc[2])")
            end

        elseif gate_type in QCO.TWO_QUBIT_GATES_ANGLE_PARAMETERS
            
            if (angle !== nothing) && (length(angle) > 0)
                
                if length(angle) == 1 
                    if (qubit_loc[1] < qubit_loc[2])
                        return QCO.kron_two_qubit_gate(num_qubits, getfield(QCO, Symbol(gate_type, "Gate"))(angle[1]), "q$(qubit_loc[1])", "q$(qubit_loc[2])")
                    else
                        return QCO.kron_two_qubit_gate(num_qubits, getfield(QCO, Symbol(gate_type, "RevGate"))(angle[1]), "q$(qubit_loc[1])", "q$(qubit_loc[2])")
                    end
                elseif length(angle) == 3
                    if (qubit_loc[1] < qubit_loc[2])
                        return QCO.kron_two_qubit_gate(num_qubits, getfield(QCO, Symbol(gate_type, "Gate"))(angle[1], angle[2], angle[3]), "q$(qubit_loc[1])", "q$(qubit_loc[2])")
                    else
                        return QCO.kron_two_qubit_gate(num_qubits, getfield(QCO, Symbol(gate_type, "RevGate"))(angle[1], angle[2], angle[3]), "q$(qubit_loc[1])", "q$(qubit_loc[2])")
                    end
                end

            else 
                Memento.error(_LOGGER, "Enter a valid angle parameter for the input $input gate")
            end

        end

    #-----------------------;
    #   Multi qubit gates   ;
    #-----------------------; 
    elseif gate_type in QCO.MULTI_QUBIT_GATES_ANGLE_PARAMETERS
        if (angle !== nothing) && (length(angle) == 2)
            return getfield(QCO, Symbol(gate_type, "Gate"))(num_qubits, angle[1], angle[2])
        else
            Memento.error(_LOGGER, "Enter a valid angle parameter for the input $input gate")
        end
    end

end

"""
    get_full_sized_kron_gate(input::String, num_qubits::Int64)

Given an input string with kronecker symbols representing the gate and number of qubits of 
the circuit, this function returns a full-sized gate with respect to the input number of qubits. 
For example, if `num_qubits = 3` and the input gate in `I_1xT_2xH_3`, then this function returns 
`IGate⨂TGate⨂HGate`, where IGate, TGate and HGate are single-qubit Identity, T and 
Hadamard gates, respectively. Two qubit gates can also be used as one of the input gates, for ex. `I_1xCV_2_3xH_4`. 
Note that this function currently does not support an input gate parametrized with Euler angles.
"""
function get_full_sized_kron_gate(input::String, num_qubits::Int64)

    kron_gates = QCO._parse_gates_with_kron_symbol(input)
    
    if length(unique(kron_gates)) !== length(kron_gates)
        Memento.error(_LOGGER, "Specify only a single gate per qubit within kron symbol gate $input")
    end

    M = 1
    gate_qubit_locs = []
    for i = 1:length(kron_gates)
        
        gate_type, qubit_loc = QCO._parse_gate_string(kron_gates[i], type = true, qubits = true)

        if !(gate_type in union(QCO.ONE_QUBIT_GATES_CONSTANTS, QCO.TWO_QUBIT_GATES_CONSTANTS))
            Memento.error(_LOGGER, "Specified $input gate is not supported in conjunction with the Kronecker product operation")
        end

        QCO._catch_input_gate_errors(gate_type, qubit_loc, num_qubits, input)
        
        if issubset(qubit_loc, gate_qubit_locs)
            Memento.error(_LOGGER, "Specified qubit(s) for $input gate is incorrect")
        else
            append!(gate_qubit_locs, range(qubit_loc[1], qubit_loc[end]))
        end

        # One-qubit gates
        if (length(qubit_loc) == 1) && (gate_type in QCO.ONE_QUBIT_GATES_CONSTANTS)
            if (gate_type == "I") || (gate_type == "Identity")
                M = kron(M, getfield(QCO, Symbol(gate_type, "Gate"))(1))
            else 
                M = kron(M, getfield(QCO, Symbol(gate_type, "Gate"))())
            end
        
        # Two-qubit gates
        elseif (length(qubit_loc) == 2) && (gate_type in QCO.TWO_QUBIT_GATES_CONSTANTS)
            gate_width = abs(qubit_loc[1] - qubit_loc[2]) 
            
            if gate_width == 1
                if (qubit_loc[1] < qubit_loc[2]) || (gate_type in QCO.TWO_QUBIT_GATES_CONSTANTS_SYMMETRIC)
                    M = kron(M, getfield(QCO, Symbol(gate_type, "Gate"))())
                else
                    M = kron(M, getfield(QCO, Symbol(gate_type, "RevGate"))())
                end
            else
                if (qubit_loc[1] < qubit_loc[2]) || (gate_type in QCO.TWO_QUBIT_GATES_CONSTANTS_SYMMETRIC)
                    kron_gate = QCO.unitary(string(gate_type, qubit_separator, 1, qubit_separator, Int(gate_width + 1)), (gate_width + 1))
                else 
                    kron_gate = QCO.unitary(string(gate_type, qubit_separator, Int(gate_width + 1), qubit_separator, 1), (gate_width + 1))
                end
                M = kron(M, kron_gate)
            end
        end
    end

    QCO._catch_kron_dimension_errors(num_qubits, size(M)[1])
    
    return M
end

"""
    get_input_circuit_dict(input_circuit::Vector{Tuple{Int64,String}}, params::Dict{String,Any})

Given the user input circuit which serves as a warm-start to the optimization model, and user input params dictionary, 
this function outputs the post-processed dictionary of the input circuit which is used by the optimization model. 
"""
function get_input_circuit_dict(
    input_circuit::Vector{Tuple{Int64,String}}, 
    params::Dict{String,Any}
    )
    
    input_circuit_dict = Dict{String, Any}()
    
    # Check if all gates are in elementary gates
    invalid_gate = findfirst(g -> !(g[2] in params["elementary_gates"]), input_circuit)
    if invalid_gate !== nothing
        Memento.warn(_LOGGER, "Neglecting the input circuit as gate $(input_circuit[invalid_gate][2]) is not in input elementary gates")
        return input_circuit_dict
    end
    
    # Process valid gates
    for i = 1:length(input_circuit)
        if i == input_circuit[i][1]
            input_circuit_dict["$i"] = Dict{String, Any}("depth" => input_circuit[i][1], "gate" => input_circuit[i][2])
            # Later: add support for universal and rotation gates here
        else
            Memento.warn(_LOGGER, "Neglecting the input circuit as multiple gates cannot be input at the same depth")
            return Dict{String, Any}()
        end
    end
    
    return input_circuit_dict
end

"""
    _catch_input_gate_errors(gate_type::String, qubit_loc::Vector{Int64}, num_qubits::Int64, input_gate::String)

Given an input gate string, number of qubits of the circuit and the qubit locations for the input gate, 
this function catches and throws any errors, should the input gate type and qubits are invalid. 
"""
function _catch_input_gate_errors(
    gate_type::String, 
    qubit_loc::Vector{Int64}, 
    num_qubits::Int64, 
    input_gate::String; 
    angle = nothing
    )

    if num_qubits <= 0
        Memento.error(_LOGGER, "Specified number of qubits has to be >= 1")
    end
    
    if (gate_type in QCO.TWO_QUBIT_GATES) && (length(qubit_loc) != 2)
        Memento.error(_LOGGER, "Specify two qubits for 2-qubit gate $gate_type in elementary gates")
    elseif (gate_type in QCO.ONE_QUBIT_GATES) && (length(qubit_loc) == 0)
        Memento.error(_LOGGER, "Specify a qubit for the 1-qubit gate $gate_type in elementary gates")
    elseif (gate_type in QCO.ONE_QUBIT_GATES) && (length(qubit_loc) >= 2)
        Memento.error(_LOGGER, "Specify only one qubit for the 1-qubit gate $gate_type in elementary gates")
    end

    if isempty(qubit_loc) && (gate_type !== "GR")
        Memento.error(_LOGGER, "Specify a valid qubit location(s) for the input $input_gate gate")

    elseif (gate_type == "GR") && (!isempty(qubit_loc))
        Memento.error(_LOGGER, "Qubit locations are not necessary for Global-R gate as it is simultaneously applied on all qubits at a depth")
            
    elseif !issubset(qubit_loc, 1:num_qubits)
        Memento.error(_LOGGER, "Specified qubit(s) for $input_gate gate ∉ {1,...,$num_qubits}")

    elseif (length(qubit_loc) == 2) && (isapprox(qubit_loc[1], qubit_loc[2]))
        Memento.error(_LOGGER, "Specified $input_gate gate cannot have identical control and target qubits") 

    elseif length(qubit_loc) > 2 
        Memento.error(_LOGGER, "Only 1- and 2-qubit elementary gates are currently supported")
    end

    if !(gate_type in union(QCO.ONE_QUBIT_GATES_ANGLE_PARAMETERS, QCO.TWO_QUBIT_GATES_ANGLE_PARAMETERS, 
         QCO.MULTI_QUBIT_GATES_ANGLE_PARAMETERS)) && (angle !== nothing)
        Memento.warn(_LOGGER, "Neglecting the angle input for gate $(gate_type) with constant parameters")
    end

end

function _get_R_XYZ_gates_idx(elementary_gates::Array{String,1})
    return findall(x -> (startswith(x, "RX") || startswith(x, "RY") || startswith(x, "RZ")) && !(occursin(kron_symbol, x)), elementary_gates)
end

function _get_R_gates_idx(elementary_gates::Array{String,1})
    return findall(x -> startswith(x, "R") && !(startswith(x, "RX") || startswith(x, "RY") || startswith(x, "RZ")) && !(occursin(kron_symbol, x)), elementary_gates)
end

function _get_GR_gates_idx(elementary_gates::Array{String,1})
    return findall(x -> (startswith(x, "GR")) && !(occursin(kron_symbol, x)), elementary_gates)
end

function _get_U3_gates_idx(elementary_gates::Array{String,1})
    return findall(x -> (startswith(x, "U3")) && !(occursin(kron_symbol, x)), elementary_gates)
end

function _get_CR_gates_idx(elementary_gates::Array{String,1})
    return findall(x -> (startswith(x, "CRX") || startswith(x, "CRY") || startswith(x, "CRZ")) && !(occursin(kron_symbol, x)), elementary_gates)
end

function _get_CU3_gates_idx(elementary_gates::Array{String, 1})
    return findall(x -> (startswith(x, "CU3")) && !(occursin(kron_symbol, x)), elementary_gates)
end

function _get_Phase_gates_idx(elementary_gates::Array{String, 1})
    return findall(x -> (startswith(x, "Phase")) && !(occursin(kron_symbol, x)), elementary_gates)
end

function _get_identity_idx(M::Array{Float64,3})
    [i for i=1:size(M)[3] if isapprox(M[:,:,i], Matrix(I, size(M)[1], size(M)[2]), atol=1E-5)]
end

function _get_gate_idx(gates_dict::Dict{String, Any}, gate_type::String)
    local is_gate_type
    
    if gate_type == "cnot_gate"
        is_gate_type = t -> (gt = QCO._parse_gate_string(t, type=true)) == "CNot" || gt == "CX"
    elseif gate_type == "T_gate"
        is_gate_type = t -> (gt = QCO._parse_gate_string(t, type=true)) == "T" || gt == "Tdagger"
    else
        Memento.error(_LOGGER, "Unsupported gate type: $gate_type")
    end

    [ parse(Int, idx_str) 
      for (idx_str, gate) in pairs(gates_dict)
      if any(
           is_gate_type,
           occursin(kron_symbol, gate["type"][1]) ?
             QCO._parse_gates_with_kron_symbol(gate["type"][1]) :
             (gate["type"][1],)
         )
    ]
end

# Backward compatibility with existing function names
_get_cnot_idx(gates_dict::Dict{String, Any}) = QCO._get_gate_idx(gates_dict, "cnot_gate")
_get_T_idx(gates_dict::Dict{String, Any}) = QCO._get_gate_idx(gates_dict, "T_gate")

function _get_angle_gates_idx(elementary_gates::Array{String,1})

    # Update this list as new gates with angle parameters are added in gates.jl
    R_XYZ_gates_idx = QCO._get_R_XYZ_gates_idx(elementary_gates)
    R_gates_idx     = QCO._get_R_gates_idx(elementary_gates)
    GR_gates_idx    = QCO._get_GR_gates_idx(elementary_gates)
    Phase_gates_idx = QCO._get_Phase_gates_idx(elementary_gates)
    CR_gates_idx    = QCO._get_CR_gates_idx(elementary_gates)
    U3_gates_idx    = QCO._get_U3_gates_idx(elementary_gates)
    CU3_gates_idx   = QCO._get_CU3_gates_idx(elementary_gates)

    one_angle_gates   = union(R_XYZ_gates_idx, Phase_gates_idx, CR_gates_idx)
    two_angle_gates   = union(R_gates_idx, GR_gates_idx)
    three_angle_gates = union(U3_gates_idx, CU3_gates_idx)

    return one_angle_gates, two_angle_gates, three_angle_gates
end

function _get_cnot_bounds!(data::Dict{String, Any}, params::Dict{String, Any})
    cnot_lb = get(params, "set_cnot_lower_bound", 0)
    cnot_ub = get(params, "set_cnot_upper_bound", data["maximum_depth"])
    max_depth = data["maximum_depth"]

    if cnot_lb > max_depth
        Memento.error(_LOGGER, "Invalid lower bound on the number of CNOT gates given the maximum depth")
    elseif cnot_lb < cnot_ub
        cnot_lb > 0 && (data["cnot_lower_bound"] = cnot_lb)
        cnot_ub < max_depth && (data["cnot_upper_bound"] = cnot_ub)
    elseif isapprox(cnot_lb, cnot_ub, atol=1E-6)
        data["cnot_lower_bound"] = cnot_lb
        data["cnot_upper_bound"] = cnot_ub
    else
        Memento.warn(_LOGGER, "Invalid CNot-gate lower/upper bound")
    end

    return data
end
