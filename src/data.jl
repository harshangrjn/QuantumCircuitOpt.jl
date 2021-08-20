import LinearAlgebra: I

function get_data(params::Dict{String, Any}; eliminate_identical_gates = true)
    
    # Number of qubits
    if "num_qubits" in keys(params)
        if params["num_qubits"] < 2 
            Memento.error(_LOGGER, "Minimum of 2 qubits is necessary")
        end
        num_qubits = params["num_qubits"]
    else
        Memento.error(_LOGGER, "Number of qubits has to be specified by the user")
    end

    # Depth
    if "depth" in keys(params)
        if params["depth"] < 2 
            Memento.error(_LOGGER, "Minimum depth of 2 is necessary")
        end
        depth = params["depth"]
    else
        Memento.error(_LOGGER, "Depth of decomposition has to be specified by the user")
    end

    # Elementary gates
    if !("elementary_gates" in keys(params)) || isempty(params["elementary_gates"])
        Memento.error(_LOGGER, "Input elementary gates is empty. Enter at least two unique unitary gates")
    end

    # Initial gate
    if "initial_gate" in keys(params)
        if params["initial_gate"] == "Identity"
            initial_gate = QCO.complex_to_real_matrix(QCO.IGate(num_qubits))
        else 
            Memento.error(_LOGGER, "Currently, only \"Identity\" is supported as an initial gate")
            # Add code here to support non-identity as an initial gate. 
        end
    else
        initial_gate = QCO.complex_to_real_matrix(QCO.IGate(num_qubits))
    end

    # Input Circuit
    if "input_circuit" in keys(params)
        input_circuit = params["input_circuit"]
    else
        # default value
        input_circuit = []
    end
    
    input_circuit_dict = Dict{String,Any}()

    if length(input_circuit) > 0 && (length(input_circuit) <= params["depth"])
        
        input_circuit_dict = QCO.get_input_circuit_dict(input_circuit, params)

    else
        (length(input_circuit) > 0) && (Memento.warn(_LOGGER, "Neglecting the input circuit as it's depth is greater than the allowable depth"))
    end

    # Decomposition type 
    if "decomposition_type" in keys(params)
        decomposition_type = params["decomposition_type"]
    else
        decomposition_type = "exact"
    end

    # Objective function
    if "objective" in keys(params)
        objective = params["objective"]
    else
        objective = "minimize_depth"
    end
    
    # Relax Integrality 
    if "relax_integrality" in keys(params)
        relax_integrality = params["relax_integrality"]
    else
        # default value
        relax_integrality = false
    end

    # Optimizer time limit (in seconds)
    if "time_limit" in keys(params) && params["time_limit"] isa Number
        time_limit = params["time_limit"]
    else
        # default value
        time_limit = 10800
    end

    elementary_gates = unique(params["elementary_gates"])
    
    if length(elementary_gates) < length(params["elementary_gates"])
        Memento.warn(_LOGGER, "Eliminating non-unique gates in the input elementary gates")
    end

    gates_dict, target_real = QCO.get_quantum_gates(params, elementary_gates)

    gates_dict_unique, M_real_unique, identity_idx, cnot_idx = eliminate_nonunique_gates(gates_dict, eliminate_identical_gates = eliminate_identical_gates)
    
    data = Dict{String, Any}("num_qubits" => num_qubits,
                             "depth" => depth,
                             "gates_dict" => gates_dict_unique,
                             "gates_real" => M_real_unique,
                             "initial_gate" => initial_gate,
                             "identity_idx" => identity_idx,
                             "cnot_idx" => cnot_idx,
                             "elementary_gates" => elementary_gates,
                             "target_gate" => target_real,
                             "objective" => objective,
                             "decomposition_type" => decomposition_type,                         
                             "relax_integrality" => relax_integrality,
                             "time_limit" => time_limit
                             )
    
    # Rotation and Universal gate angle discretizations
    data = QCO._populate_data_angle_discretization!(data, params)

    # Input circuit
    if length(keys(input_circuit_dict)) > 0
        data["input_circuit"] = input_circuit_dict
    end

    # Slack Penalty
    if decomposition_type == "approximate"
        if "slack_penalty" in keys(params)
            data["slack_penalty"] = params["slack_penalty"]
        else
            # default value
            data["slack_penalty"] = 1E3
        end
    end
                         
    return data
end

"""
    eliminate_nonunique_gates(gates_dict::Dict{String, Any})

"""
function eliminate_nonunique_gates(gates_dict::Dict{String, Any}; eliminate_identical_gates = false)

    num_gates = length(keys(gates_dict))

    M_real = zeros(2*size(gates_dict["1"]["matrix"])[1], 2*size(gates_dict["1"]["matrix"])[2], num_gates)
    
    for i=1:num_gates
        M_real[:,:,i] = complex_to_real_matrix(gates_dict["$i"]["matrix"])
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

    return gates_dict_unique, M_real_unique, identity_idx, cnot_idx

end

function _get_identity_idx(M::Array{Float64,3})
    
    identity_idx = Int64[]
    
    for i=1:size(M)[3]
        if isapprox(M[:,:,i], Matrix(I, size(M)[1], size(M)[2]), atol=1E-5)   
            push!(identity_idx, i)
        end
    end

    return identity_idx
end

function _get_cnot_idx(gates_dict::Dict{String, Any})
    
    cnot_idx = Int64[]

    # Note: The below objective minimizes both cnot_12 and cnot_21 in the decomposition
    for i in keys(gates_dict)
        if !isempty(findall(x -> startswith(x, "CNot"), gates_dict[i]["type"]))
            push!(cnot_idx, parse(Int64, i))
        end
    end
    
    return cnot_idx
end

function _get_R_gates_idx(elementary_gates::Array{String,1})

    return findall(x -> (startswith(x, "RX") || startswith(x, "RY") || startswith(x, "RZ")) && !(occursin(kron_symbol, x)), elementary_gates)
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

function _get_kron_gates_idx(elementary_gates::Array{String, 1})

    return findall(x -> occursin(kron_symbol, x), elementary_gates)
end

function _populate_data_angle_discretization!(data::Dict{String, Any}, params::Dict{String, Any})

    R_gates_idx   = QCO._get_R_gates_idx(data["elementary_gates"])
    U3_gates_idx  = QCO._get_U3_gates_idx(data["elementary_gates"])
    CR_gates_idx  = QCO._get_CR_gates_idx(data["elementary_gates"])
    CU3_gates_idx = QCO._get_CU3_gates_idx(data["elementary_gates"])
    
    if !isempty(union(R_gates_idx, U3_gates_idx, CR_gates_idx, CU3_gates_idx)) 
        
        data["discretization"] = Dict{String, Any}()
    
        if !isempty(R_gates_idx) 
            for i in R_gates_idx
                if startswith(data["elementary_gates"][i], "RX")
                    data["discretization"]["RX"] = params["RX_discretization"]
                elseif startswith(data["elementary_gates"][i], "RY")
                    data["discretization"]["RY"] = params["RY_discretization"]
                elseif startswith(data["elementary_gates"][i], "RZ")
                    data["discretization"]["RZ"] = params["RZ_discretization"]
                end
            end
        end

        if !isempty(U3_gates_idx)
            for i in U3_gates_idx
                if startswith(data["elementary_gates"][i], "U3")
                    data["discretization"]["U3_θ"] = params["U_θ_discretization"]
                    data["discretization"]["U3_ϕ"] = params["U_ϕ_discretization"]
                    data["discretization"]["U3_λ"] = params["U_λ_discretization"]
                end
            end
        end

        if !isempty(CR_gates_idx) 
            for i in CR_gates_idx
                if startswith(data["elementary_gates"][i], "CRX")
                    data["discretization"]["CRX"] = params["CRX_discretization"]
                elseif startswith(data["elementary_gates"][i], "CRY")
                    data["discretization"]["CRY"] = params["CRY_discretization"]
                elseif startswith(data["elementary_gates"][i], "CRZ")
                    data["discretization"]["CRZ"] = params["CRZ_discretization"]
                end
            end
        end

        if !isempty(CU3_gates_idx)
            for i in CU3_gates_idx
                if startswith(data["elementary_gates"][i], "CU3")
                    data["discretization"]["CU3_θ"] = params["CU_θ_discretization"]
                    data["discretization"]["CU3_ϕ"] = params["CU_ϕ_discretization"]
                    data["discretization"]["CU3_λ"] = params["CU_λ_discretization"]
                end
            end
        end
    end

    return data
end


"""
    get_quantum_gates(params::Dict{String, Any}, elementary_gates::Array{String,1})

Given a vector of input with the names of gates (see examples folder), `get_quantum_gates` function 
returns the corresponding elementary gates in the three-dimensional complex matrix form. 
""" 
function get_quantum_gates(params::Dict{String, Any}, elementary_gates::Array{String,1})

    num_qubits = params["num_qubits"]

    gates_dict = QCO.get_all_gates_dictionary(params, elementary_gates)

    if !("target_gate" in keys(params)) || isempty(params["target_gate"])
        Memento.error(_LOGGER, "Target gate not found in the input data")
    end 

    if (size(params["target_gate"])[1] != size(params["target_gate"])[2]) || (size(params["target_gate"])[1] != 2^params["num_qubits"])
        Memento.error(_LOGGER, "Dimensions of target gate do not match the input num_qubits")
    end
 
    return gates_dict, complex_to_real_matrix(params["target_gate"])
end

function get_all_gates_dictionary(params::Dict{String, Any}, elementary_gates::Array{String,1})

    num_qubits = params["num_qubits"]

    R_gates_idx   = QCO._get_R_gates_idx(elementary_gates)
    U3_gates_idx  = QCO._get_U3_gates_idx(elementary_gates)
    CR_gates_idx  = QCO._get_CR_gates_idx(elementary_gates)
    CU3_gates_idx = QCO._get_CU3_gates_idx(elementary_gates)
    kron_gates_idx = QCO._get_kron_gates_idx(elementary_gates)

    R_complex_dict = Dict{}
    if !isempty(R_gates_idx)
        R_complex_dict = QCO.get_all_R_gates(params, elementary_gates, R_gates_idx)
    end
    
    U3_complex_dict = Dict{}
    if !isempty(U3_gates_idx)
        U3_complex_dict = QCO.get_all_U3_gates(params, elementary_gates, U3_gates_idx)
    end

    CR_complex_dict = Dict{}
    if !isempty(CR_gates_idx)
        CR_complex_dict = QCO.get_all_CR_gates(params, elementary_gates, CR_gates_idx)
    end

    CU3_complex_dict = Dict{}
    if !isempty(CU3_gates_idx)
        CU3_complex_dict = QCO.get_all_CU3_gates(params, elementary_gates, CU3_gates_idx)
    end

    gates_dict = Dict{String, Any}()

    counter = 1

    for i=1:length(elementary_gates)

        if i in union(R_gates_idx, U3_gates_idx, CR_gates_idx, CU3_gates_idx)
            M_elementary_dict = Dict{}

            if i in R_gates_idx
                M_elementary_dict = R_complex_dict
            elseif i in U3_gates_idx
                M_elementary_dict = U3_complex_dict
            elseif i in CR_gates_idx
                M_elementary_dict = CR_complex_dict
            elseif i in CU3_gates_idx
                M_elementary_dict = CU3_complex_dict
            end

            for j in keys(M_elementary_dict) # Gate-type
                if (j == elementary_gates[i])
                    for k in keys(M_elementary_dict[j]) # Angle
                        for l in keys(M_elementary_dict[j][k]["$(num_qubits)qubit_rep"]) # qubits (which will now be 1)
                            
                            M_sqrd = M_elementary_dict[j][k]["$(num_qubits)qubit_rep"][l]^2

                            gates_dict["$counter"] = Dict{String, Any}("type" => [j],
                                                                       "angle" => Any,
                                                                       "qubit_loc" => l,
                                                                       "matrix" => M_elementary_dict[j][k]["$(num_qubits)qubit_rep"][l],
                                                                       "isInvolutory" => isapprox(M_sqrd, Matrix(LA.I, size(M_sqrd)[1], size(M_sqrd)[2]), atol=1E-6))

                            if i in union(R_gates_idx, CR_gates_idx)
                                gates_dict["$counter"]["angle"] = M_elementary_dict[j][k]["angle"]

                            elseif i in union(U3_gates_idx, CU3_gates_idx)
                                gates_dict["$counter"]["angle"] = Dict{String, Any}("θ" => M_elementary_dict[j][k]["θ"],
                                                                                    "ϕ" => M_elementary_dict[j][k]["ϕ"],
                                                                                    "λ" => M_elementary_dict[j][k]["λ"],)
                            end

                            counter += 1

                        end
                    end
                end
            end
        
        # Elementary gates which contain kronecker symbols
        elseif i in kron_gates_idx
            M = QCO.get_full_sized_kron_symbol_gate(elementary_gates[i], num_qubits)
            M_sqrd = M^2

            gates_dict["$counter"] = Dict{String, Any}("type" => [elementary_gates[i]],
                                                       "matrix" => M,
                                                       "isInvolutory" => isapprox(M_sqrd, Matrix(LA.I, size(M_sqrd)[1], size(M_sqrd)[2]), atol=1E-6))
            counter += 1

        else 
            
            M = QCO.get_full_sized_gate(elementary_gates[i], num_qubits)
            M_sqrd = M^2

            gates_dict["$counter"] = Dict{String, Any}("type" => [elementary_gates[i]],
                                                       "matrix" => M,
                                                       "isInvolutory" => isapprox(M_sqrd, Matrix(LA.I, size(M_sqrd)[1], size(M_sqrd)[2]), atol=1E-6))
            counter += 1

        end

    end

    return gates_dict
    
end

function get_all_R_gates(params::Dict{String, Any}, elementary_gates::Array{String,1}, R_gates_idx::Vector{Int64})

    R_complex = Dict{String, Any}()

    if length(R_gates_idx) >= 1 
        for i=1:length(R_gates_idx)

            gate_type = elementary_gates[R_gates_idx[i]]
            
            qubits_string_1 = string.(1:params["num_qubits"])

            if !(gate_type in union("RX_" .* qubits_string_1, "RY_" .* qubits_string_1, "RZ_" .* qubits_string_1))
                Memento.error(_LOGGER, "Input R gate type ($(gate_type)) is not supported.")
            end
            
            # This implies that discretizations are the same for all gates, checks only for RX/RY/RZ discretization
            if isempty(params[string(gate_type[1:2],"_discretization")])
                Memento.error(_LOGGER, "Empty discretization angles for $(gate_type) gate. Input at least one angle")
            end        

            R_complex["$(gate_type)"] = Dict{String, Any}()    
            R_complex["$(gate_type)"] = QCO.get_discretized_R_gates(gate_type, R_complex[gate_type], collect(params[string(gate_type[1:2],"_discretization")]), params["num_qubits"])

        end
    end
    
    return R_complex    
end

function get_discretized_R_gates(gate_type::String, R::Dict{String, Any}, discretization::Array{Float64,1}, num_qubits::Int64)

    if length(discretization) >= 1

        for i=1:length(discretization)
            R["angle_$i"] = Dict{String, Any}("angle" => discretization[i],
                                             "$(num_qubits)qubit_rep" => Dict{String, Any}() )
            
            R["angle_$i"]["$(num_qubits)qubit_rep"]["qubit_$(gate_type[4])"] = QCO.get_full_sized_gate(gate_type, num_qubits, angle = discretization[i])
        end
    end 

    return R 
end

function get_all_CR_gates(params::Dict{String, Any}, elementary_gates::Array{String,1}, CR_gates_idx::Vector{Int64})

    CR_complex = Dict{String, Any}()

    if length(CR_gates_idx) >= 1 
        for i=1:length(CR_gates_idx)

            gate_type = elementary_gates[CR_gates_idx[i]]
            # CRP_12
            if !(startswith(gate_type, "CRX") || startswith(gate_type, "CRY") || startswith(gate_type, "CRZ"))
                Memento.error(_LOGGER, "Input controlled-rotation (CR) gate type ($(gate_type)) is not supported.")
            end
            
            # This implies that discretizations are the same for all gates, checks only for RX/RY/RZ discretization
            if isempty(params[string(gate_type[1:3],"_discretization")])
                Memento.error(_LOGGER, "Empty discretization angles for $(gate_type) gate. Input a rotation angle")
            end        

            CR_complex["$(gate_type)"] = Dict{String, Any}()    
            CR_complex["$(gate_type)"] = QCO.get_discretized_CR_gates(gate_type, CR_complex[gate_type], collect(params[string(gate_type[1:3],"_discretization")]), params["num_qubits"])

        end
    end

    return CR_complex
end

function get_discretized_CR_gates(gate_type::String, CR::Dict{String, Any}, discretization::Array{Float64,1}, num_qubits::Int64)

    if length(discretization) >= 1

        for i=1:length(discretization)
            angle = discretization[i]
            CR["angle_$i"] = Dict{String, Any}("angle" => discretization[i],
                                               "$(num_qubits)qubit_rep" => Dict{String, Any}()
                                              )

            CR["angle_$i"]["$(num_qubits)qubit_rep"]["qubit_$(gate_type[5:6])"] = QCO.get_full_sized_gate(gate_type, num_qubits, angle = angle)
        end
    end 

    return CR 
end

function get_all_U3_gates(params::Dict{String, Any}, elementary_gates::Array{String,1}, U3_gates_idx::Vector{Int64})

    U3_complex = Dict{String, Any}()

    if length(U3_gates_idx) >= 1 
        
        for i=1:length(U3_gates_idx)
            gate_name = elementary_gates[U3_gates_idx[i]]
            if startswith(gate_name, "U3")        
                U3_complex[gate_name] = Dict{String, Any}()    
                
                for angle in ["θ", "ϕ", "λ"]
                    if isempty(params["U_$(angle)_discretization"])
                        Memento.error(_LOGGER, "Empty $(angle) discretization angle for U3 gate. Input at least one angle")
                    end
                end

                U3_complex[gate_name] = QCO.get_discretized_U3_gates(gate_name, U3_complex[gate_name], collect(float(params["U_θ_discretization"])), collect(float(params["U_ϕ_discretization"])), collect(float(params["U_λ_discretization"])), params["num_qubits"])
            end

        end
    end
    
    return U3_complex    
end

function get_discretized_U3_gates(gate_type::String, U3::Dict{String, Any}, θ_discretization::Array{Float64,1}, ϕ_discretization::Array{Float64,1}, λ_discretization::Array{Float64,1}, num_qubits::Int64) 

    counter = 1
    for i=1:length(θ_discretization)
        for j=1:length(ϕ_discretization)
            for k=1:length(λ_discretization)
                
                U3["angle_$(counter)"] = Dict{String, Any}("θ" => θ_discretization[i],
                                                           "ϕ" => ϕ_discretization[j],
                                                           "λ" => λ_discretization[k],
                                                           "$(num_qubits)qubit_rep" => Dict{String, Any}()
                                                          )
                
                U3["angle_$(counter)"]["$(num_qubits)qubit_rep"]["qubit_$(gate_type[4])"] = QCO.get_full_sized_gate(gate_type, num_qubits, angle = [θ_discretization[i], ϕ_discretization[j], λ_discretization[k]])

                counter += 1
            end
        end
    end
    
    return U3
end

function get_all_CU3_gates(params::Dict{String, Any}, elementary_gates::Array{String,1}, CU3_gates_idx::Vector{Int64})

    CU3_complex = Dict{String, Any}()

    if length(CU3_gates_idx) >= 1 
        
        for i=1:length(CU3_gates_idx)
            gate_name = elementary_gates[CU3_gates_idx[i]]
            if startswith(gate_name, "CU3")        
                CU3_complex[gate_name] = Dict{String, Any}()    
                
                for angle in ["θ", "ϕ", "λ"]
                    if isempty(params["CU_$(angle)_discretization"])
                        Memento.error(_LOGGER, "Empty $(angle) discretization angle for CU3 gate. Input at least one angle")
                    end
                end

                CU3_complex[gate_name] = QCO.get_discretized_CU3_gates(gate_name, CU3_complex[gate_name], collect(float(params["CU_θ_discretization"])), collect(float(params["CU_ϕ_discretization"])), collect(float(params["CU_λ_discretization"])), params["num_qubits"])
            end

        end
    end
    
    return CU3_complex  
end

function get_discretized_CU3_gates(gate_type::String, CU3::Dict{String, Any}, θ_discretization::Array{Float64,1}, ϕ_discretization::Array{Float64,1}, λ_discretization::Array{Float64,1}, num_qubits::Int64)
    
    counter = 1 

    for i=1:length(θ_discretization)
        for j=1:length(ϕ_discretization)
            for k=1:length(λ_discretization)

                angles = [θ_discretization[i], ϕ_discretization[j], λ_discretization[k]]

                CU3["angle_$(counter)"] = Dict{String, Any}("θ" => θ_discretization[i],
                                                            "ϕ" => ϕ_discretization[j],
                                                            "λ" => λ_discretization[k],
                                                            "$(num_qubits)qubit_rep" => Dict{String, Any}()
                                                           )
                
                CU3["angle_$(counter)"]["$(num_qubits)qubit_rep"]["qubit_$(gate_type[5:6])"] = QCO.get_full_sized_gate(gate_type, num_qubits, angle = angles)

                counter += 1
            end
        end
    end
    
    return CU3
end

"""
    get_full_sized_gate(input::String, num_qubits::Int64; target_gate = nothing)

For a given string and number of qubits in the input specified input, this function returns a full 
sized gate with respect to the input number of qubits. 
"""
function get_full_sized_gate(input::String, num_qubits::Int64; angle = nothing)

    if input == "Identity"
        return QCO.IGate(num_qubits)
    end

    if num_qubits > 5
        Memento.error(_LOGGER, "Gates with greater than 5 qubits are not currently supported")
    end

    qubits_string_1, qubits_string_2 = QCO._get_qubit_strings(num_qubits)
    
    #----------------------;
    #   One qubit gates    ;
    #----------------------; 
    if input in "H_" .* qubits_string_1
        return QCO.kron_single_qubit_gate(num_qubits, QCO.HGate(), "q$(input[end])")

    elseif input in "T_" .* qubits_string_1
        return QCO.kron_single_qubit_gate(num_qubits, QCO.TGate(), "q$(input[end])")

    elseif input in "Tdagger_" .* qubits_string_1
        return QCO.kron_single_qubit_gate(num_qubits, QCO.TdaggerGate(), "q$(input[end])")

    elseif input in "S_" .* qubits_string_1
        return QCO.kron_single_qubit_gate(num_qubits, QCO.SGate(), "q$(input[end])")
    
    elseif input in "Sdagger_" .* qubits_string_1
        return QCO.kron_single_qubit_gate(num_qubits, QCO.SdaggerGate(), "q$(input[end])")

    elseif input in "SX_" .* qubits_string_1
        return QCO.kron_single_qubit_gate(num_qubits, QCO.SXGate(), "q$(input[end])")

    elseif input in "SXdagger_" .* qubits_string_1
        return QCO.kron_single_qubit_gate(num_qubits, QCO.SXdaggerGate(), "q$(input[end])")

    elseif input in "X_" .* qubits_string_1
        return QCO.kron_single_qubit_gate(num_qubits, QCO.XGate(), "q$(input[end])")

    elseif input in "Y_" .* qubits_string_1
        return QCO.kron_single_qubit_gate(num_qubits, QCO.YGate(), "q$(input[end])")

    elseif input in "Z_" .* qubits_string_1
        return QCO.kron_single_qubit_gate(num_qubits, QCO.ZGate(), "q$(input[end])") 

    # Gates with continuous angle parameters
    elseif input in "RX_" .* qubits_string_1
        return QCO.kron_single_qubit_gate(num_qubits, QCO.RXGate(angle), "q$(input[end])")

    elseif input in "RY_" .* qubits_string_1
        return QCO.kron_single_qubit_gate(num_qubits, QCO.RYGate(angle), "q$(input[end])")

    elseif input in "RZ_" .* qubits_string_1
        return QCO.kron_single_qubit_gate(num_qubits, QCO.RZGate(angle), "q$(input[end])")

    elseif input in "U3_" .* qubits_string_1
        return QCO.kron_single_qubit_gate(num_qubits, QCO.U3Gate(angle[1], angle[2], angle[3]), "q$(input[end])")

    #----------------------;
    #   Two qubit gates    ;
    #----------------------;
    elseif input in "CNot_" .* qubits_string_2
        c_qubit = parse(Int, input[end-1])
        t_qubit = parse(Int, input[end])

        if c_qubit < t_qubit
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CNotGate(), "q$c_qubit", "q$t_qubit")
        else
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CNotRevGate(), "q$c_qubit", "q$t_qubit")
        end

    elseif input in "CRX_" .* qubits_string_2
        c_qubit = parse(Int, input[end-1])
        t_qubit = parse(Int, input[end])

        if c_qubit < t_qubit
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CRXGate(angle), "q$c_qubit", "q$t_qubit")
        else
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CRXRevGate(angle), "q$c_qubit", "q$t_qubit")
        end

    elseif input in "CRY_" .* qubits_string_2
        c_qubit = parse(Int, input[end-1])
        t_qubit = parse(Int, input[end])

        if c_qubit < t_qubit
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CRYGate(angle), "q$c_qubit", "q$t_qubit")
        else
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CRYRevGate(angle), "q$c_qubit", "q$t_qubit")
        end

    elseif input in "CRZ_" .* qubits_string_2
        c_qubit = parse(Int, input[end-1])
        t_qubit = parse(Int, input[end])

        if c_qubit < t_qubit
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CRZGate(angle), "q$c_qubit", "q$t_qubit")
        else
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CRZRevGate(angle), "q$c_qubit", "q$t_qubit")
        end

    elseif input in "CU3_" .* qubits_string_2
        c_qubit = parse(Int, input[end-1])
        t_qubit = parse(Int, input[end])

        if c_qubit < t_qubit
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CU3Gate(angle[1], angle[2], angle[3]), "q$c_qubit", "q$t_qubit")
        else
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CU3RevGate(angle[1], angle[2], angle[3]), "q$c_qubit", "q$t_qubit")
        end

    elseif input in "CV_" .* qubits_string_2
        c_qubit = parse(Int, input[end-1])
        t_qubit = parse(Int, input[end])

        if c_qubit < t_qubit
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CVGate(), "q$c_qubit", "q$t_qubit")
        else
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CVRevGate(), "q$c_qubit", "q$t_qubit")
        end

    elseif input in "CVdagger_" .* qubits_string_2
        c_qubit = parse(Int, input[end-1])
        t_qubit = parse(Int, input[end])

        if c_qubit < t_qubit
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CVdaggerGate(), "q$c_qubit", "q$t_qubit")
        else
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CVRevdaggerGate(), "q$c_qubit", "q$t_qubit")
        end

    elseif input in "CH_" .* qubits_string_2
        c_qubit = parse(Int, input[end-1])
        t_qubit = parse(Int, input[end])

        if c_qubit < t_qubit
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CHGate(), "q$c_qubit", "q$t_qubit")
        else
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CHRevGate(), "q$c_qubit", "q$t_qubit")
        end

    elseif input in "CX_" .* qubits_string_2
        c_qubit = parse(Int, input[end-1])
        t_qubit = parse(Int, input[end])

        if c_qubit < t_qubit
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CXGate(), "q$c_qubit", "q$t_qubit")
        else
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CXRevGate(), "q$c_qubit", "q$t_qubit")
        end

    elseif input in "CY_" .* qubits_string_2
        c_qubit = parse(Int, input[end-1])
        t_qubit = parse(Int, input[end])

        if c_qubit < t_qubit
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CYGate(), "q$c_qubit", "q$t_qubit")
        else
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CYRevGate(), "q$c_qubit", "q$t_qubit")
        end

    elseif input in "CZ_" .* qubits_string_2
        c_qubit = parse(Int, input[end-1])
        t_qubit = parse(Int, input[end])

        if c_qubit < t_qubit
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CZGate(), "q$c_qubit", "q$t_qubit")
        else
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CZRevGate(), "q$c_qubit", "q$t_qubit")
        end

    elseif input in "CSX_" .* qubits_string_2
        c_qubit = parse(Int, input[end-1])
        t_qubit = parse(Int, input[end])

        if c_qubit < t_qubit
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CSXGate(), "q$c_qubit", "q$t_qubit")
        else
            return QCO.kron_double_qubit_gate(num_qubits, QCO.CSXRevGate(), "q$c_qubit", "q$t_qubit")
        end

    elseif input in "Swap_" .* qubits_string_2
        c_qubit = parse(Int, input[end-1])
        t_qubit = parse(Int, input[end])

        return QCO.kron_double_qubit_gate(num_qubits, QCO.SwapGate(), "q$c_qubit", "q$t_qubit")

    elseif input in "iSwap_" .* qubits_string_2
        c_qubit = parse(Int, input[end-1])
        t_qubit = parse(Int, input[end])

        return QCO.kron_double_qubit_gate(num_qubits, QCO.iSwapGate(), "q$c_qubit", "q$t_qubit")

    else 
        Memento.error(_LOGGER, "Specified input elementary gates or the target gate does not exist in the predefined set of gates.")
    end

end

function get_full_sized_kron_symbol_gate(input::String, num_qubits::Int64)

    qubits_string_1, qubits_string_2 = QCO._get_qubit_strings(num_qubits)

    kron_gates = QCO._parse_gates_with_kron_symbol(input)
    
    M = 1

    for i = 1:length(kron_gates)
        
        if kron_gates[i] in "I_" .* qubits_string_1
            M = kron(M, QCO.IGate(1))
            
        elseif kron_gates[i] in "H_" .* qubits_string_1
            M = kron(M, QCO.HGate())

        elseif kron_gates[i] in "T_" .* qubits_string_1
            M = kron(M, QCO.TGate())
    
        elseif kron_gates[i] in "Tdagger_" .* qubits_string_1
            M = kron(M, QCO.TdaggerGate())
    
        elseif kron_gates[i] in "S_" .* qubits_string_1
            M = kron(M, QCO.SGate())
        
        elseif kron_gates[i] in "Sdagger_" .* qubits_string_1
            M = kron(M, QCO.SdaggerGate())
    
        elseif kron_gates[i] in "SX_" .* qubits_string_1
            M = kron(M, QCO.SXGate())
    
        elseif kron_gates[i] in "SXdagger_" .* qubits_string_1
            M = kron(M, QCO.SXdaggerGate())
    
        elseif kron_gates[i] in "X_" .* qubits_string_1
            M = kron(M, QCO.XGate())
    
        elseif kron_gates[i] in "Y_" .* qubits_string_1
            M = kron(M, QCO.YGate())
    
        elseif kron_gates[i] in "Z_" .* qubits_string_1
            M = kron(M, QCO.ZGate())

        elseif kron_gates[i] in "CNot_" .* qubits_string_2
            c_qubit = parse(Int, kron_gates[i][end-1])
            t_qubit = parse(Int, kron_gates[i][end])
    
            if c_qubit < t_qubit
                M = kron(M, QCO.CNotGate())
            else
                M = kron(M, QCO.CNotRevGate())
            end
    
        elseif kron_gates[i] in "CV_" .* qubits_string_2
            c_qubit = parse(Int, kron_gates[i][end-1])
            t_qubit = parse(Int, kron_gates[i][end])
    
            if c_qubit < t_qubit
                M = kron(M, QCO.CVGate())
            else
                M = kron(M, QCO.CVRevGate())
            end

        elseif kron_gates[i] in "CVdagger_" .* qubits_string_2
            c_qubit = parse(Int, kron_gates[i][end-1])
            t_qubit = parse(Int, kron_gates[i][end])
    
            if c_qubit < t_qubit
                M = kron(M, QCO.CVdaggerGate())
            else
                M = kron(M, QCO.CVdaggerRevGate())
            end

        elseif kron_gates[i] in "CH_" .* qubits_string_2
            c_qubit = parse(Int, kron_gates[i][end-1])
            t_qubit = parse(Int, kron_gates[i][end])
    
            if c_qubit < t_qubit
                M = kron(M, QCO.CHGate())
            else
                M = kron(M, QCO.CHRevGate())
            end

        elseif kron_gates[i] in "CX_" .* qubits_string_2
            c_qubit = parse(Int, kron_gates[i][end-1])
            t_qubit = parse(Int, kron_gates[i][end])
    
            if c_qubit < t_qubit
                M = kron(M, QCO.CXGate())
            else
                M = kron(M, QCO.CXRevGate())
            end

        elseif kron_gates[i] in "CY_" .* qubits_string_2
            c_qubit = parse(Int, kron_gates[i][end-1])
            t_qubit = parse(Int, kron_gates[i][end])
    
            if c_qubit < t_qubit
                M = kron(M, QCO.CYGate())
            else
                M = kron(M, QCO.CYRevGate())
            end

        elseif kron_gates[i] in "CZ_" .* qubits_string_2
            c_qubit = parse(Int, kron_gates[i][end-1])
            t_qubit = parse(Int, kron_gates[i][end])
    
            if c_qubit < t_qubit
                M = kron(M, QCO.CZGate())
            else
                M = kron(M, QCO.CZRevGate())
            end

        elseif kron_gates[i] in "CSX_" .* qubits_string_2
            c_qubit = parse(Int, kron_gates[i][end-1])
            t_qubit = parse(Int, kron_gates[i][end])
    
            if c_qubit < t_qubit
                M = kron(M, QCO.CSXGate())
            else
                M = kron(M, QCO.CSXRevGate())
            end

        elseif kron_gates[i] in "Swap_" .* qubits_string_2
                M = kron(M, QCO.SwapGate())

        elseif kron_gates[i] in "iSwap_" .* qubits_string_2
                M = kron(M, QCO.iSwapGate())

        else 
            # Gate with angle parameters is not supported yet
            Memento.error(_LOGGER, "Specified $(kron_gates[i]) gate is not supported in conjunction with the Kronecker symbol in $(input) elementary gate")
        end

    end

    if size(M)[1] == 2^(num_qubits)
        return M 
    else 
        Memento.error(_LOGGER, "Dimensions mismatch in evaluation of elementary gates with Kronecker symbol")
    end

end

function _get_qubit_strings(num_qubits::Int)
    qubits_string_1 = string.(1:num_qubits)
    qubits_string_2 = []
    
    for i=1:num_qubits
        for j=1:num_qubits
            (i != j) && (push!(qubits_string_2, qubits_string_1[i] .* qubits_string_1[j]))
        end
    end

    return qubits_string_1, qubits_string_2    
end

function get_input_circuit_dict(input_circuit::Vector{Tuple{Int64,String}}, params::Dict{String,Any})

    input_circuit_dict = Dict{String, Any}()

    status = true
    gate_type = []
    
    for i = 1:length(input_circuit)
        if !(input_circuit[i][2] in params["elementary_gates"])
            status = false
            gate_type = input_circuit[i][2]
        end
    end

    if status  
        for i = 1:length(input_circuit)
    
            if i == input_circuit[i][1]
                input_circuit_dict["$i"] = Dict{String, Any}("depth" => input_circuit[i][1],
                                                             "gate" => input_circuit[i][2])
                # Later: add support for universal and rotation gates here
            else
                input_circuit_dict = Dict{String, Any}()          
                Memento.warn(_LOGGER, "Neglecting the input circuit as multiple gates cannot be input at the same depth")
                break
            end

        end
    else
        Memento.warn(_LOGGER, "Neglecting the input circuit as gate $gate_type is not in input elementary gates")
    end
    
    return input_circuit_dict
end 
