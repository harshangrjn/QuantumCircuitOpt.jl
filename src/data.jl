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

    R_gates_ids = findall(x -> startswith(x, "R"), data["elementary_gates"])
    U3_gates_ids = findall(x -> startswith(x, "U3"), data["elementary_gates"])
    CR_gates_ids = findall(x -> startswith(x, "CR"), data["elementary_gates"])
    CU3_gates_ids = findall(x -> startswith(x, "CU3"), data["elementary_gates"])
    
    if !isempty(R_gates_ids) || !isempty(U3_gates_ids) || !isempty(CR_gates_ids) || !isempty(CU3_gates_ids)
        data["discretization"] = Dict{String, Any}()
    end

    if !isempty(R_gates_ids) 
        for i in R_gates_ids
            if startswith(data["elementary_gates"][i], "RX")
                data["discretization"]["RX"] = params["RX_discretization"]
            elseif startswith(data["elementary_gates"][i], "RY")
                data["discretization"]["RY"] = params["RY_discretization"]
            elseif startswith(data["elementary_gates"][i], "RZ")
                data["discretization"]["RZ"] = params["RZ_discretization"]
            end
        end
    end

    if !isempty(U3_gates_ids)
        for i in U3_gates_ids
            if startswith(data["elementary_gates"][i], "U3")
                data["discretization"]["U3_θ"] = params["U_θ_discretization"]
                data["discretization"]["U3_ϕ"] = params["U_ϕ_discretization"]
                data["discretization"]["U3_λ"] = params["U_λ_discretization"]
            end
        end
    end

    if !isempty(CR_gates_ids) 
        for i in CR_gates_ids
            if startswith(data["elementary_gates"][i], "CRX")
                data["discretization"]["CRX"] = params["CRX_discretization"]
            elseif startswith(data["elementary_gates"][i], "CRY")
                data["discretization"]["CRY"] = params["CRY_discretization"]
            elseif startswith(data["elementary_gates"][i], "CRZ")
                data["discretization"]["CRZ"] = params["CRZ_discretization"]
            end
        end
    end

    if !isempty(CU3_gates_ids)
        for i in CU3_gates_ids
            if startswith(data["elementary_gates"][i], "CU3")
                data["discretization"]["CU3_θ"] = params["CU_θ_discretization"]
                data["discretization"]["CU3_ϕ"] = params["CU_ϕ_discretization"]
                data["discretization"]["CU3_λ"] = params["CU_λ_discretization"]
            end
        end
    end

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

    R_gates_ids = findall(x -> startswith(x, "RX") || startswith(x, "RY") || startswith(x, "RZ"), elementary_gates)
    U3_gates_ids = findall(x -> startswith(x, "U3"), elementary_gates)
    CR_gates_ids = findall(x -> startswith(x, "CRX") || startswith(x, "CRY") || startswith(x, "CRZ"), elementary_gates)
    CU3_gates_ids = findall(x -> startswith(x, "CU3"), elementary_gates)

    R_complex_dict = Dict{}
    if !isempty(R_gates_ids)
        R_complex_dict = QCO.get_all_R_gates(params, elementary_gates)
    end
    
    U3_complex_dict = Dict{}
    if !isempty(U3_gates_ids)
        U3_complex_dict = QCO.get_all_U3_gates(params, elementary_gates)
    end

    CR_complex_dict = Dict{}
    if !isempty(CR_gates_ids)
        CR_complex_dict = QCO.get_all_CR_gates(params, elementary_gates)
    end

    CU3_complex_dict = Dict{}
    if !isempty(CU3_gates_ids)
        CU3_complex_dict = QCO.get_all_CU3_gates(params, elementary_gates)
    end

    gates_dict = Dict{String, Any}()

    counter = 1

    for i=1:length(elementary_gates)

        if startswith(elementary_gates[i], "R") || startswith(elementary_gates[i], "U3") || startswith(elementary_gates[i], "CR") || startswith(elementary_gates[i], "CU3")
            M_elementary_dict = Dict{}

            if startswith(elementary_gates[i], "R")
                M_elementary_dict = R_complex_dict
            elseif startswith(elementary_gates[i], "U3")
                M_elementary_dict = U3_complex_dict
            elseif startswith(elementary_gates[i], "CR")
                M_elementary_dict = CR_complex_dict
            elseif startswith(elementary_gates[i], "CU3")
                M_elementary_dict = CU3_complex_dict
            end

            for j in keys(M_elementary_dict) # Gate type
                if (j == elementary_gates[i])
                    for k in keys(M_elementary_dict[j]) # Angle
                        for l in keys(M_elementary_dict[j][k]["$(num_qubits)qubit_rep"]) # qubits (which will now be 1)
                            
                            M_sqrd = M_elementary_dict[j][k]["$(num_qubits)qubit_rep"][l]^2

                            gates_dict["$counter"] = Dict{String, Any}("type" => [j],
                                                                       "angle" => Any,
                                                                       "qubit_loc" => l,
                                                                       "matrix" => M_elementary_dict[j][k]["$(num_qubits)qubit_rep"][l],
                                                                       "isInvolutory" => isapprox(M_sqrd, Matrix(LA.I, size(M_sqrd)[1], size(M_sqrd)[2]), atol=1E-6))

                            if startswith(elementary_gates[i], "R") || startswith(elementary_gates[i], "CR")
                                gates_dict["$counter"]["angle"] = M_elementary_dict[j][k]["angle"]

                            elseif startswith(elementary_gates[i], "U3") || startswith(elementary_gates[i], "CU3")
                                gates_dict["$counter"]["angle"] = Dict{String, Any}("θ" => M_elementary_dict[j][k]["θ"],
                                                                                    "ϕ" => M_elementary_dict[j][k]["ϕ"],
                                                                                    "λ" => M_elementary_dict[j][k]["λ"],)
                            end

                            counter += 1

                        end
                    end
                end
            end
            
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

function get_all_CR_gates(params::Dict{String, Any}, elementary_gates::Array{String,1})
    
    CR_gates_ids = findall(x -> startswith(x, "CRX") || startswith(x, "CRY") || startswith(x, "CRZ"), elementary_gates)

    CR_complex = Dict{String, Any}()

    if length(CR_gates_ids) >= 1 
        for i=1:length(CR_gates_ids)

            gate_type = elementary_gates[CR_gates_ids[i]]
            # CRP_12
            if !(startswith(gate_type, "CRX") || startswith(gate_type, "CRY") || startswith(gate_type, "CRZ"))
                Memento.error(_LOGGER, "Input controlled-rotation (CR) gate type ($(gate_type)) is not supported.")
            end
            
            # This implies that discretizations are the same for all gates, checks only for RX/RY/RZ discretization
            if isempty(params[string(gate_type[1:3],"_discretization")])
                Memento.error(_LOGGER, "Empty discretization angles for $(gate_type) gate. Input at least one angle")
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
            if startswith(gate_type, "CRX")
                if gate_type[5] < gate_type[6]
                    CR_discrete = QCO.CRXGate(discretization[i])
                else
                    CR_discrete = QCO.CRXRevGate(discretization[i])
                end
            elseif startswith(gate_type, "CRY")
                if gate_type[5] < gate_type[6]
                    CR_discrete = QCO.CRYGate(discretization[i])
                else
                    CR_discrete = QCO.CRYRevGate(discretization[i])
                end
            elseif startswith(gate_type, "CRZ")
                if gate_type[5] < gate_type[6]
                    CR_discrete = QCO.CRZGate(discretization[i])
                else
                    CR_discrete = QCO.CRZRevGate(discretization[i])
                end
            end

            angle = discretization[i]
            CR["angle_$i"] = Dict{String, Any}("angle" => discretization[i],
                                             "$(num_qubits)qubit_rep" => Dict{String, Any}()
                                            )

            CR["angle_$i"]["$(num_qubits)qubit_rep"]["qubit_$(gate_type[5:6])"] = QCO.get_full_sized_gate(gate_type, num_qubits, matrix = CR_discrete, qubit_location = "q$(gate_type[5:6])", target_angle = angle)

        end
    end 

    return CR 
end

function get_all_R_gates(params::Dict{String, Any}, elementary_gates::Array{String,1})

    R_gates_ids = findall(x -> startswith(x, "RX") || startswith(x, "RY") || startswith(x, "RZ"), elementary_gates)

    R_complex = Dict{String, Any}()

    if length(R_gates_ids) >= 1 
        for i=1:length(R_gates_ids)

            gate_type = elementary_gates[R_gates_ids[i]]

            if !(gate_type in ["RX_1", "RX_2", "RX_3", "RY_1", "RY_2", "RY_3", "RZ_1", "RZ_2", "RZ_3"])
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
            if startswith(gate_type, "RX")
                R_discrete = QCO.RXGate(discretization[i])
            elseif startswith(gate_type, "RY")
                R_discrete = QCO.RYGate(discretization[i])
            elseif startswith(gate_type, "RZ")
                R_discrete = QCO.RZGate(discretization[i])
            end

            R["angle_$i"] = Dict{String, Any}("angle" => discretization[i],
                                             "1qubit_rep" => R_discrete,
                                             "$(num_qubits)qubit_rep" => Dict{String, Any}()
                                            )
            
            R["angle_$i"]["$(num_qubits)qubit_rep"]["qubit_$(gate_type[4])"] = QCO.get_full_sized_gate(gate_type[1:2], num_qubits, matrix = R_discrete, qubit_location = "q$(gate_type[4])")

        end
    end 

    return R 
end

function get_all_CU3_gates(params::Dict{String, Any}, elementary_gates::Array{String,1})

    CU3_gates_ids = findall(x -> startswith(x, "CU3"), elementary_gates)

    CU3_complex = Dict{String, Any}()

    if length(CU3_gates_ids) >= 1 
        
        for i=1:length(CU3_gates_ids)
            gate_name = elementary_gates[CU3_gates_ids[i]]
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

                if gate_type[5] < gate_type[6]
                    CU3_discrete = QCO.CU3Gate(θ_discretization[i], ϕ_discretization[j], λ_discretization[k])
                else
                    CU3_discrete = QCO.CU3RevGate(θ_discretization[i], ϕ_discretization[j], λ_discretization[k])
                end

                angles = [θ_discretization[i], ϕ_discretization[j], λ_discretization[k]]

                CU3["angle_$(counter)"] = Dict{String, Any}("θ" => θ_discretization[i],
                                                            "ϕ" => ϕ_discretization[j],
                                                            "λ" => λ_discretization[k],
                                                            "$(num_qubits)qubit_rep" => Dict{String, Any}()
                                                           )
                
                CU3["angle_$(counter)"]["$(num_qubits)qubit_rep"]["qubit_$(gate_type[5:6])"] = QCO.get_full_sized_gate(gate_type, num_qubits, matrix = CU3_discrete, qubit_location = "q$(gate_type[5:6])", target_angle = angles)

                counter += 1
            end
        end
    end
    
    return CU3
end

function get_all_U3_gates(params::Dict{String, Any}, elementary_gates::Array{String,1})

    U3_gates_ids = findall(x -> startswith(x, "U3"), elementary_gates)

    U3_complex = Dict{String, Any}()

    if length(U3_gates_ids) >= 1 
        
        for i=1:length(U3_gates_ids)
            gate_name = elementary_gates[U3_gates_ids[i]]
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
                
                U3_discrete = QCO.U3Gate(θ_discretization[i], ϕ_discretization[j], λ_discretization[k])
                
                U3["angle_$(counter)"] = Dict{String, Any}("θ" => θ_discretization[i],
                                                           "ϕ" => ϕ_discretization[j],
                                                           "λ" => λ_discretization[k],
                                                           "1qubit_rep" => U3_discrete,
                                                           "$(num_qubits)qubit_rep" => Dict{String, Any}()
                                                          )
                
                U3["angle_$(counter)"]["$(num_qubits)qubit_rep"]["qubit_$(gate_type[4])"] = QCO.get_full_sized_gate(gate_type[1:2], num_qubits, matrix = U3_discrete, qubit_location = "q$(gate_type[4])")

                counter += 1
            end
        end
    end
    
    return U3
end

"""
    get_full_sized_gate(input::String, num_qubits::Int64; matrix = nothing, qubit_location = nothing, target_gate = nothing)

For a given string and number of qubits in the input specified input, this function returns a full 
sized gate with respect to the input number of qubits. 
"""
function get_full_sized_gate(input::String, num_qubits::Int64; matrix = nothing, qubit_location = nothing, target_angle = nothing)

    if input == "Identity"
        return QCO.IGate(num_qubits)
    end

    if num_qubits > 3
        Memento.error(_LOGGER, "Gates with greater than 3 qubits are not currently supported")
    end
    
    if input == "H_1"
        return QCO.kron_single_gate(num_qubits, QCO.HGate(), "q1")

    elseif input == "H_2"
        return QCO.kron_single_gate(num_qubits, QCO.HGate(), "q2")

    elseif input == "T_1"
        return QCO.kron_single_gate(num_qubits, QCO.TGate(), "q1")

    elseif input == "T_2"
        return QCO.kron_single_gate(num_qubits, QCO.TGate(), "q2")

    elseif input == "Tdagger_1"
        return QCO.kron_single_gate(num_qubits, QCO.TdaggerGate(), "q1")

    elseif input == "Tdagger_2"
        return QCO.kron_single_gate(num_qubits, QCO.TdaggerGate(), "q2")   

    elseif input == "S_1"
        return QCO.kron_single_gate(num_qubits, QCO.SGate(), "q1") 

    elseif input == "S_2"
        return QCO.kron_single_gate(num_qubits, QCO.SGate(), "q2")  

    elseif input == "Sdagger_1"
        return QCO.kron_single_gate(num_qubits, QCO.SdaggerGate(), "q1") 

    elseif input == "Sdagger_2"
        return QCO.kron_single_gate(num_qubits, QCO.SdaggerGate(), "q2")   

    elseif input == "SX_1"
        return QCO.kron_single_gate(num_qubits, QCO.SXGate(), "q1") 

    elseif input == "SX_2"
        return QCO.kron_single_gate(num_qubits, QCO.SXGate(), "q2")  

    elseif input == "SXdagger_1"
        return QCO.kron_single_gate(num_qubits, QCO.SXdaggerGate(), "q1")

    elseif input == "SXdagger_2"
        return QCO.kron_single_gate(num_qubits, QCO.SXdaggerGate(), "q2")   

    elseif input == "X_1"
        return QCO.kron_single_gate(num_qubits, QCO.XGate(), "q1") 

    elseif input == "X_2"
        return QCO.kron_single_gate(num_qubits, QCO.XGate(), "q2")    

    elseif input == "Y_1"
        return QCO.kron_single_gate(num_qubits, QCO.YGate(), "q1")  

    elseif input == "Y_2"
        return QCO.kron_single_gate(num_qubits, QCO.YGate(), "q2")   

    elseif input == "Z_1"
        return QCO.kron_single_gate(num_qubits, QCO.ZGate(), "q1") 

    elseif input == "Z_2"
        return QCO.kron_single_gate(num_qubits, QCO.ZGate(), "q2")     

    # Gates with continuous angle parameters
    elseif input in ["RX", "RY", "RZ", "U3"] 
        return QCO.kron_single_gate(num_qubits, matrix, qubit_location)
    
    end

    # All 2-qubit full-sized gates
    if num_qubits == 2

        if input == "CNot_12"
            return QCO.CNotGate()

        elseif input == "CNot_21"
            return QCO.CNotRevGate()

        elseif input == "CNotSwap"
            return QCO.CNotGate() * QCO.CNotRevGate()

        elseif input == "H_1⊗H_2"
            return kron(QCO.HGate(), QCO.HGate())   
        
        elseif input == "X_1⊗X_2"
            return kron(QCO.XGate(), QCO.XGate())

        elseif input == "Y_1⊗Y_2"
            return kron(QCO.YGate(), QCO.YGate())

        elseif input == "Z_1⊗Z_2"
            return kron(QCO.ZGate(), QCO.ZGate())
            
        elseif input == "S_1⊗S_2"
            return kron(QCO.SGate(), QCO.SGate())

        elseif input == "SX_1⊗SX_2"
            return kron(QCO.SXGate(), QCO.SXGate())

        elseif input == "T_1⊗T_2"
            return kron(QCO.TGate(), QCO.TGate())
            
        elseif input == "CZ_12"
            return QCO.CZGate()

        elseif input == "CH_12"
            return QCO.CHGate()

        elseif input == "CV_12"
            return QCO.CVGate()

        elseif input == "Swap"
            return QCO.SwapGate()

        elseif input == "iSwap"
            return QCO.iSwapGate()

        elseif input == "M_12"
            return QCO.MGate()

        elseif input == "QFT_12"
            return QCO.QFT2Gate()

        elseif input == "CSX_12"
            return QCO.C2SXGate()
        
        elseif input == "W_12"
            return QCO.WGate()
        
        elseif input == "HCoin"
            return QCO.HCoinGate()
        
        elseif input in ["CRX_12", "CRY_12", "CRZ_12", "CU3_12", "CRX_21", "CRY_21", "CRZ_21", "CU3_21"]
            return matrix
        
        else
            
            Memento.error(_LOGGER, "Specified input elementary gates or the target gate does not exist in the predefined set of gates.")
        end
    end

    # All 3-qubit full-sized gates
    if num_qubits == 3
        
        if input == "H_3"
            return QCO.kron_single_gate(num_qubits, QCO.HGate(), "q3")
    
        elseif input == "T_3"
            return QCO.kron_single_gate(num_qubits, QCO.TGate(), "q3")
    
        elseif input == "Tdagger_3"
            return QCO.kron_single_gate(num_qubits, QCO.TdaggerGate(), "q3")   
    
        elseif input == "S_3"
            return QCO.kron_single_gate(num_qubits, QCO.SGate(), "q3") 
    
        elseif input == "Sdagger_3"
            return QCO.kron_single_gate(num_qubits, QCO.SdaggerGate(), "q3") 
    
        elseif input == "SX_3"
            return QCO.kron_single_gate(num_qubits, QCO.SXGate(), "q3") 
    
        elseif input == "SXdagger_3"
            return QCO.kron_single_gate(num_qubits, QCO.SXdaggerGate(), "q3")
    
        elseif input == "X_3"
            return QCO.kron_single_gate(num_qubits, QCO.XGate(), "q3") 
    
        elseif input == "Y_3"
            return QCO.kron_single_gate(num_qubits, QCO.YGate(), "q3")  
     
        elseif input == "Z_3"
            return QCO.kron_single_gate(num_qubits, QCO.ZGate(), "q3") 
     
        # Gates with continuous angle parameters
        elseif input in ["RX", "RY", "RZ", "U3"] 
            return QCO.kron_single_gate(num_qubits, matrix, qubit_location)

        elseif input in ["CRX_12", "CRY_12", "CRZ_12", "CU3_12", "CRX_21", "CRY_21", "CRZ_21", "CU3_21"] 
            return kron(matrix, QCO.IGate(1))

        elseif input in ["CRX_23", "CRY_23", "CRZ_23", "CU3_23", "CRX_32", "CRY_32", "CRZ_32", "CU3_32"] 
            return kron(QCO.IGate(1), matrix)

        elseif input in ["CRX_13", "CRY_13", "CRZ_13", "CU3_13"] 
            if startswith(input, "CRX")
                target_gate = QCO.RXGate(target_angle)
            elseif startswith(input, "CRY")
                target_gate = QCO.RYGate(target_angle)
            elseif startswith(input, "CRZ")
                target_gate = QCO.RZGate(target_angle)
            elseif startswith(input, "CU3")
                target_gate = QCO.U3Gate(target_angle[1], target_angle[2], target_angle[3])
            end
            # |0⟩⟨0| ⊗ I ⊗ I 
            control_0 = kron(Array{Complex{Float64},2}([1 0; 0 0]) , kron(QCO.IGate(1), QCO.IGate(1)))
            # |1⟩⟨1| ⊗ I ⊗ R
            control_1 = kron(Array{Complex{Float64},2}([0 0; 0 1]) , kron(QCO.IGate(1), target_gate))
            return control_0 + control_1

        elseif input in ["CRX_31", "CRY_31", "CRZ_31", "CU3_31"] 
            if startswith(input, "CRX")
                target_gate = QCO.RXGate(target_angle)
            elseif startswith(input, "CRY")
                target_gate = QCO.RYGate(target_angle)
            elseif startswith(input, "CRZ")
                target_gate = QCO.RZGate(target_angle)
            elseif startswith(input, "CU3")
                target_gate = QCO.U3Gate(target_angle[1], target_angle[2], target_angle[3])
            end
            # I ⊗ I ⊗ |0⟩⟨0| 
            control_0 = kron(QCO.IGate(1), kron(QCO.IGate(1), Array{Complex{Float64},2}([1 0; 0 0])))
            # R ⊗ I ⊗ |1⟩⟨1|
            control_1 = kron(kron(target_gate, QCO.IGate(1)), Array{Complex{Float64},2}([0 0; 0 1]))
            return control_0 + control_1
        
        elseif input == "Toffoli"
            return QCO.ToffoliGate()

        elseif input == "CSwap"
            return QCO.CSwapGate()

        elseif input == "CCZ"
            return QCO.CCZGate()

        elseif input == "Peres"
            return QCO.PeresGate()

        elseif input == "CNot_12"
            return kron(QCO.CNotGate(), QCO.IGate(1))

        elseif input == "CNot_23"
            return kron(QCO.IGate(1), QCO.CNotGate())

        elseif input == "CNot_21"
            return kron(QCO.CNotRevGate(), QCO.IGate(1))

        elseif input == "CNot_32"
            return kron(QCO.IGate(1), QCO.CNotRevGate())

        elseif input == "CNot_13"
            # |0⟩⟨0| ⊗ I ⊗ I 
            control_0 = kron(Array{Complex{Float64},2}([1 0; 0 0]) , kron(QCO.IGate(1), QCO.IGate(1)))
            # |1⟩⟨1| ⊗ I ⊗ X 
            control_1 = kron(Array{Complex{Float64},2}([0 0; 0 1]) , kron(QCO.IGate(1), QCO.XGate()))
            return control_0 + control_1

        elseif input == "CNot_31"
            # I ⊗ I ⊗ |0⟩⟨0| 
            control_0 = kron(QCO.IGate(1), kron(QCO.IGate(1), Array{Complex{Float64},2}([1 0; 0 0])))
            # X ⊗ I ⊗ |1⟩⟨1|
            control_1 = kron(kron(QCO.XGate(), QCO.IGate(1)), Array{Complex{Float64},2}([0 0; 0 1]))
            return control_0 + control_1

        elseif input == "CV_12"
            return kron(QCO.CVGate(), QCO.IGate(1))

        elseif input == "CV_23"
            return kron(QCO.IGate(1), QCO.CVGate())

        elseif input == "CV_21"
            return kron(QCO.CVRevGate(), QCO.IGate(1))

        elseif input == "CV_32"
            return kron(QCO.IGate(1), QCO.CVRevGate())

        elseif input == "CV_13"
            # |0⟩⟨0| ⊗ I ⊗ I 
            control_0 = kron(Array{Complex{Float64},2}([1 0; 0 0]) , kron(QCO.IGate(1), QCO.IGate(1)))
            # |1⟩⟨1| ⊗ I ⊗ V 
            control_1 = kron(Array{Complex{Float64},2}([0 0; 0 1]) , kron(QCO.IGate(1), QCO.SXGate()))
            return control_0 + control_1

        elseif input == "CV_31"
            # I ⊗ I ⊗ |0⟩⟨0| 
            control_0 = kron(QCO.IGate(1), kron(QCO.IGate(1), Array{Complex{Float64},2}([1 0; 0 0])))
            # V ⊗ I ⊗ |1⟩⟨1|
            control_1 = kron(kron(QCO.SXGate(), QCO.IGate(1)), Array{Complex{Float64},2}([0 0; 0 1]))
            return control_0 + control_1

        elseif input == "CVdagger_12"
            return kron(QCO.CVdaggerGate(), QCO.IGate(1))

        elseif input == "CVdagger_23"
            return kron(QCO.IGate(1), QCO.CVdaggerGate())

        elseif input == "CVdagger_21"
            return kron(QCO.CVRevdaggerGate(), QCO.IGate(1))

        elseif input == "CVdagger_32"
            return kron(QCO.IGate(1), QCO.CVRevdaggerGate())

        elseif input == "CVdagger_13"
            # |0⟩⟨0| ⊗ I ⊗ I 
            control_0 = kron(Array{Complex{Float64},2}([1 0; 0 0]) , kron(QCO.IGate(1), QCO.IGate(1)))
            # |1⟩⟨1| ⊗ I ⊗ V 
            control_1 = kron(Array{Complex{Float64},2}([0 0; 0 1]) , kron(QCO.IGate(1), QCO.SXdaggerGate()))
            return control_0 + control_1

        elseif input == "CVdagger_31"
            # I ⊗ I ⊗ |0⟩⟨0| 
            control_0 = kron(QCO.IGate(1), kron(QCO.IGate(1), Array{Complex{Float64},2}([1 0; 0 0])))
            # V ⊗ I ⊗ |1⟩⟨1|
            control_1 = kron(kron(QCO.SXdaggerGate(), QCO.IGate(1)), Array{Complex{Float64},2}([0 0; 0 1]))
            return control_0 + control_1

        else
            Memento.error(_LOGGER, "Specified input elementary gates or the target gate does not exist in the predefined set of gates.")
        end
    
    end

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
